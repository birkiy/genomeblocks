#!/usr/bin/env Rscript
# The R competitors of the landscape benchmark, one task per call.
#
#   Rscript worker.R <task> <outdir> <args...>
#
# Prints one JSON object per line per (engine, op): {"engine", "op", "seconds":
# [...], "n_out", "out": <file>, ...}. Timings cover the operation only (inputs
# are read once, outside the clock; a warm-up run precedes the timed runs),
# except where the operation *is* reading (parse_gtf, bigwig matrices). Each
# engine writes its answer to <outdir>/<file> in a normalised form that
# bench_landscape.py compares with genomeblocks' answer. An engine that fails
# emits {"error": ...} instead of killing the task.
suppressPackageStartupMessages({library(jsonlite); library(GenomicRanges); library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
task <- args[1]; outdir <- args[2]; rest <- args[-(1:2)]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
REPEAT <- as.integer(Sys.getenv("GB_REPEAT", "3"))
SLOW <- as.integer(Sys.getenv("GB_REPEAT_SLOW", "1"))

emit <- function(engine, op, secs, n_out = NA, out = NULL, ...) {
  cat(toJSON(c(list(engine = engine, op = op, seconds = secs, n_out = n_out, out = out), list(...)),
             auto_unbox = TRUE, null = "null", na = "null", digits = NA), "\n", sep = "")
}
bench <- function(fn, times = REPEAT) {
  res <- fn()                                   # warm-up (and the answer)
  secs <- numeric(times)
  for (i in seq_len(times)) {
    gc(); t0 <- proc.time()[["elapsed"]]; res <- fn(); secs[i] <- proc.time()[["elapsed"]] - t0
  }
  list(secs = secs, res = res)
}
run <- function(engine, op, fn, times = REPEAT) {
  tryCatch({
    b <- bench(fn, times)
    out <- b$res
    emit(engine, op, b$secs, n_out = out$n, out = out$file)
  }, error = function(e) emit(engine, op, NULL, error = conditionMessage(e)))
}
read_bed <- function(p) fread(p, header = FALSE, select = 1:3, col.names = c("chrom", "start", "end"))
to_gr <- function(df, keep_idx = TRUE) {
  g <- GRanges(df$chrom, IRanges(df$start + 1L, df$end))
  if (keep_idx) mcols(g)$idx <- seq_len(nrow(df)) - 1L
  g
}
gr_df <- function(gr) data.table(chrom = as.character(seqnames(gr)), start = start(gr) - 1L, end = end(gr))
save <- function(df, name) { p <- file.path(outdir, name); fwrite(df, p, sep = "\t", col.names = FALSE); list(n = nrow(df), file = p) }

# ── intervals: intersect (A rows overlapping any B), merge (A), nearest (B per A row), and the jobs added
#    2026-10-23: overlap_count (B hits per A row), coverage (A bases covered by B), genome_depth (depth track
#    of A), subtract (A minus B), complement (the genome minus A; args: <A> <B> <ops> <chrom.sizes>) ──
if (task == "intervals") {
  Adf <- read_bed(rest[1]); Bdf <- read_bed(rest[2]); ops <- strsplit(rest[3], ",")[[1]]
  A <- to_gr(Adf); B <- to_gr(Bdf, keep_idx = FALSE)        # only A carries idx: join_nearest must not suffix it
  mcols(B)$bstart <- start(B) - 1L; mcols(B)$bend <- end(B)
  # the chromosome sizes (complement needs them): A2 is A with its seqlengths set, outside the clock
  sizes <- if (length(rest) >= 4) fread(rest[4], header = FALSE, col.names = c("chrom", "size")) else NULL
  A2 <- NULL
  if (!is.null(sizes)) { A2 <- granges(A); seqlevels(A2) <- sizes$chrom; seqlengths(A2) <- sizes$size }
  runs_df <- function(g) {                                  # a coverage GRanges -> (chrom, start, end, depth > 0)
    g <- g[mcols(g)$score > 0]
    data.table(chrom = as.character(seqnames(g)), start = start(g) - 1L, end = end(g), depth = mcols(g)$score)
  }
  # GenomicRanges
  if ("intersect" %in% ops) run("GenomicRanges", "intersect", function() save(gr_df(A[overlapsAny(A, B)]), "gr_intersect.tsv"))
  if ("merge" %in% ops) run("GenomicRanges", "merge", function() save(gr_df(reduce(A)), "gr_merge.tsv"))
  if ("nearest" %in% ops) run("GenomicRanges", "nearest", function() {
    h <- distanceToNearest(A, B)
    save(data.table(idx = queryHits(h) - 1L, bstart = start(B)[subjectHits(h)] - 1L, bend = end(B)[subjectHits(h)]), "gr_nearest.tsv")
  })
  # the jobs added 2026-10-23 (no GenomicRanges function gives coverage of A by B, so that one is not run);
  # setdiff / gaps are namespaced because dplyr masks them once valr is loaded
  if ("overlap_count" %in% ops) run("GenomicRanges", "overlap_count", function() save(cbind(gr_df(A), count = countOverlaps(A, B)), "gr_count.tsv"))
  if ("genome_depth" %in% ops) run("GenomicRanges", "genome_depth", function() save(runs_df(as(coverage(A), "GRanges")), "gr_depth.tsv"))
  if ("subtract" %in% ops) run("GenomicRanges", "subtract", function() save(gr_df(GenomicRanges::setdiff(A, B)), "gr_subtract.tsv"))
  if ("complement" %in% ops && !is.null(A2)) run("GenomicRanges", "complement", function() {
    g <- GenomicRanges::gaps(A2)
    save(gr_df(g[strand(g) == "*"]), "gr_complement.tsv")       # gaps() also returns whole chromosomes on + and -
  })
  # plyranges
  if (requireNamespace("plyranges", quietly = TRUE)) {
    suppressPackageStartupMessages(library(plyranges))
    if ("intersect" %in% ops) run("plyranges", "intersect", function() save(gr_df(filter_by_overlaps(A, B)), "ply_intersect.tsv"))
    if ("merge" %in% ops) run("plyranges", "merge", function() save(gr_df(reduce_ranges(A)), "ply_merge.tsv"))
    if ("nearest" %in% ops) run("plyranges", "nearest", function() {
      j <- join_nearest(A, B)
      save(data.table(idx = mcols(j)$idx, bstart = mcols(j)$bstart, bend = mcols(j)$bend), "ply_nearest.tsv")
    })
    if ("overlap_count" %in% ops) run("plyranges", "overlap_count", function() save(cbind(gr_df(A), count = count_overlaps(A, B)), "ply_count.tsv"))
    if ("genome_depth" %in% ops) run("plyranges", "genome_depth", function() save(runs_df(compute_coverage(A)), "ply_depth.tsv"))
    if ("subtract" %in% ops) run("plyranges", "subtract", function() save(gr_df(setdiff_ranges(A, B)), "ply_subtract.tsv"))
    if ("complement" %in% ops && !is.null(A2)) run("plyranges", "complement", function() save(gr_df(complement_ranges(A2)), "ply_complement.tsv"))
  }
  # valr
  if (requireNamespace("valr", quietly = TRUE)) {
    suppressPackageStartupMessages({library(valr); library(dplyr)})
    At <- tibble::as_tibble(Adf); At$id <- seq_len(nrow(At)) - 1L; Bt <- tibble::as_tibble(Bdf)
    if ("intersect" %in% ops) run("valr", "intersect", function() {
      r <- bed_intersect(At, Bt) %>% distinct(chrom, start.x, end.x)
      save(data.table(chrom = r$chrom, start = r$start.x, end = r$end.x), "valr_intersect.tsv")
    })
    if ("merge" %in% ops) run("valr", "merge", function() {
      r <- bed_merge(bed_sort(At))
      save(data.table(chrom = r$chrom, start = r$start, end = r$end), "valr_merge.tsv")
    })
    if ("nearest" %in% ops) {
      if (nrow(At) > 500000) {
        # bed_closest on 1 M x 1 M sat for 4.5 h at 0.2% CPU on the cluster: not run above 500k rows
        emit("valr", "nearest", NULL, error = sprintf("not run at n = %d: bed_closest did not finish at 1 M (4.5 h, 0.2%% CPU); timed at 100k", nrow(At)))
      } else run("valr", "nearest", function() {
        r <- bed_closest(bed_sort(At), bed_sort(Bt)) %>% distinct(id.x, .keep_all = TRUE)
        save(data.table(idx = r$id.x, bstart = r$start.y, bend = r$end.y), "valr_nearest.tsv")
      })
    }
    # bed_coverage reports both the overlap count (.ints) and the covered bases (.cov); valr has no
    # genome-wide depth function, so genome_depth is not run
    if ("overlap_count" %in% ops) run("valr", "overlap_count", function() {
      r <- bed_coverage(At, Bt)
      save(data.table(chrom = r$chrom, start = r$start, end = r$end, count = r$.ints), "valr_count.tsv")
    })
    if ("coverage" %in% ops) run("valr", "coverage", function() {
      r <- bed_coverage(At, Bt)
      save(data.table(chrom = r$chrom, start = r$start, end = r$end, bases = r$.cov), "valr_coverage.tsv")
    })
    if ("subtract" %in% ops) run("valr", "subtract", function() {
      r <- bed_subtract(At, Bt)
      save(data.table(chrom = r$chrom, start = r$start, end = r$end), "valr_subtract.tsv")
    })
    if ("complement" %in% ops && !is.null(sizes)) {
      gen <- tibble::tibble(chrom = sizes$chrom, size = sizes$size)
      run("valr", "complement", function() {
        r <- bed_complement(At, gen)
        save(data.table(chrom = r$chrom, start = r$start, end = r$end), "valr_complement.tsv")
      })
    }
  }
  # data.table::foverlaps (intersect only)
  if ("intersect" %in% ops) {
    Ad <- copy(Adf); Bd <- copy(Bdf); Ad[, start := start + 1L]; Bd[, start := start + 1L]   # closed [s+1, e]
    setkey(Bd, chrom, start, end)
    run("data.table::foverlaps", "intersect", function() {
      ov <- foverlaps(Ad, Bd, type = "any", nomatch = NULL, which = TRUE)
      keep <- sort(unique(ov$xid))
      save(Adf[keep], "dt_intersect.tsv")
    })
  }
}

# ── annotation: GTF parsing, region classes, nearest TSS ─────────────────────
if (task == "annotation") {
  pk <- read_bed(rest[1]); gtf <- rest[2]; r <- as.integer(rest[3])
  P <- to_gr(pk); names(P) <- as.character(mcols(P)$idx)
  suppressPackageStartupMessages({library(txdbmaker); library(GenomicFeatures)})
  txdb <- NULL
  run("txdbmaker::makeTxDbFromGFF", "parse_gtf", function() { txdb <<- makeTxDbFromGFF(gtf, format = "gtf"); list(n = length(transcripts(txdb)), file = NULL) }, SLOW)
  if (requireNamespace("rtracklayer", quietly = TRUE))
    run("rtracklayer::import", "parse_gtf", function() list(n = length(rtracklayer::import(gtf)), file = NULL), SLOW)
  tx <- transcripts(txdb)
  gn <- genes(txdb, single.strand.genes.only = TRUE)
  # TSS ± r windows (2r+1 bases, like genomeblocks) at two levels: one TSS per gene (genomeblocks'
  # definition) and one per transcript (ChIPseeker's, HOMER's, UROPA's)
  prom_gene <- promoters(gn, upstream = r, downstream = r + 1L)
  prom_tx <- promoters(tx, upstream = r, downstream = r + 1L)
  # GenomicFeatures, the DIY way: the same sets and precedence as genomeblocks.labels
  ex <- reduce(exons(txdb)); u5 <- reduce(unlist(fiveUTRsByTranscript(txdb))); u3 <- reduce(unlist(threeUTRsByTranscript(txdb)))
  body <- reduce(unlist(transcriptsBy(txdb, "gene")))
  diy <- function(promr, file) {
    lab <- rep("Intergenic", length(P))
    lab[overlapsAny(P, body)] <- "Intronic"
    lab[overlapsAny(P, ex)] <- "Exonic"
    lab[overlapsAny(P, u3)] <- "3UTR"
    lab[overlapsAny(P, u5)] <- "5UTR"
    lab[overlapsAny(P, promr)] <- "Promoter-TSS"
    save(data.table(idx = mcols(P)$idx, label = lab), file)
  }
  run("GenomicFeatures DIY labels (gene-level TSS)", "labels", function() diy(reduce(prom_gene), "gf_labels_gene.tsv"))
  run("GenomicFeatures DIY labels (transcript-level TSS)", "labels", function() diy(reduce(prom_tx), "gf_labels_tx.tsv"))
  run("GenomicRanges distanceToNearest (gene-level TSS windows)", "nearest_tss", function() {
    h <- distanceToNearest(P, prom_gene)
    save(data.table(idx = queryHits(h) - 1L, gap = mcols(h)$distance), "gr_nearest_tss_gene.tsv")
  })
  run("GenomicRanges distanceToNearest (transcript-level TSS windows)", "nearest_tss", function() {
    h <- distanceToNearest(P, prom_tx)
    save(data.table(idx = queryHits(h) - 1L, gap = mcols(h)$distance), "gr_nearest_tss_tx.tsv")
  })
  # ChIPseeker
  if (requireNamespace("ChIPseeker", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ChIPseeker))
    run("ChIPseeker annotatePeak", "labels", function() {
      a <- as.data.frame(annotatePeak(P, tssRegion = c(-r, r), TxDb = txdb, level = "transcript", verbose = FALSE))
      tss <- ifelse(a$geneStrand == 1, a$geneStart, a$geneEnd)   # 1-based TSS of the assigned transcript
      save(data.table(idx = a$idx, label = a$annotation, tss1 = tss), "chipseeker_labels.tsv")
    }, SLOW)
  }
  # ChIPpeakAnno
  if (requireNamespace("ChIPpeakAnno", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ChIPpeakAnno))
    names(tx) <- mcols(tx)$tx_name
    run("ChIPpeakAnno annotatePeakInBatch", "nearest_tss", function() {
      a <- annotatePeakInBatch(P, AnnotationData = tx, output = "nearestLocation", FeatureLocForDistance = "TSS", PeakLocForDistance = "middle")
      tss <- ifelse(as.character(mcols(a)$feature_strand) == "+", mcols(a)$start_position, mcols(a)$end_position)
      save(data.table(idx = as.integer(mcols(a)$peak), tss1 = tss), "chippeakanno_nearest.tsv")
    }, SLOW)
  }
}

# ── signal: bigWig -> (regions x 200 bins) matrices ───────────────────────────
if (task == "signal") {
  pk <- read_bed(rest[1]); bw <- rest[2]; flank <- as.integer(rest[3]); nb <- as.integer(rest[4])
  ctr <- (pk$start + pk$end) %/% 2L                                     # 0-based centre, as genomeblocks
  win <- GRanges(pk$chrom, IRanges(ctr - flank + 1L, ctr + flank))      # [c - flank, c + flank) in 0-based
  w <- (2L * flank) %/% nb
  suppressPackageStartupMessages(library(rtracklayer))
  run("rtracklayer::import (bigWig, which=windows)", "bigwig_read", function() list(n = length(import(BigWigFile(bw), which = win)), file = NULL))
  if (requireNamespace("EnrichedHeatmap", quietly = TRUE)) {
    suppressPackageStartupMessages(library(EnrichedHeatmap))
    tgt <- GRanges(pk$chrom, IRanges(ctr + 1L, ctr + 1L))
    run("rtracklayer + EnrichedHeatmap::normalizeToMatrix", "region_matrix", function() {
      sig <- import(BigWigFile(bw), which = win)
      m <- normalizeToMatrix(sig, tgt, value_column = "score", extend = flank, w = w, mean_mode = "weighted",
                             background = 0, smooth = FALSE)
      p <- file.path(outdir, "enrichedheatmap_matrix.tsv"); fwrite(as.data.table(unclass(m)), p, sep = "\t", col.names = FALSE)
      list(n = nrow(m), file = p)
    })
  }
  if (requireNamespace("genomation", quietly = TRUE)) {
    suppressPackageStartupMessages(library(genomation))
    run("genomation::ScoreMatrixBin", "region_matrix", function() {
      sm <- ScoreMatrixBin(target = bw, windows = win, bin.num = nb, bin.op = "mean", type = "bigWig", strand.aware = FALSE)
      p <- file.path(outdir, "genomation_matrix.tsv"); fwrite(as.data.table(sm@.Data), p, sep = "\t", col.names = FALSE)
      list(n = nrow(sm), file = p)
    })
  }
}

# ── motifs: the same log-odds matrices through every R scanner ────────────────
if (task == "motifs") {
  fa <- rest[1]; pssm <- rest[2]; thr <- as.numeric(rest[3]); pval <- as.numeric(rest[4])
  suppressPackageStartupMessages(library(Biostrings))
  seqs <- readDNAStringSet(fa)
  # pssm file: blocks of ">name" then 4 rows A C G T of w numbers (log2-odds)
  lines <- readLines(pssm); starts <- grep("^>", lines)
  mats <- lapply(seq_along(starts), function(i) {
    m <- do.call(rbind, lapply(lines[starts[i] + 1:4], function(l) as.numeric(strsplit(trimws(l), "\\s+")[[1]])))
    rownames(m) <- c("A", "C", "G", "T"); m
  })
  names(mats) <- sub("^>", "", lines[starts])
  counts_file <- function(cnt, name) { p <- file.path(outdir, name); fwrite(data.table(motif = names(cnt), hits = as.integer(cnt)), p, sep = "\t", col.names = FALSE); list(n = sum(cnt), file = p) }
  if (requireNamespace("universalmotif", quietly = TRUE)) {
    suppressPackageStartupMessages(library(universalmotif))
    um <- lapply(names(mats), function(n) create_motif(mats[[n]], type = "PWM", name = n))
    run("universalmotif::scan_sequences", "scan_logodds", function() {
      h <- scan_sequences(um, seqs, threshold = thr, threshold.type = "logodds.abs", RC = FALSE, nthreads = 1, verbose = 0)
      cnt <- table(factor(h$motif, levels = names(mats)))
      counts_file(cnt, "universalmotif_hits.tsv")
    }, SLOW)
  }
  run("Biostrings::countPWM", "scan_logodds", function() {
    cnt <- vapply(mats, function(m) sum(vapply(seq_along(seqs), function(i) countPWM(m, seqs[[i]], min.score = thr), 1L)), 1L)
    counts_file(cnt, "biostrings_hits.tsv")
  }, SLOW)
  if (requireNamespace("motifmatchr", quietly = TRUE)) {
    suppressPackageStartupMessages({library(motifmatchr); library(TFBSTools)})
    pwml <- do.call(PWMatrixList, lapply(names(mats), function(n) PWMatrix(ID = n, name = n, profileMatrix = mats[[n]])))
    run("motifmatchr::matchMotifs (p-value cutoff)", "scan_pvalue", function() {
      x <- matchMotifs(pwml, seqs, bg = "even", out = "scores", p.cutoff = pval)
      cnt <- Matrix::colSums(motifCounts(x)); names(cnt) <- names(mats)
      counts_file(cnt, "motifmatchr_hits.tsv")
    })
  }
}

# ── cooccupancy: LOLA over the 500 tracks (1 kb-tile universe), regioneR on a few ──
if (task == "cooccupancy") {
  query <- rest[1]; tracks_dir <- rest[2]; cs <- rest[3]; n_regioner <- as.integer(rest[4]); nperm <- as.integer(rest[5])
  sizes <- fread(cs, header = FALSE, col.names = c("chrom", "len"))
  Q <- to_gr(read_bed(query))
  tracks <- sort(list.files(tracks_dir, pattern = "\\.bed$", full.names = TRUE))
  if (requireNamespace("LOLA", quietly = TRUE)) {
    suppressPackageStartupMessages(library(LOLA))
    db <- file.path(outdir, "lola_db", "landscape", "regions"); dir.create(db, showWarnings = FALSE, recursive = TRUE)
    for (t in tracks) if (!file.exists(file.path(db, basename(t)))) file.symlink(normalizePath(t), file.path(db, basename(t)))
    sl <- setNames(sizes$len, sizes$chrom)
    universe <- tileGenome(sl, tilewidth = 1000, cut.last.tile.in.chrom = TRUE)
    regionDB <- NULL
    run("LOLA::loadRegionDB (500 tracks)", "index_collection", function() {
      regionDB <<- loadRegionDB(file.path(outdir, "lola_db"), useCache = FALSE, collections = "landscape")
      list(n = length(regionDB$regionGRL), file = NULL)
    }, SLOW)
    run("LOLA::runLOLA (1 kb-tile universe)", "enrichment_search", function() {
      res <- runLOLA(Q, universe, regionDB, cores = 1)
      p <- file.path(outdir, "lola_scores.tsv")
      fwrite(res[, .(name = sub("\\.bed$", "", filename), support, pValueLog, oddsRatio)], p, sep = "\t", col.names = FALSE)
      list(n = nrow(res), file = p)
    }, SLOW)
  }
  if (requireNamespace("regioneR", quietly = TRUE) && n_regioner > 0) {
    suppressPackageStartupMessages(library(regioneR))
    genome <- data.frame(chr = sizes$chrom, start = 1, end = sizes$len)
    run(sprintf("regioneR::overlapPermTest (%d perms, per track)", nperm), "permutation_null", function() {
      rows <- lapply(tracks[seq_len(n_regioner)], function(t) {
        B <- to_gr(read_bed(t), keep_idx = FALSE)
        pt <- overlapPermTest(A = Q, B = B, ntimes = nperm, genome = genome, verbose = FALSE, count.once = TRUE)
        data.table(name = sub("\\.bed$", "", basename(t)), observed = pt$numOverlaps$observed, z = pt$numOverlaps$zscore, p = pt$numOverlaps$pval)
      })
      p <- file.path(outdir, "regioner_scores.tsv"); fwrite(rbindlist(rows), p, sep = "\t", col.names = FALSE)
      list(n = n_regioner, file = p)
    }, SLOW)
  }
}

# ── graph: BEDPE loops vs regions, and the graph algorithms in igraph ─────────
if (task == "graph") {
  bedpe <- rest[1]; cre <- rest[2]; edges_file <- rest[3]; n_vertices <- as.integer(rest[4])
  C <- to_gr(read_bed(cre), keep_idx = FALSE)
  if (requireNamespace("GenomicInteractions", quietly = TRUE)) {
    suppressPackageStartupMessages({library(GenomicInteractions); library(InteractionSet)})
    gi <- NULL
    run("GenomicInteractions::makeGenomicInteractionsFromFile", "read_bedpe", function() { gi <<- makeGenomicInteractionsFromFile(bedpe, type = "bedpe"); list(n = length(gi), file = NULL) })
    run("InteractionSet::overlapsAny (either anchor)", "loop_overlap_either", function() {
      m <- overlapsAny(gi, C, use.region = "both"); save(data.table(idx = which(m) - 1L), "is_either.tsv")
    })
    run("InteractionSet::overlapsAny (both anchors)", "loop_overlap_both", function() {
      m <- overlapsAny(gi, C, use.region = "first") & overlapsAny(gi, C, use.region = "second"); save(data.table(idx = which(m) - 1L), "is_both.tsv")
    })
  }
  if (requireNamespace("igraph", quietly = TRUE)) {
    suppressPackageStartupMessages(library(igraph))
    ed <- fread(edges_file, header = FALSE, col.names = c("src", "tgt", "w"))
    g <- NULL
    run("igraph (R)", "build_graph", function() {
      g <<- graph_from_data_frame(data.frame(from = ed$src, to = ed$tgt, w = ed$w), directed = FALSE, vertices = data.frame(name = 0:(n_vertices - 1L)))
      list(n = ecount(g), file = NULL)
    })
    run("igraph (R)", "components", function() save(data.table(member = components(g)$membership), "igraph_components.tsv"))
    run("igraph (R)", "pagerank", function() save(data.table(pr = page_rank(g, weights = E(g)$w, damping = 0.85)$vector), "igraph_pagerank.tsv"))
  }
}

# ── expression_de: the DE family on simulated counts with a known truth (args: <counts.tsv> <coldata.tsv>;
#    design ~ batch + group, B vs A). Each engine's answer: gene, p, log2 FC (de_test); gene, dispersion;
#    gene, shrunk log2 FC; a gene x sample matrix (vst; batch-corrected counts or log2 CPM) ──
if (task == "expression_de") {
  suppressPackageStartupMessages({library(DESeq2); library(edgeR); library(limma)})
  cts <- as.matrix(read.delim(rest[1], row.names = 1, check.names = FALSE)); mode(cts) <- "integer"
  cd <- read.delim(rest[2], row.names = 1, colClasses = "character")
  cd$group <- relevel(factor(cd$group), "A"); cd$batch <- factor(cd$batch)
  design <- model.matrix(~ batch + group, cd)
  savemat <- function(m, name) {
    p <- file.path(outdir, name)
    fwrite(cbind(data.table(gene = rownames(m)), as.data.table(m)), p, sep = "\t", col.names = FALSE)
    list(n = nrow(m), file = p)
  }
  dds0 <- DESeqDataSetFromMatrix(cts, cd, ~ batch + group)
  dds <- NULL
  run("DESeq2", "de_test", function() {
    dds <<- DESeq(dds0, quiet = TRUE)
    r <- results(dds, name = "group_B_vs_A")
    save(data.table(gene = rownames(r), p = r$pvalue, lfc = r$log2FoldChange), "deseq2_de.tsv")
  }, SLOW)
  run("DESeq2", "dispersion", function() {
    d <- estimateDispersions(estimateSizeFactors(dds0), quiet = TRUE)
    save(data.table(gene = rownames(d), disp = dispersions(d)), "deseq2_disp.tsv")
  }, SLOW)
  if (!is.null(dds)) for (ty in c("apeglm", "ashr", "normal")) if (ty == "normal" || requireNamespace(ty, quietly = TRUE))
    run(sprintf("DESeq2 lfcShrink (%s)", ty), "lfc_shrinkage", function() {
      r <- suppressMessages(lfcShrink(dds, coef = "group_B_vs_A", type = ty, quiet = TRUE))
      save(data.table(gene = rownames(r), lfc = r$log2FoldChange), sprintf("deseq2_shrink_%s.tsv", ty))
    }, SLOW)
  run("DESeq2 vst", "vst", function() savemat(assay(vst(dds0, blind = TRUE)), "deseq2_vst.tsv"), SLOW)
  run("edgeR QL", "de_test", function() {
    y <- estimateDisp(calcNormFactors(DGEList(cts)), design)
    tt <- glmQLFTest(glmQLFit(y, design), coef = "groupB")$table
    save(data.table(gene = rownames(tt), p = tt$PValue, lfc = tt$logFC), "edger_de.tsv")
  }, SLOW)
  run("edgeR estimateDisp", "dispersion", function() {
    y <- estimateDisp(calcNormFactors(DGEList(cts)), design)
    save(data.table(gene = rownames(y), disp = y$tagwise.dispersion), "edger_disp.tsv")
  }, SLOW)
  run("limma-voom", "de_test", function() {
    v <- voom(calcNormFactors(DGEList(cts)), design)
    tt <- topTable(eBayes(lmFit(v, design)), coef = "groupB", number = Inf, sort.by = "none")
    save(data.table(gene = rownames(tt), p = tt$P.Value, lfc = tt$logFC), "limma_de.tsv")
  }, SLOW)
  if (requireNamespace("sva", quietly = TRUE))
    run("sva ComBat_seq", "batch_correction", function()
      savemat(suppressMessages(sva::ComBat_seq(cts, batch = cd$batch, group = cd$group)), "sva_combatseq_counts.tsv"), SLOW)
  run("limma removeBatchEffect", "batch_correction", function() {
    lc <- cpm(calcNormFactors(DGEList(cts)), log = TRUE, prior.count = 1)
    savemat(removeBatchEffect(lc, batch = cd$batch, design = model.matrix(~ group, cd)), "limma_rbe_log2cpm.tsv")
  }, SLOW)
}

# ── bam_counts: reads overlapping each region in each BAM (args: <regions.saf> <bam>...); answer: chrom, start
#    (0-based), end, one count column per BAM ──
if (task == "bam_counts") {
  saf <- read.delim(rest[1]); bams <- rest[-1]
  key <- data.table(chrom = saf$Chr, start = saf$Start - 1L, end = saf$End)
  if (requireNamespace("Rsubread", quietly = TRUE))
    run("Rsubread featureCounts (SAF, allowMultiOverlap)", "multi_bam_counts", function() {
      fc <- suppressMessages(Rsubread::featureCounts(bams, annot.ext = saf, allowMultiOverlap = TRUE, nthreads = 1, verbose = FALSE))
      save(cbind(key, as.data.table(fc$counts)), "rsubread_counts.tsv")
    })
  suppressPackageStartupMessages({library(GenomicAlignments); library(Rsamtools)})
  feats <- GRanges(saf$Chr, IRanges(saf$Start, saf$End))
  run("GenomicAlignments summarizeOverlaps (Union, inter.feature = FALSE)", "multi_bam_counts", function() {
    se <- summarizeOverlaps(feats, BamFileList(bams), mode = "Union", inter.feature = FALSE, ignore.strand = TRUE, singleEnd = TRUE)
    save(cbind(key, as.data.table(assay(se))), "summarizeoverlaps_counts.tsv")
  })
}

# ── count_bam: reads per gene from one RNA-seq BAM (args: <bam> <gtf>); answer: gene, count ──
if (task == "count_bam") {
  bam <- rest[1]; gtf <- rest[2]
  if (requireNamespace("Rsubread", quietly = TRUE))
    run("Rsubread featureCounts (GTF)", "read_counts_bam", function() {
      fc <- suppressMessages(Rsubread::featureCounts(bam, annot.ext = gtf, isGTFAnnotationFile = TRUE, GTF.featureType = "exon",
                                                     GTF.attrType = "gene_id", nthreads = 1, verbose = FALSE))
      save(data.table(gene = rownames(fc$counts), count = fc$counts[, 1]), "rsubread_gene_counts.tsv")
    }, SLOW)
  suppressPackageStartupMessages({library(GenomicAlignments); library(Rsamtools); library(txdbmaker)})
  ebg <- exonsBy(suppressMessages(makeTxDbFromGFF(gtf, format = "gtf")), by = "gene")          # outside the clock
  run("GenomicAlignments summarizeOverlaps (Union)", "read_counts_bam", function() {
    se <- summarizeOverlaps(ebg, BamFile(bam, yieldSize = 2000000L), mode = "Union", ignore.strand = TRUE, singleEnd = TRUE)
    save(data.table(gene = rownames(se), count = assay(se)[, 1]), "summarizeoverlaps_gene_counts.tsv")
  }, SLOW)
}

# ── track_tmm: TMM-normalised counts per million of a region x sample count table (args: <counts.tsv>) ──
if (task == "track_tmm") {
  suppressPackageStartupMessages(library(edgeR))
  m <- as.matrix(read.delim(rest[1]))
  savem <- function(x, name) { p <- file.path(outdir, name); fwrite(as.data.table(x), p, sep = "\t", col.names = FALSE); list(n = nrow(x), file = p) }
  run("edgeR calcNormFactors + cpm", "track_tmm", function() savem(cpm(calcNormFactors(DGEList(m))), "edger_tmm_cpm.tsv"))
  if (requireNamespace("csaw", quietly = TRUE)) {
    suppressPackageStartupMessages(library(SummarizedExperiment))
    se <- SummarizedExperiment(list(counts = m), colData = S4Vectors::DataFrame(totals = colSums(m)))
    run("csaw normFactors (TMM)", "track_tmm", function() {
      f <- csaw::normFactors(se, se.out = FALSE)
      savem(t(t(m) / (colSums(m) * f)) * 1e6, "csaw_tmm_cpm.tsv")
    })
  }
}

# ── write_bigwig: a bigWig from a sorted bedGraph (args: <bedGraph> <chrom.sizes>); read outside the clock ──
if (task == "write_bigwig") {
  suppressPackageStartupMessages(library(rtracklayer))
  bg <- fread(rest[1], header = FALSE, col.names = c("chrom", "start", "end", "score"))
  cs <- fread(rest[2], header = FALSE, col.names = c("chrom", "size"))
  gr <- GRanges(bg$chrom, IRanges(bg$start + 1L, bg$end), score = bg$score, seqinfo = Seqinfo(cs$chrom, cs$size))
  p <- file.path(outdir, "rtracklayer.bw")
  run("rtracklayer export.bw", "write_bigwig", function() { export.bw(gr, p); list(n = length(gr), file = p) }, SLOW)
}

# ── genes_window: (peak, gene) pairs with the gene's TSS base in [start - w, end + w) (args: <peaks.bed>
#    <tss.bed (1 bp, name = gene)> <w>); GRanges' maxgap counts the bases between, so w - 1 gives that rule ──
if (task == "genes_window") {
  pk <- read_bed(rest[1]); w <- as.integer(rest[3])
  ts <- fread(rest[2], header = FALSE, col.names = c("chrom", "start", "end", "gene"))
  P <- to_gr(pk, keep_idx = FALSE); Tg <- GRanges(ts$chrom, IRanges(ts$start + 1L, ts$end), gene = ts$gene)
  pairs <- function(qh, sh, name) save(data.table(chrom = pk$chrom[qh], start = pk$start[qh], end = pk$end[qh], gene = ts$gene[sh]), name)
  run("GenomicRanges findOverlaps (maxgap)", "genes_in_window", function() {
    h <- findOverlaps(P, Tg, maxgap = w - 1L); pairs(queryHits(h), subjectHits(h), "gr_window.tsv")
  })
  if (requireNamespace("plyranges", quietly = TRUE)) {
    suppressPackageStartupMessages(library(plyranges))
    mcols(P)$pi <- seq_along(P)
    run("plyranges join_overlap_inner (maxgap)", "genes_in_window", function() {
      j <- join_overlap_inner(P, Tg, maxgap = w - 1L); save(data.table(chrom = pk$chrom[j$pi], start = pk$start[j$pi], end = pk$end[j$pi], gene = j$gene), "ply_window.tsv")
    })
  }
  if (requireNamespace("valr", quietly = TRUE)) {
    suppressPackageStartupMessages({library(valr); library(dplyr)})
    At <- tibble::as_tibble(pk); Bt <- tibble::as_tibble(ts)
    gen <- At %>% group_by(chrom) %>% summarise(size = max(end) + w + 1L) %>% left_join(Bt %>% group_by(chrom) %>% summarise(s2 = max(end) + 1L), by = "chrom") %>%
      mutate(size = pmax(size, s2, na.rm = TRUE)) %>% select(chrom, size)
    run("valr bed_window (both)", "genes_in_window", function() {
      r <- bed_window(At, Bt, gen, both = w)
      save(data.table(chrom = r$chrom, start = r$start.x, end = r$end.x, gene = r$gene.y), "valr_window.tsv")
    })
  }
  if (requireNamespace("ChIPpeakAnno", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ChIPpeakAnno))
    Pn <- P; names(Pn) <- paste0("p", seq_along(Pn)); Tn <- Tg; names(Tn) <- ts$gene
    run("ChIPpeakAnno annotatePeakInBatch (overlapping, maxgap)", "genes_in_window", function() {
      a <- suppressMessages(annotatePeakInBatch(Pn, AnnotationData = Tn, output = "overlapping", maxgap = w - 1L, select = "all"))
      a <- a[!is.na(a$feature)]
      i <- as.integer(sub("^p", "", a$peak))
      save(data.table(chrom = pk$chrom[i], start = pk$start[i], end = pk$end[i], gene = a$feature), "chippeakanno_window.tsv")
    }, SLOW)
  }
}

# ── universe: enrichment against a custom universe (args: <query as universe regions> <universe.bed> <tracks dir>
#    <n tracks for regioneR> <permutations>); answer: track name, score (LOLA -log10 p; regioneR z) ──
if (task == "universe") {
  Qu <- to_gr(read_bed(rest[1]), keep_idx = FALSE); U <- to_gr(read_bed(rest[2]), keep_idx = FALSE)
  tracks <- sort(list.files(rest[3], pattern = "\\.bed$", full.names = TRUE)); n_regioner <- as.integer(rest[4]); nperm <- as.integer(rest[5])
  if (requireNamespace("LOLA", quietly = TRUE)) {
    suppressPackageStartupMessages(library(LOLA))
    db <- file.path(outdir, "lola_db", "landscape", "regions"); dir.create(db, showWarnings = FALSE, recursive = TRUE)
    for (t in tracks) if (!file.exists(file.path(db, basename(t)))) file.symlink(normalizePath(t), file.path(db, basename(t)))
    regionDB <- loadRegionDB(file.path(outdir, "lola_db"), useCache = FALSE, collections = "landscape")      # outside the clock
    run("LOLA::runLOLA (userUniverse)", "universe_enrichment", function() {
      res <- runLOLA(Qu, U, regionDB, cores = 1)
      p <- file.path(outdir, "lola_universe.tsv")
      fwrite(res[, .(name = sub("\\.bed$", "", filename), score = pValueLog)], p, sep = "\t", col.names = FALSE)
      list(n = nrow(res), file = p)
    }, SLOW)
  }
  if (requireNamespace("regioneR", quietly = TRUE) && n_regioner > 0) {
    suppressPackageStartupMessages(library(regioneR))
    run(sprintf("regioneR::permTest (resampleRegions in the universe, %d perms, per track)", nperm), "universe_enrichment", function() {
      rows <- lapply(tracks[seq_len(n_regioner)], function(t) {
        B <- to_gr(read_bed(t), keep_idx = FALSE)
        pt <- permTest(A = Qu, B = B, ntimes = nperm, randomize.function = resampleRegions, universe = U,
                       evaluate.function = numOverlaps, count.once = TRUE, verbose = FALSE)
        data.table(name = sub("\\.bed$", "", basename(t)), score = pt$numOverlaps$zscore)
      })
      p <- file.path(outdir, "regioner_universe.tsv"); fwrite(rbindlist(rows), p, sep = "\t", col.names = FALSE)
      list(n = n_regioner, file = p)
    }, SLOW)
  }
}

# ── pairwise: base-pair Jaccard matrix of k tracks (args: <track.bed> ...); answer: the k x k matrix ──
if (task == "pairwise") {
  tr <- rest
  if (requireNamespace("valr", quietly = TRUE)) {
    suppressPackageStartupMessages(library(valr))
    tb <- lapply(tr, function(p) tibble::as_tibble(read_bed(p)))
    run(sprintf("valr bed_jaccard x %d pairs", length(tr) * (length(tr) - 1) / 2), "pairwise_cooccurrence", function() {
      k <- length(tb); m <- diag(k)
      for (i in seq_len(k - 1)) for (j in (i + 1):k) m[i, j] <- m[j, i] <- bed_jaccard(tb[[i]], tb[[j]])$jaccard
      p <- file.path(outdir, "valr_jaccard.tsv"); fwrite(as.data.table(m), p, sep = "\t", col.names = FALSE); list(n = k, file = p)
    }, SLOW)
  }
}

# ── motif_compare: query motifs against a database (args: <queries.meme> <db.meme>); answer: query, target, score ──
if (task == "motif_compare") {
  if (requireNamespace("universalmotif", quietly = TRUE)) {
    suppressPackageStartupMessages(library(universalmotif))
    q <- read_meme(rest[1]); db <- read_meme(rest[2]); nq <- length(q)
    run("universalmotif compare_motifs (PCC)", "db_comparison", function() {
      r <- as.data.frame(compare_motifs(c(q, db), compare.to = seq_len(nq), method = "PCC", min.overlap = 5, tryRC = TRUE,
                                        max.p = Inf, max.e = Inf, nthreads = 1))
      save(data.table(query = r$subject, target = r$target, score = r$score), "universalmotif_compare.tsv")   # subject = compare.to
    }, SLOW)
  }
}

# ── chromvar: per-sample motif deviations (args: <peaks.fa, one record per peak> <counts.tsv> <motifs.jaspar>
#    <groups.txt>); answer: motif, one deviation z-score per sample ──
if (task == "chromvar") {
  if (requireNamespace("chromVAR", quietly = TRUE) && requireNamespace("motifmatchr", quietly = TRUE)) {
    suppressPackageStartupMessages({library(chromVAR); library(motifmatchr); library(TFBSTools); library(Biostrings); library(SummarizedExperiment)})
    seqs <- readDNAStringSet(rest[1]); m <- as.matrix(read.delim(rest[2], row.names = 1, check.names = FALSE))
    pfm <- readJASPARMatrix(rest[3], matrixClass = "PFM")
    names(pfm) <- vapply(pfm, ID, character(1))                                 # deviations named by matrix ID
    rr <- GRanges(names(seqs), IRanges(1, width(seqs)))
    se <- SummarizedExperiment(list(counts = m[names(seqs), ]), rowRanges = rr)
    rowData(se)$bias <- letterFrequency(seqs, "GC", as.prob = TRUE)[, 1]            # chromVAR's GC-bias column
    run("chromVAR matchMotifs + computeDeviations", "motif_activity", function() {
      ix <- matchMotifs(pfm, seqs, p.cutoff = 5e-5, bg = "even")
      bg <- getBackgroundPeaks(se, niterations = 50)
      z <- deviationScores(computeDeviations(object = se, annotations = ix, background_peaks = bg))
      p <- file.path(outdir, "chromvar_z.tsv"); fwrite(cbind(data.table(motif = rownames(z)), as.data.table(z)), p, sep = "\t", col.names = FALSE)
      list(n = nrow(z), file = p)
    }, SLOW)
  }
}
