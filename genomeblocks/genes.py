"""Transcript/Gene/Genes definitions and GTF parsing helpers."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union

# local imports delayed where necessary to avoid circular refs
from .locus import Exon, CDS, UTR, Locus


def _tss_locus(chrom: str, start: int, end: int, strand: str) -> Locus:
    """TSS as a 1-bp Locus, strand-aware ('-' keeps the historical
    ``Locus(chrom, end, end-1)`` orientation)."""
    return Locus(chrom, start, start + 1) if strand == '+' else Locus(chrom, end, end - 1)

@dataclass
class Transcript(Locus):
    transcript_id: str = ""
    exons: Optional[List[Exon]] = field(default_factory=list)
    cds: Optional[List[CDS]] = field(default_factory=list)
    utr: Optional[List[UTR]] = field(default_factory=list)
    # filled by Genes.select_isoforms(); None = never evaluated
    tss_score: Optional[float] = None
    tss_support: Optional[bool] = None

    @property
    def tss(s) -> Locus: return _tss_locus(s.chrom, s.start, s.end, s.strand)

    def __post_init__(s):
        # Loci import delayed to avoid circular import
        from .loci import Loci
        if s.exons is not None: s.exons = Loci(s.exons)
        if s.cds is not None: s.cds = Loci(s.cds)
        if s.utr is not None: s.utr = Loci(s.utr)

    def add_exon(s, e: Exon) -> None:
        if s.exons is not None:
            s.exons.append(e)
    def add_cds(s, c: CDS) -> None:
        if s.cds is not None:
            s.cds.append(c)
    def add_utr(s, u: UTR) -> None:
        if s.utr is not None:
            s.utr.append(u)


@dataclass
class Gene(Locus):
    gene_id: str = ""
    gene_name: Optional[str] = ""
    gene_type: Optional[str] = ""
    transcripts: Optional[Dict[str, Transcript]] = field(default_factory=dict)
    # transcript key chosen by Genes.select_isoforms(); None = not selected
    canonical: Optional[str] = None

    def __post_init__(s):
        if s.transcripts is not None: s.transcripts = {k: v for k,v in s.transcripts.items()}
        s.tss = _tss_locus(s.chrom, s.start, s.end, s.strand)

    @property
    def canonical_transcript(s) -> Optional[Transcript]:
        return (s.transcripts or {}).get(s.canonical) if s.canonical else None

    def set_span(s, start: int, end: int) -> None:
        """Move the gene body to [start, end) and refresh its TSS."""
        s.start, s.end = start, end
        s.tss = _tss_locus(s.chrom, s.start, s.end, s.strand)

    def add_transript(s, t_id, t: Transcript) -> None:
        if s.transcripts is not None:
            s.transcripts[t_id] = t

@dataclass
class Genes(dict):
    """a dictionary of genes (lightweight container)
    """
    filename: Optional[str] = None
    _annot: Optional[dict] = None
    _promoter_r: Optional[int] = 1000
        

    def _build_annot(s):
        # delay Loci import to avoid circular import
        from .loci import Loci
        # Loci functions are used here, slop
        s._annot = {
            'body' : Loci(s.values()),
            'prom' : Loci(s.get_tss().values()).slop(s._promoter_r).sort().merge(), 
            'exon' : Loci(e for g in s.values() for t in g.transcripts.values() for e in t.exons).sort().merge(),
            'utr5' : Loci(u for g in s.values() for t in g.transcripts.values() for u in t.utr if u.type ==  "5'").sort().merge(),
            'utr3' : Loci(u for g in s.values() for t in g.transcripts.values() for u in t.utr if u.type ==  "3'").sort().merge()
        }
        return

    @property
    def annot(s):
        if s._annot is None: s._build_annot()
        return s._annot

    def table(self):
        lines = []
        if len(self) > 0:
            lines.append(f"{'Name':<20} {'Count':<10}")
            lines.append("-" * 30)
            lines.append(f"{'Transcripts':<20} {sum(len(g.transcripts) for g in self.values()):<10}")
            lines.append(f"{'Exons':<20} {sum(len(t.exons) for g in self.values() for t in g.transcripts.values()):<10}")
            lines.append(f"{'CDS':<20} {sum(len(t.cds) for g in self.values() for t in g.transcripts.values()):<10}")
            lines.append(f"{'UTR':<20} {sum(len(t.utr) for g in self.values() for t in g.transcripts.values()):<10}")
            n_canon = sum(1 for g in self.values() if g.canonical)
            if n_canon: lines.append(f"{'Canonical':<20} {n_canon:<10}")
        return "\n".join(lines)

    def get_tss(s, gene_type: Optional[Union[str, List[str]]] = None):
        tss = {}
        for g in s.values():
            if gene_type is None or (isinstance(gene_type, str) and g.gene_type == gene_type) or (isinstance(gene_type, list) and g.gene_type in gene_type):
                tss[g.gene_name] = g.tss
        return tss


# Helper: parse attributes column

def _parse_attributes(attr_str: str) -> Dict[str, str]:
    attrs = {}
    for attr in attr_str.split(";"):
        if not attr.strip(): continue
        key_value = attr.strip().split(" ", 1)
        if len(key_value) == 2:
            key, value = key_value
            attrs[key] = value.strip('"')
    return attrs


# GTF/GFF make function (as classmethod)
def make(cls, filename, gene_name_key='gene_name', gene_type_key='gene_type', chr_map=None, promoter_r=1000,
         cre=None, bw=None, r=None, kw=None):
    """Parse a GTF/GFF file into a Genes object.

    Args:
        filename: GTF/GFF path.
        gene_name_key/gene_type_key: attribute keys to read gene name/type from.
        chr_map: optional {old_chrom: new_chrom} renaming applied while parsing.
        promoter_r: promoter half-window (bp) used by ``annot['prom']``.
        cre: ATAC/DNase peaks — a BED/narrowPeak path, a ``Loci`` / ``Locus``,
            or a list of any of those — used to keep only the isoforms whose
            TSS sits in open chromatin.
        bw: ATAC bigwig path (or list of paths) scored at each TSS.
        r: TSS half-window (bp) for both of the above; defaults to
            ``promoter_r``.
        kw: extra options forwarded to :meth:`Genes.select_isoforms`
            (``agg``, ``min_signal``, ``min_frac``, ``rank``, ``collapse``,
            ``verbose``).

    Passing ``cre`` and/or ``bw`` runs :meth:`Genes.select_isoforms` after
    parsing, so each gene's body and TSS follow its longest *ATAC-supported*
    isoform instead of its longest annotated one. Peaks/bigwigs must use the
    same chromosome names as the parsed genes (i.e. post ``chr_map``).
    """
    from tqdm import tqdm

    genes = cls(filename=filename, _promoter_r=promoter_r)
    unmapped = []

    with open(filename) as f:
        for line in tqdm(f, desc='[INFO] Parsing GTF/GFF file 🧩', mininterval=30):
            if line.startswith("#") or not line.strip(): continue
            fields = line.strip().split("\t")
            if len(fields) == 9:
                chrom, source, feature_type, start, end, score, strand, phase, attributes = fields
            elif len(fields) == 8:
                chrom, source, feature_type, start, end, score, strand, attributes = fields
            else:
                # Skip malformed lines that do not conform to GTF/GFF column counts
                continue

            if chr_map is not None and chrom in chr_map: chrom = chr_map[chrom]
            start, end = int(start), int(end)
            attrs = _parse_attributes(attributes)
            gene_id = attrs['gene_id']

            if feature_type == "gene":
                gene_name = attrs.get(gene_name_key, gene_id)
                gene_type = attrs.get(gene_type_key, None)
                genes[gene_id] = Gene(chrom, start, end, strand, gene_id=gene_id, gene_name=gene_name, gene_type=gene_type)

            elif feature_type == "transcript":
                t_id = attrs['transcript_id']
                t = Transcript(chrom, start, end, strand, t_id)
                if gene_id not in genes:
                    genes[gene_id] = Gene(chrom, start, end, strand, gene_id=gene_id)
                genes[gene_id].add_transript(t_id, t)

            elif feature_type == "exon":
                t_id = attrs['transcript_id']
                e_number = int(attrs.get('exon_number', 0))
                e = Exon(chrom, start, end, strand, exon_number=e_number)
                genes[gene_id].transcripts[t_id].add_exon(e)

            elif feature_type == "CDS":
                t_id = attrs['transcript_id']
                e_number = int(attrs.get('exon_number', 0))
                c = CDS(chrom, start, end, strand, exon_number=e_number)
                genes[gene_id].transcripts[t_id].add_cds(c)

            elif feature_type in ('five_prime_UTR','three_prime_UTR', 'UTR'):
                t_id = attrs['transcript_id']
                e_number = int(attrs.get('exon_number', 0))
                if feature_type == 'five_prime_UTR': utr_type = "5'"
                elif feature_type == 'three_prime_UTR': utr_type = "3'"
                else:
                    # Generic 'UTR': infer 5'/3' from the CDS start. Guard the
                    # empty-CDS case (non-coding transcript, or a UTR line that
                    # precedes its CDS) so parsing never crashes on IndexError.
                    cds = genes[gene_id].transcripts[t_id].cds
                    cds_start = cds[0].start if cds else None
                    if cds_start is None: utr_type = "5'"
                    elif strand == '+': utr_type = "5'" if end <= cds_start else "3'"
                    elif strand == '-': utr_type = "3'" if end <= cds_start else "5'"
                    else: utr_type = "5'"  # default fallback
                u = UTR(chrom, start, end, strand, exon_number=e_number, type=utr_type)
                genes[gene_id].transcripts[t_id].add_utr(u)
            else:
                unmapped.append(feature_type)
    unmapped = list(set(unmapped))
    if len(unmapped) > 0: print('[INFO] Unmapped feature types:', ', '.join(unmapped))
    if cre is not None or bw is not None:
        genes.select_isoforms(cre, bw, r=r, **(kw or {}))
    return genes

# Attach helper functions to Genes to preserve original API assignment
Genes.make = classmethod(make)


# UCSC RefSeq (.ucsc) parser.
# Column layout (tab-separated, '#'-prefixed header line):
#   0:bin 1:name(refseqId) 2:chrom 3:strand 4:txStart 5:txEnd
#   6:cdsStart 7:cdsEnd 8:exonCount 9:exonStarts 10:exonEnds
#   11:score 12:name2(geneSymbol) 13:cdsStartStat 14:cdsEndStat 15:exonFrames
# Each row = one RefSeq transcript; multiple rows can share name2 (gene symbol).
def make_ucsc(cls, filename, chr_map=None, promoter_r=1000, keep_alt_contigs=False,
              cre=None, bw=None, r=None, kw=None):
    """Parse a UCSC RefSeq .ucsc table into a Genes object.

    Args:
        cre/bw/r/kw: open-chromatin isoform selection, exactly as in
            :meth:`Genes.make`.
        keep_alt_contigs: if False (default), rows whose chrom differs from a
            gene's primary chrom are skipped (alt-contig haplotypes like
            chr6_GL000xxx_alt). If True, those rows land in a separate Gene
            keyed `f"{gene_name}__{chrom}"`. Same-chrom transcript-id collisions
            (MHC paralogs etc.) get a `__N` suffix in either mode so no row is
            silently dropped — `Gene.transcripts` is a dict.
    """
    from tqdm import tqdm

    genes = cls(filename=filename, _promoter_r=promoter_r)

    with open(filename) as f:
        for line in tqdm(f, desc='[INFO] Parsing UCSC RefSeq file 🧩', mininterval=30):
            if line.startswith('#') or not line.strip(): continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 13: continue

            t_id       = fields[1]
            chrom      = fields[2]
            strand     = fields[3]
            tx_start   = int(fields[4])
            tx_end     = int(fields[5])
            cds_start  = int(fields[6])
            cds_end    = int(fields[7])
            exon_count = int(fields[8])
            exon_starts = [int(x) for x in fields[9].rstrip(',').split(',') if x]
            exon_ends   = [int(x) for x in fields[10].rstrip(',').split(',') if x]
            gene_name = fields[12] if len(fields) > 12 and fields[12] else t_id

            if chr_map is not None and chrom in chr_map: chrom = chr_map[chrom]

            # One Gene per (gene_name, chrom). The primary chrom keeps the bare
            # name; alt-contig copies (when kept) get a __<chrom> suffix.
            gene_id = gene_name
            if gene_id in genes and chrom != genes[gene_id].chrom:
                if not keep_alt_contigs: continue
                gene_id = f"{gene_name}__{chrom}"

            if gene_id not in genes:
                genes[gene_id] = Gene(chrom, tx_start, tx_end, strand,
                                      gene_id=gene_id, gene_name=gene_name, gene_type=None)
            else:
                g = genes[gene_id]
                if tx_start < g.start: g.start = tx_start
                if tx_end   > g.end:   g.end   = tx_end
                g.tss = _tss_locus(g.chrom, g.start, g.end, g.strand)

            # Disambiguate same-chrom transcript-id collisions (MHC paralogs).
            t_key = t_id
            n = 2
            while t_key in genes[gene_id].transcripts:
                t_key = f"{t_id}__{n}"; n += 1
            t = Transcript(chrom, tx_start, tx_end, strand, t_id)
            genes[gene_id].add_transript(t_key, t)

            coding = cds_end > cds_start
            for i in range(exon_count):
                es, ee = exon_starts[i], exon_ends[i]
                # exon_number is 1-based along the direction of transcription
                e_number = i + 1 if strand == '+' else exon_count - i
                t.add_exon(Exon(chrom, es, ee, strand, exon_number=e_number))

                if not coding: continue
                cs, ce = max(es, cds_start), min(ee, cds_end)
                if cs < ce:
                    t.add_cds(CDS(chrom, cs, ce, strand, exon_number=e_number))
                if es < cds_start:
                    utype = "5'" if strand == '+' else "3'"
                    t.add_utr(UTR(chrom, es, min(ee, cds_start), strand,
                                  exon_number=e_number, type=utype))
                if ee > cds_end:
                    utype = "3'" if strand == '+' else "5'"
                    t.add_utr(UTR(chrom, max(es, cds_end), ee, strand,
                                  exon_number=e_number, type=utype))
    if cre is not None or bw is not None:
        genes.select_isoforms(cre, bw, r=r, **(kw or {}))
    return genes

Genes.make_ucsc = classmethod(make_ucsc)


# ── ATAC-supported isoform selection ─────────────────────────────────────────
# A GTF gene spans the union of its isoforms, so its body and TSS follow the
# longest *annotated* transcript — frequently a long isoform that is silent in
# the cell type at hand (TGFBR3 is a textbook case). Given open-chromatin
# evidence — ATAC/DNase peaks, a bigwig, or both — we can keep the isoforms
# whose TSS is actually accessible and let the gene follow the longest of those.

def _tss_window(t: Transcript, r: int) -> Locus:
    """The ±``r`` bp window around a transcript's TSS (strand-aware)."""
    pos = t.start if t.strand == '+' else t.end
    return Locus(t.chrom, max(0, pos - r), pos + r, t.strand)


def _as_loci(cre):
    """Coerce a BED path / Loci / Locus / iterable of those into one Loci."""
    from .loci import Loci
    if cre is None: return None
    if isinstance(cre, str): return Loci.make(cre)
    if isinstance(cre, Loci): return cre
    if isinstance(cre, Locus): return Loci([cre])
    parts = [_as_loci(p) for p in cre]
    parts = [p for p in parts if p is not None]
    if not parts: return None
    if len(parts) == 1: return parts[0]
    out = Loci()
    for part in parts: out.extend(part)
    return out.sort().merge()


def _as_paths(bw) -> List[str]:
    if bw is None: return []
    if isinstance(bw, str): return [bw]
    return [str(p) for p in bw]


def _rank_key(t: Transcript, rank: str):
    length, score = t.end - t.start, (t.tss_score or 0.0)
    return (length, score) if rank == 'longest' else (score, length)


def select_isoforms(s, cre=None, bw=None, *, r=None, agg='max',
                    min_signal=0.0, min_frac=0.5, rank='longest',
                    collapse=True, verbose=True):
    """Pick one open-chromatin-supported isoform per gene.

    A transcript is *supported* when it passes every criterion given:

    * ``cre`` — its TSS window (±``r`` bp) overlaps an ATAC/DNase peak;
    * ``bw`` — its TSS window scores ``> min_signal`` **and** reaches
      ``min_frac`` of the best TSS score among that gene's own candidate
      isoforms (a within-gene relative cut, so no absolute threshold has to be
      guessed). With several bigwigs a TSS keeps its highest score, i.e. a
      TSS open in *any* of the samples counts as supported.

    The winner among the supported isoforms is the longest one
    (``rank='longest'``), or the one with the strongest TSS signal
    (``rank='signal'``); the other value breaks ties. Genes with no supported
    isoform fall back to the same ranking over *all* their isoforms, so every
    gene keeps a canonical transcript.

    Results land on the objects: ``Transcript.tss_score`` /
    ``Transcript.tss_support`` per isoform, ``Gene.canonical`` (and
    ``Gene.canonical_transcript``) per gene. With ``collapse=True`` (default)
    the gene body and TSS are moved onto the chosen isoform, which is what
    ``annot``, ``get_tss()``, ``annotations()`` and ``nearest_genes()`` read;
    no transcript is ever dropped.

    Args:
        cre: BED/narrowPeak path, ``Loci``, ``Locus``, or a list of any of
            those. A ``Loci`` is used as-is — its interval index is reused.
        bw: bigwig path or list of paths.
        r: TSS half-window in bp; defaults to the object's ``promoter_r``.
        agg: how to summarise the signal across the window
            ('max' | 'mean' | 'sum' | 'min' | 'std' | 'coverage').
        min_signal: absolute floor a TSS score must exceed (default 0.0).
        min_frac: fraction of the gene's best TSS score a supported isoform
            must reach (default 0.5; use 0 to disable the relative cut).
        rank: 'longest' (default) or 'signal'.
        collapse: move each gene's body/TSS onto its canonical isoform.
        verbose: print the support summary.

    Returns:
        self (chainable).
    """
    from .loci import Loci

    if rank not in ('longest', 'signal'):
        raise ValueError(f"rank must be 'longest' or 'signal', got {rank!r}")

    if r is None: r = s._promoter_r
    peak_loci = _as_loci(cre)
    bigwigs = _as_paths(bw)
    if peak_loci is None and not bigwigs:
        raise ValueError("select_isoforms() needs `cre` (BED/narrowPeak/Loci) and/or `bw` (bigwig).")

    tx = [t for g in s.values() for t in (g.transcripts or {}).values()]
    for t in tx: t.tss_score, t.tss_support = None, None

    # 1. peak filter (cheap) — only the survivors are worth a bigwig query
    if peak_loci is not None:
        cand = []
        for t in tx:
            w = _tss_window(t, r)
            if any(True for *_, _ in peak_loci.cgr.overlap(w.chrom, w.start, w.end)):
                cand.append(t)
    else:
        cand = list(tx)

    # 2. score the candidate TSSs; isoforms sharing a TSS share one query
    if bigwigs and cand:
        seen, at = {}, []
        for t in cand:
            w = _tss_window(t, r)
            key = (w.chrom, w.start, w.end)
            if key not in seen: seen[key] = len(seen)
            at.append(seen[key])
        wins = Loci(Locus(c, a, b) for c, a, b in seen)      # dict keeps insertion order
        cube = wins.signal(bigwigs, n_bins=1, span=True, agg=agg,
                           progress=verbose, verbose=verbose)
        scores = cube[:, :, 0].max(axis=1)                   # open in any sample
        for i, t in enumerate(cand): t.tss_score = float(scores[at[i]])

    for t in cand: t.tss_support = True
    for t in tx:
        if t.tss_support is None: t.tss_support = False

    # 3. within-gene signal cut, then pick a canonical isoform per gene
    n_genes = n_supported = n_fallback = 0
    for g in s.values():
        items = list((g.transcripts or {}).items())
        if not items: continue
        n_genes += 1

        if bigwigs:
            open_tx = [t for _, t in items if t.tss_support]
            best_score = max((t.tss_score for t in open_tx), default=0.0)
            for t in open_tx:
                t.tss_support = (t.tss_score > min_signal
                                 and t.tss_score >= min_frac * best_score)

        pool = [(k, t) for k, t in items if t.tss_support]
        if pool:
            n_supported += 1
        else:
            pool = items          # nothing open here: fall back to every isoform
            n_fallback += 1

        pool = sorted(pool, key=lambda kt: kt[0])            # deterministic ties
        pool.sort(key=lambda kt: _rank_key(kt[1], rank), reverse=True)
        k, t = pool[0]
        g.canonical = k
        if collapse: g.set_span(t.start, t.end)

    s._annot = None                                          # spans moved: rebuild lazily
    if verbose:
        n_tx = sum(1 for t in tx if t.tss_support)
        print(f"[INFO] Isoform support: {n_supported}/{n_genes} genes with an open TSS "
              f"({n_tx}/{len(tx)} isoforms), {n_fallback} fell back to the longest isoform.")
        if n_supported == 0:
            print("[WARN] No isoform was supported — check that the peaks/bigwig use the "
                  "same chromosome names as the annotation.")
    return s


Genes.select_isoforms = select_isoforms


def annotations(s: Genes, L):
    import pandas as pd
    annot = s.annot
    labels = []
    for l in L:
        lab = "Intergenic"
        if any(True for _ in annot['prom'].overlaps(l)):
            lab = "Promoter-TSS"
        elif any(True for _ in annot['utr5'].overlaps(l)):
            lab = "5UTR"
        elif any(True for _ in annot['utr3'].overlaps(l)):
            lab = "3UTR"
        elif any(True for _ in annot['exon'].overlaps(l)):
            lab = "Exonic"
        elif any(True for _ in annot['body'].overlaps(l)):
            lab = "Intronic"
        labels.append({"uid": l.uid, "annotation": lab})
    return pd.DataFrame(labels)


def nearest_genes(s, loci):
    from .loci import Loci
    l_tss = Loci(s.get_tss().values()).slop(s._promoter_r)
    g_names = list(s.get_tss().keys())
    return loci.nearest(l_tss, o_names=g_names)


Genes.annotations = annotations
Genes.nearest_genes = nearest_genes
