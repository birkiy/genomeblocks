"""Transcript/Gene/Genes definitions and GTF parsing helpers."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union, ClassVar, Callable

# local imports delayed where necessary to avoid circular refs
from .locus import Exon, CDS, UTR, Locus

@dataclass
class Transcript(Locus):
    transcript_id: str = ""
    exons: Optional[List[Exon]] = field(default_factory=list)
    cds: Optional[List[CDS]] = field(default_factory=list)
    utr: Optional[List[UTR]] = field(default_factory=list)

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

    def __post_init__(s):
        if s.transcripts is not None: s.transcripts = {k: v for k,v in s.transcripts.items()}
        s.tss = Locus(s.chrom, s.start, s.start+1) if s.strand == '+' else Locus(s.chrom, s.end, s.end-1)

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
def make(cls, filename, gene_name_key='gene_name', gene_type_key='gene_type', chr_map=None, promoter_r=1000):
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
                    if strand == '+': utr_type = "5'" if end <= genes[gene_id].transcripts[t_id].cds[0].start else "3'"
                    elif strand == '-': utr_type = "3'" if end <= genes[gene_id].transcripts[t_id].cds[0].start else "5'"
                    else: utr_type = "5'"  # default fallback
                u = UTR(chrom, start, end, strand, exon_number=e_number, type=utr_type)
                genes[gene_id].transcripts[t_id].add_utr(u)
            else:
                unmapped.append(feature_type)
    unmapped = list(set(unmapped))
    if len(unmapped) > 0: print('[INFO] Unmapped feature types:', ', '.join(unmapped))
    return genes

# Attach helper functions to Genes to preserve original API assignment
Genes.make = classmethod(make)


# UCSC RefSeq (.ucsc) parser.
# Column layout (tab-separated, '#'-prefixed header line):
#   0:bin 1:name(refseqId) 2:chrom 3:strand 4:txStart 5:txEnd
#   6:cdsStart 7:cdsEnd 8:exonCount 9:exonStarts 10:exonEnds
#   11:score 12:name2(geneSymbol) 13:cdsStartStat 14:cdsEndStat 15:exonFrames
# Each row = one RefSeq transcript; multiple rows can share name2 (gene symbol).
def make_ucsc(cls, filename, chr_map=None, promoter_r=1000, keep_alt_contigs=False):
    """Parse a UCSC RefSeq .ucsc table into a Genes object.

    Args:
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
                g.tss = Locus(g.chrom, g.start, g.start+1) if g.strand == '+' \
                        else Locus(g.chrom, g.end, g.end-1)

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
    return genes

Genes.make_ucsc = classmethod(make_ucsc)


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


def get_tss_transcripts(s, gene_type: Optional[Union[str, List[str]]] = None):
    """Return dict transcript_id -> TSS Locus. Useful for alt-promoter-aware
    proximity analysis (TP53, CDKN2A, TCF7L2 etc. have multiple functional TSSs).
    """
    tss = {}
    for g in s.values():
        if gene_type is not None:
            if isinstance(gene_type, str) and g.gene_type != gene_type: continue
            if isinstance(gene_type, list) and g.gene_type not in gene_type: continue
        if g.transcripts is None: continue
        for t_key, t in g.transcripts.items():
            p = t.start if t.strand == '+' else t.end - 1
            tss[t_key] = Locus(t.chrom, p, p + 1, t.strand)
    return tss


def nearest_transcripts(s, loci):
    """Nearest transcript (alt-promoter-aware) per query locus."""
    from .loci import Loci
    tss = s.get_tss_transcripts()
    l_tss = Loci(tss.values()).slop(s._promoter_r)
    t_names = list(tss.keys())
    return loci.nearest(l_tss, o_names=t_names)


def enhancer_to_genes(s, loci, prox: int = 50_000, level: str = 'gene'):
    """Bundle ROSE-style enhancer→gene mapping into one call.

    For each query locus returns overlap / proximal / closest gene (or transcript)
    sets. Unlike ROSE this uses real interval-to-interval distance for `closest`
    (via pyranges.nearest), not enhancer-center → TSS.

    Args:
        loci: Loci to annotate (enhancers, peaks, etc.).
        prox: proximal window around TSS in bp (±prox). Default 50 kb (ROSE).
        level: 'gene' (default) uses gene body + gene TSS; 'transcript' uses
            transcript body + per-transcript TSS for alt-promoter resolution.

    Returns:
        pandas.DataFrame with columns: uid, overlap, proximal, closest.
        `overlap`/`proximal` are comma-separated name lists; `closest` is a
        single name (may be empty if no feature on that chrom).
    """
    import pandas as pd
    import pyranges as pr
    from .loci import Loci

    if level == 'gene':
        bodies = {g.gene_name: g for g in s.values()}
        tss_map = s.get_tss()
    elif level == 'transcript':
        bodies = {}
        for g in s.values():
            if g.transcripts is None: continue
            for t_key, t in g.transcripts.items():
                bodies[t_key] = t
        tss_map = s.get_tss_transcripts()
    else:
        raise ValueError(f"level must be 'gene' or 'transcript', got {level!r}")

    body_names = list(bodies.keys())
    body_loci = Loci(bodies.values())
    tss_names = list(tss_map.keys())
    tss_loci = Loci(tss_map.values())

    q_df = loci.to_frame().rename(columns={'Chr': 'Chromosome'})
    q_df['rid'] = range(len(q_df))
    q_pr = pr.PyRanges(q_df[['Chromosome', 'Start', 'End', 'rid']])

    body_pr = body_loci.to_pyranges(names=body_names)
    tss_pr = tss_loci.to_pyranges(names=tss_names)

    q_slop = q_df.copy()
    q_slop['Start'] = (q_slop['Start'] - prox).clip(lower=0)
    q_slop['End'] = q_slop['End'] + prox
    q_slop_pr = pr.PyRanges(q_slop[['Chromosome', 'Start', 'End', 'rid']])

    ov = q_pr.join(body_pr, suffix='_b').df
    px = q_slop_pr.join(tss_pr, suffix='_t').df

    def _group(df):
        if len(df) == 0: return {}
        return {rid: sub['Name'].tolist() for rid, sub in df.groupby('rid')}
    ov_g, px_g = _group(ov), _group(px)

    near = loci.nearest(tss_loci, o_names=tss_names)
    # pyranges.nearest puts the neighbor's name in Name_b (suffix default '_b')
    neighbor_col = 'Name_b' if 'Name_b' in near.columns else 'Name'
    near_map = dict(zip(near['Name'], near[neighbor_col]))

    rows = []
    for rid, uid in enumerate(q_df['Name']):
        ov_ids = list(dict.fromkeys(ov_g.get(rid, [])))
        px_ids_raw = list(dict.fromkeys(px_g.get(rid, [])))
        ov_set = set(ov_ids)
        px_ids = [x for x in px_ids_raw if x not in ov_set]
        rows.append({
            'uid': uid,
            'overlap': ','.join(ov_ids),
            'proximal': ','.join(px_ids),
            'closest': near_map.get(uid, ''),
        })
    return pd.DataFrame(rows)


Genes.annotations = annotations
Genes.nearest_genes = nearest_genes
Genes.get_tss_transcripts = get_tss_transcripts
Genes.nearest_transcripts = nearest_transcripts
Genes.enhancer_to_genes = enhancer_to_genes
