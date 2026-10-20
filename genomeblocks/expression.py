"""Expression: a features x samples matrix whose rows are the rows of a
:class:`~genomeblocks.Genes` table (or of any :class:`~genomeblocks.Loci`).

    import genomeblocks as gb
    genes = gb.Genes.make("gencode.gtf")
    E = gb.Expression.make("salmon/*/quant.sf", genes=genes)        # tximport -> aligned to genes.genes
    E.counts[i], E.tpm[i]                                            # gene row i, one value per sample
    E.samples["condition"] = ["WT", "WT", "KO", "KO"]
    tss = genes.get_tss(); tss["tpm"] = E.tpm.mean(axis=1)           # no join: the row is the key

Layers (``counts`` always; ``tpm``, ``length`` when the quantifier gives
them) are float64 arrays of shape (n features, m samples). With ``genes=``
the rows are exactly the rows of ``genes.genes`` (``level="gene"``) or
``genes.transcripts`` (``level="transcript"``): an input feature that
matches no row is dropped and counted (``E.unmatched``), a row the input did
not quantify is NaN in every layer (``E.missing``). With ``loci=`` the rows
are the rows of that ``Loci`` (a featureCounts matrix over peaks). Without
either, rows are the input's features in file order and :meth:`align` lines
them up later.

Inputs: salmon ``quant.sf`` / ``quant.genes.sf``, kallisto ``abundance.tsv``,
RSEM ``.genes.results`` / ``.isoforms.results``, featureCounts, STAR
``ReadsPerGene.out.tab``, HTSeq-count, nf-core/rnaseq matrices, any
id + samples matrix (TSV / CSV / parquet), any frame (:meth:`from_frame`), an
AnnData (:meth:`from_anndata`). Formats are sniffed from the header; the
``tables`` backend (polars, else pandas) parses.

Transcript-level inputs are summarised to genes over ``transcripts['gene']``
exactly as ``tximport::summarizeToGene`` (:meth:`to_genes`). Normalisation
(:meth:`cpm`, :meth:`tmm`, :meth:`size_factors`, :meth:`filter_by_expr`) is
pure numpy from :mod:`genomeblocks._norm`. Differential expression is the
``de`` backend family (phase 2).
"""
from __future__ import annotations

import glob as _glob
import gzip
import os
import re
from typing import Dict, List, Optional, Sequence

import numpy as np

from ._table import TableMixin

LEVELS = ("gene", "transcript", "loci")
HOW = ("counts", "scaled_tpm", "length_scaled_tpm")

_SALMON = ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]
_KALLISTO = ["target_id", "length", "eff_length", "est_counts", "tpm"]
_RSEM_GENE = ["gene_id", "transcript_id(s)", "length", "effective_length", "expected_count", "TPM", "FPKM"]
_RSEM_TX = ["transcript_id", "gene_id", "length", "effective_length", "expected_count", "TPM", "FPKM", "IsoPct"]
_STAR_ROWS = ("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
_PER_SAMPLE_FILES = {"quant.sf": "parent", "quant.genes.sf": "parent", "abundance.tsv": "parent"}
_NAME_COLUMNS = ("gene_name", "gene_symbol", "symbol", "genesymbol", "gene_names", "external_gene_name")
_DROP_COLUMNS = ("chr", "start", "end", "strand", "gene_type", "gene_biotype", "biotype", "transcript_type",
                 "gene_id", "transcript_id", "tx", "transcript")
_ID_COLUMNS = ("gene_id", "geneid", "gene", "genes", "id", "ids", "feature", "feature_id", "transcript_id",
               "tx", "transcript", "target_id", "name", "ensembl_gene_id", "ensembl_id", "uid")
_VERSION_RE = re.compile(r"\.\d+(_PAR_Y)?$")
_REGION_RE = re.compile(r"^(\S+?)[:_\-](\d+)[\-_:](\d+)(?:\([+\-.]\))?$")
_BAM_RE = re.compile(r"(\.(bam|sam|cram))$", re.IGNORECASE)
_BAM_TAILS = (".sorted", ".markdup", ".mkdup", ".dedup", ".filtered", ".rmdup", ".mLb.clN", ".mLb.mkD",
              ".Lb.clN", "Aligned.sortedByCoord.out", "Aligned.out", ".sortedByCoord.out")
_FILE_TAILS = (".gz", ".txt", ".tsv", ".csv", ".tab", ".counts", ".count", ".genes.results", ".isoforms.results",
               ".htseq", ".featureCounts", ".parquet")


# ── text tables ─────────────────────────────────────────────────────────────

def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _head_lines(path, n: int = 8) -> List[str]:
    out = []
    with _open_text(path) as f:
        for line in f:
            out.append(line.rstrip("\r\n"))
            if len(out) >= n:
                break
    return out


def _sep_of(path, line: str) -> str:
    if str(path).lower().removesuffix(".gz").endswith(".csv"):
        return ","
    if "\t" in line:
        return "\t"
    return "," if "," in line else "\t"


def _read_delim(path, *, header: bool, skip: int = 0, sep: str = "\t", numeric=(), numeric_after=None,
                backend=None):
    """A delimited text table as ``(column names, {name: numpy array})``:
    ``numeric`` columns (by name), and every column from position
    ``numeric_after`` on, as float64; the rest as object strings."""
    from .backends import resolve
    b = resolve("tables", backend)
    numeric = set(numeric)
    if b == "polars":
        import polars as pl
        try:
            df = pl.read_csv(path, separator=sep, has_header=header, skip_rows=skip, infer_schema_length=0,
                             quote_char=None, truncate_ragged_lines=True)
        except Exception as e:                       # noqa: BLE001 — polars has many error classes
            raise ValueError(f"{path}: cannot read as a delimited table: {e}") from None
        names = [f"c{i}" for i in range(df.width)] if not header else list(df.columns)
        if not header:
            df.columns = names
        if numeric_after is not None:
            numeric |= set(names[numeric_after:])
        out = {}
        for c in names:
            s = df[c]
            out[c] = (s.cast(pl.Float64, strict=False).to_numpy() if c in numeric
                      else s.to_numpy().astype(object))
        return names, out
    import pandas as pd
    try:
        df = pd.read_csv(path, sep=sep, header=0 if header else None, skiprows=skip, dtype=str, quoting=3,
                         engine="c")
    except (ValueError, TypeError) as e:
        raise ValueError(f"{path}: cannot read as a delimited table: {e}") from None
    names = [f"c{i}" for i in range(df.shape[1])] if not header else [str(c) for c in df.columns]
    df.columns = names
    if numeric_after is not None:
        numeric |= set(names[numeric_after:])
    out = {}
    for c in names:
        out[c] = (pd.to_numeric(df[c], errors="coerce").to_numpy(np.float64) if c in numeric
                  else df[c].to_numpy(object))
    return names, out


def _is_number(s) -> bool:
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return str(s).strip().upper() in ("NA", "NAN", "")


# ── the formats ─────────────────────────────────────────────────────────────

class _Table:
    """What one input file holds: ids, optional names, layers (n x k), k sample names."""

    def __init__(self, ids, layers, samples, *, names=None, level=None, fmt="", notes=()):
        self.ids = np.asarray(ids, dtype=object)
        self.names = None if names is None else np.asarray(names, dtype=object)
        self.layers = {k: np.asarray(v, dtype=np.float64) for k, v in layers.items()}
        self.samples = [str(s) for s in samples]
        self.level = level
        self.fmt = fmt
        self.notes = list(notes)
        n, k = len(self.ids), len(self.samples)
        for name, v in self.layers.items():
            if v.shape != (n, k):
                raise ValueError(f"{fmt}: layer {name!r} has shape {v.shape}, expected {(n, k)}")


def sniff_format(path) -> dict:
    """What a quantification file is: ``{'format', 'header', 'skip', 'sep', 'level'}``.

    Recognised by the header line: salmon, kallisto, RSEM, featureCounts,
    STAR ``ReadsPerGene``, HTSeq-count (two columns), nf-core / generic
    matrices (an id column, then samples).
    """
    lines = _head_lines(path)
    if not lines:
        raise ValueError(f"{path}: empty file")
    base = os.path.basename(str(path)).removesuffix(".gz")
    skip = 0
    if lines[0].startswith("# Program:featureCounts"):
        skip = 1
    while skip < len(lines) and (not lines[skip].strip() or (lines[skip].startswith("#") and skip > 0)):
        skip += 1
    if skip >= len(lines):
        raise ValueError(f"{path}: no data lines")
    first = lines[skip]
    sep = _sep_of(path, first)
    fields = [f.strip() for f in first.split(sep)]
    low = [f.lower() for f in fields]
    info = {"format": "matrix", "header": True, "skip": skip, "sep": sep, "level": None}
    if fields[:5] == _SALMON:
        info.update(format="salmon", level="gene" if base == "quant.genes.sf" else "transcript")
    elif fields[:5] == _KALLISTO:
        info.update(format="kallisto", level="transcript")
    elif fields[:7] == _RSEM_GENE:
        info.update(format="rsem", level="gene")
    elif fields[:5] == _RSEM_TX[:5]:
        info.update(format="rsem", level="transcript")
    elif lines[0].startswith("# Program:featureCounts") or low[:6] == ["geneid", "chr", "start", "end", "strand",
                                                                        "length"]:
        info.update(format="featurecounts", level=None)
    elif fields[0] in _STAR_ROWS:
        info.update(format="star", header=False)
    elif len(fields) == 2 and _is_number(fields[1]) and not _is_number(fields[0]):
        info.update(format="htseq", header=False)
    else:
        if len(fields) < 2:
            raise ValueError(f"{path}: expected an id column followed by sample columns, got one column")
        # a header-less matrix: every field but the first is a number
        if all(_is_number(f) for f in fields[1:]) and not _is_number(fields[0]) and len(lines) > skip + 1:
            info["header"] = False
        if low[0] in ("tx", "transcript_id", "transcript", "target_id"):
            info["level"] = "transcript"
        elif low[0] in ("gene_id", "geneid", "gene", "ensembl_gene_id", "gene_name"):
            info["level"] = "gene"
    return info


def _layer_name_from_path(path) -> str:
    low = os.path.basename(str(path)).lower()
    if "tpm" in low:
        return "tpm"
    if "fpkm" in low:
        return "fpkm"
    if "length_scaled" in low:
        return "counts_length_scaled"
    if "counts_scaled" in low or "scaled" in low:
        return "counts_scaled"
    return "counts"


def _bam_stem(name: str) -> str:
    s = os.path.basename(str(name))
    s = _BAM_RE.sub("", s)
    changed = True
    while changed:
        changed = False
        for t in _BAM_TAILS:
            if s.endswith(t) and len(s) > len(t):
                s = s[: -len(t)]
                changed = True
    return s.rstrip("._")


def sample_name(path) -> str:
    """The sample a per-sample file belongs to: the directory for salmon /
    kallisto layouts, else the file name without its format suffix."""
    p = str(path)
    base = os.path.basename(p)
    if base in _PER_SAMPLE_FILES or base.removesuffix(".gz") in _PER_SAMPLE_FILES:
        return os.path.basename(os.path.dirname(os.path.abspath(p))) or base
    if base.endswith("ReadsPerGene.out.tab"):
        return base[: -len("ReadsPerGene.out.tab")].rstrip("._") or base
    s = base
    changed = True
    while changed:
        changed = False
        for t in _FILE_TAILS:
            if s.lower().endswith(t) and len(s) > len(t):
                s = s[: -len(t)]
                changed = True
    return s


def read_table(path, *, strand: str = "auto", sample: Optional[str] = None, backend=None) -> _Table:
    """One quantification file as a :class:`_Table` (internal; :meth:`Expression.make` is the API)."""
    info = sniff_format(path)
    fmt, sep, skip = info["format"], info["sep"], info["skip"]
    name = sample if sample is not None else sample_name(path)
    if fmt == "salmon":
        names, d = _read_delim(path, header=True, skip=skip, sep=sep, numeric=_SALMON[1:], backend=backend)
        return _Table(d["Name"], {"counts": d["NumReads"][:, None], "tpm": d["TPM"][:, None],
                                  "length": d["EffectiveLength"][:, None]}, [name], level=info["level"], fmt=fmt)
    if fmt == "kallisto":
        names, d = _read_delim(path, header=True, skip=skip, sep=sep, numeric=_KALLISTO[1:], backend=backend)
        return _Table(d["target_id"], {"counts": d["est_counts"][:, None], "tpm": d["tpm"][:, None],
                                       "length": d["eff_length"][:, None]}, [name], level="transcript", fmt=fmt)
    if fmt == "rsem":
        numeric = ["length", "effective_length", "expected_count", "TPM", "FPKM"]
        names, d = _read_delim(path, header=True, skip=skip, sep=sep, numeric=numeric, backend=backend)
        idc = "gene_id" if info["level"] == "gene" else "transcript_id"
        return _Table(d[idc], {"counts": d["expected_count"][:, None], "tpm": d["TPM"][:, None],
                               "length": d["effective_length"][:, None]}, [name], level=info["level"], fmt=fmt)
    if fmt == "featurecounts":
        names, d = _read_delim(path, header=True, skip=skip, sep=sep, numeric_after=5, backend=backend)
        samples = names[6:]
        counts = np.column_stack([d[c] for c in samples]) if samples else np.zeros((len(d[names[0]]), 0))
        length = d[names[5]]
        return _Table(d[names[0]], {"counts": counts, "length": np.repeat(length[:, None], len(samples), axis=1)},
                      [_bam_stem(c) for c in samples], level=None, fmt=fmt)
    if fmt == "star":
        names, d = _read_delim(path, header=False, skip=skip, sep=sep, numeric=["c1", "c2", "c3"], backend=backend)
        ids = d["c0"]
        body = ~np.isin(ids.astype(str), _STAR_ROWS)
        tot = {k: np.nansum(d[k][body]) for k in ("c1", "c2", "c3")}
        notes = []
        pick = {"none": "c1", "unstranded": "c1", "forward": "c2", "yes": "c2", "reverse": "c3"}.get(str(strand).lower())
        if pick is None:
            if strand != "auto":
                raise ValueError(f"strand must be 'auto', 'none', 'forward' or 'reverse', got {strand!r}")
            hi, lo = max(tot["c2"], tot["c3"]), min(tot["c2"], tot["c3"])
            if hi > 2 * max(lo, 1):
                pick = "c2" if tot["c2"] >= tot["c3"] else "c3"
                notes.append(f"STAR strand=auto: {'forward' if pick == 'c2' else 'reverse'}-stranded counts "
                             f"(column {3 if pick == 'c2' else 4}; the other strand has far fewer reads)")
            else:
                pick = "c1"
                notes.append("STAR strand=auto: the two stranded columns are alike; using the unstranded counts "
                             "(column 2)")
        return _Table(ids[body], {"counts": d[pick][body][:, None]}, [name], level=None, fmt=fmt, notes=notes)
    if fmt == "htseq":
        names, d = _read_delim(path, header=False, skip=skip, sep=sep, numeric=["c1"], backend=backend)
        ids = d["c0"]
        special = np.char.startswith(ids.astype(str), "__")
        notes = [f"dropped the {int(special.sum())} HTSeq summary rows (__no_feature, __ambiguous, ...)"] \
            if special.any() else []
        return _Table(ids[~special], {"counts": d["c1"][~special][:, None]}, [name], level=None, fmt=fmt, notes=notes)
    # a matrix: id column, optional name column, numeric sample columns
    names, d = _read_delim(path, header=info["header"], skip=skip, sep=sep, backend=backend)
    t = _matrix_table(names, d, path=path, level=info["level"], header=info["header"],
                      layer=_layer_name_from_path(path))
    if not info["header"] and sample is not None and len(t.samples) == 1:
        t.samples = [name]
    return t


def _to_float(a) -> np.ndarray:
    if np.asarray(a).dtype.kind == "f":
        return np.asarray(a, np.float64)
    import pandas as pd
    return pd.to_numeric(pd.Series(np.asarray(a, dtype=object)), errors="coerce").to_numpy(np.float64)


def _matrix_table(names, d, *, path="", level=None, header=True, ids=None, name_col=None, layer=None) -> _Table:
    """A generic id + samples matrix (from a text file or a frame)."""
    cols = list(names)
    if ids is None:
        ids = cols[0]
    elif ids not in cols:
        raise KeyError(f"{path}: id column {ids!r} not found; columns are {cols}")
    rest = [c for c in cols if c != ids]
    low = {c: str(c).lower() for c in cols}
    if name_col is None:
        name_col = next((c for c in rest if low[c] in _NAME_COLUMNS), None)
    elif name_col not in cols:
        raise KeyError(f"{path}: name column {name_col!r} not found; columns are {cols}")
    rest = [c for c in rest if c != name_col]
    numeric, dropped, length = [], [], None
    for c in rest:
        v = d[c]
        if low[c] in ("length", "gene_length", "eff_length", "effectivelength", "effective_length"):
            length = _to_float(v)
            continue
        if low[c] in _DROP_COLUMNS:
            dropped.append(c)
            continue
        f = _to_float(v)
        fin = np.isfinite(f)
        # a sample column: every value a number, or at least half of them (NA / blanks allowed)
        if fin.all() or (len(f) and fin.sum() >= 0.5 * len(f)):
            numeric.append(c)
            d[c] = f
        else:
            dropped.append(c)
    if not numeric:
        raise ValueError(f"{path}: no numeric sample columns after the id column {ids!r}; columns are {cols}")
    counts = np.column_stack([d[c] for c in numeric])
    layers = {layer or "counts": counts}
    if length is not None:
        layers["length"] = np.repeat(length[:, None], len(numeric), axis=1)
    notes = [f"ignored non-numeric column(s) {dropped}"] if dropped else []
    samples = [_bam_stem(c) if _BAM_RE.search(str(c)) else str(c) for c in numeric] if header else \
        [f"sample{i + 1}" for i in range(len(numeric))]
    return _Table(d[ids], layers, samples, names=None if name_col is None else d[name_col], level=level,
                  fmt="matrix", notes=notes)


# ── combining files ─────────────────────────────────────────────────────────

def _combine(tables: List[_Table], paths: List[str], verbose: bool):
    """Per-sample files side by side (ids identical, else the union with NaN);
    multi-sample matrices stacked as layers."""
    per_sample = [t for t in tables if len(t.samples) == 1 and t.fmt != "matrix"]
    matrices = [t for t in tables if t not in per_sample]
    seen: Dict[str, int] = {}
    for t in tables:
        for n in t.notes:
            seen[n] = seen.get(n, 0) + 1
    notes = [n if c == 1 else f"{n} — {c} files" for n, c in seen.items()]
    if per_sample and matrices:
        raise ValueError("cannot mix per-sample quantification files with multi-sample matrices in one call")
    if per_sample:
        first = per_sample[0]
        same = all(len(t.ids) == len(first.ids) and np.array_equal(t.ids, first.ids) for t in per_sample[1:])
        layer_names = [k for k in first.layers if all(k in t.layers for t in per_sample)]
        if same:
            ids = first.ids
            layers = {k: np.concatenate([t.layers[k] for t in per_sample], axis=1) for k in layer_names}
        else:
            seen, ids = set(), []
            for t in per_sample:
                for x in t.ids.tolist():
                    if x not in seen:
                        seen.add(x)
                        ids.append(x)
            ids = np.asarray(ids, dtype=object)
            pos = {x: i for i, x in enumerate(ids.tolist())}
            layers = {k: np.full((len(ids), len(per_sample)), np.nan) for k in layer_names}
            for j, t in enumerate(per_sample):
                rows = np.fromiter((pos[x] for x in t.ids.tolist()), np.int64, len(t.ids))
                for k in layer_names:
                    layers[k][rows, j] = t.layers[k][:, 0]
            notes.append(f"the files do not list the same features: {len(ids)} in the union, NaN where absent")
        level = first.level
        fmt = first.fmt if len({t.fmt for t in per_sample}) == 1 else "mixed"
        return _Table(ids, layers, [t.samples[0] for t in per_sample], level=level, fmt=fmt, notes=notes)
    first = matrices[0]
    layers = dict(first.layers)
    for t, p in zip(matrices[1:], paths[1:]):
        if not np.array_equal(t.ids, first.ids):
            raise ValueError(f"{p}: its features differ from {paths[0]}; matrices passed together must list "
                             f"the same features in the same order")
        order = [t.samples.index(s) for s in first.samples] if set(t.samples) == set(first.samples) else None
        if order is None:
            raise ValueError(f"{p}: samples {t.samples[:3]}... differ from {paths[0]}")
        for k, v in t.layers.items():
            key = k if k not in layers else f"{k}_{len(layers)}"
            layers[key] = v[:, order]
    return _Table(first.ids, layers, first.samples, names=first.names, level=first.level, fmt="matrix", notes=notes)


# ── alignment to a Genes / Loci table ───────────────────────────────────────

def _strip_version(a) -> np.ndarray:
    return np.asarray([_VERSION_RE.sub(r"\1", str(x)) for x in a], dtype=object)


def _pipe_cut(a) -> np.ndarray:
    return np.asarray([str(x).split("|", 1)[0] for x in a], dtype=object)


def _region_key(a) -> np.ndarray:
    out = []
    for x in a:
        m = _REGION_RE.match(str(x))
        out.append(f"{m.group(1)}:{m.group(2)}-{m.group(3)}" if m else str(x))
    return np.asarray(out, dtype=object)


def _first_index(keys) -> dict:
    ix = {}
    for i, k in enumerate(keys.tolist()):
        ix.setdefault(k, i)
    return ix


def _candidates(rows, level, key: str):
    """(rule name, table keys, transform of the input ids) to try, in order."""
    from .genes import Genes
    if isinstance(rows, Genes):
        if level == "transcript":
            tid = rows.transcripts.cols["transcript_id"]
            base = [("transcript_id", tid, None), ("transcript_id, pipe header cut", tid, _pipe_cut),
                    ("transcript_id, version stripped", _strip_version(tid), _strip_version),
                    ("transcript_id, pipe header cut + version stripped", _strip_version(tid),
                     lambda a: _strip_version(_pipe_cut(a)))]
            by_name = []
        else:
            gid = rows.genes.cols["gene_id"]
            base = [("gene_id", gid, None), ("gene_id, pipe header cut", gid, _pipe_cut),
                    ("gene_id, version stripped", _strip_version(gid), _strip_version),
                    ("gene_id, pipe header cut + version stripped", _strip_version(gid),
                     lambda a: _strip_version(_pipe_cut(a)))]
            by_name = [("gene_name", rows.genes.cols["gene_name"], None)]
        if key == "auto":
            return base + by_name
        if key in ("gene_id", "transcript_id", "id"):
            return base
        if key == "gene_name":
            if level == "transcript":
                raise ValueError("key='gene_name' matches gene rows; use level='gene'")
            return by_name
        raise ValueError(f"key must be 'auto', 'gene_id', 'transcript_id' or 'gene_name' for a Genes table, got {key!r}")
    L = rows
    cands = []
    if "name" in L.cols:
        cands.append(("name", L.cols["name"], None))
    cands += [("uid", L.uid, None), ("region (chrom:start-end)", _region_key(L.uid), _region_key)]
    if key == "auto":
        return cands
    if key in ("name", "uid", "region"):
        pick = {"name": "name", "uid": "uid", "region": "region (chrom:start-end)"}[key]
        sel = [c for c in cands if c[0] == pick]
        if not sel:
            raise KeyError(f"the Loci has no {key!r} column")
        return sel
    raise ValueError(f"key must be 'auto', 'name', 'uid' or 'region' for a Loci, got {key!r}")


def _match(ids, rows, level, key: str):
    """(row index per input id, -1 = none; rule used). Every rule is tried and
    the one matching most ids wins (the first on a tie), so a file mixing
    versioned and pipe-headed ids still lands on one rule."""
    best = None
    for rule, table_keys, transform in _candidates(rows, level, key):
        q = ids if transform is None else transform(ids)
        ix = _first_index(np.asarray(table_keys, dtype=object))
        get = ix.get
        idx = np.fromiter((get(x, -1) for x in q.tolist()), np.int64, len(q))
        n_hit = int((idx >= 0).sum())
        if best is None or n_hit > best[2]:
            best = (idx, rule, n_hit)
        if n_hit == len(ids):
            break
    return best[0], best[1], best[2]


def _summarise(layers: Dict[str, np.ndarray], target: np.ndarray, n_targets: int, *, how: str = "counts",
               fill: str = "plain"):
    """Collapse feature rows onto ``target`` rows (``-1`` = drop): counts and
    TPM sum; the length is the abundance-weighted mean of the members'
    lengths, as ``tximport::summarizeToGene``. Where a target has no
    abundance in a sample, ``fill="plain"`` (aligning input rows to table
    rows) takes the plain mean of the members' lengths in that sample — for
    a one-to-one match, the member's own length — and ``fill="tximport"``
    (transcripts to genes) applies tximport's rule: the geometric mean of the
    target's lengths in the samples that have abundance, or, with none
    anywhere, the plain mean over members and samples. ``how`` rescales the
    counts from TPM (``scaled_tpm`` / ``length_scaled_tpm``). A target no
    quantified feature points at is NaN."""
    ok = target >= 0
    t = target[ok]
    main = "counts" if "counts" in layers else next(iter(layers))
    m = layers[main].shape[1]
    present = {k: ~np.isnan(v[ok]) for k, v in layers.items()}
    hit = np.zeros((n_targets, m))
    np.add.at(hit, t, present[main].astype(float))

    def total(k):
        v = np.where(present[k], layers[k][ok], 0.0)
        out = np.zeros((n_targets, m))
        np.add.at(out, t, v)
        return out

    out = {}
    for k in layers:
        if k != "length":
            out[k] = total(k)
    if "length" in layers:
        ln = np.where(present["length"], layers["length"][ok], 0.0)
        n_len = np.zeros((n_targets, m))
        np.add.at(n_len, t, present["length"].astype(float))
        plain = np.zeros((n_targets, m))
        np.add.at(plain, t, ln)
        with np.errstate(invalid="ignore", divide="ignore"):
            plain = plain / n_len                                  # plain mean of the members' lengths, per sample
        if "tpm" in layers:
            w = np.where(present["tpm"] & present["length"], layers["tpm"][ok], 0.0)
            num = np.zeros((n_targets, m))
            np.add.at(num, t, w * ln)
            den = np.zeros((n_targets, m))
            np.add.at(den, t, w)
            with np.errstate(invalid="ignore", divide="ignore"):
                length = num / den                                 # NaN where the target has no abundance
                gap = ~(den > 0) & (hit > 0)
                if gap.any() and fill == "plain":
                    length = np.where(gap, plain, length)
                elif gap.any():
                    # tximport's replaceMissingLength: a target with abundance in other samples gets the
                    # geometric mean of its lengths there; one with none anywhere gets the plain mean of its
                    # members' lengths (each averaged over the samples first)
                    has = (den > 0).any(axis=1)
                    geo = np.full(n_targets, np.nan)
                    if has.any():
                        lg = np.where(den[has] > 0, np.log(length[has]), 0.0)
                        geo[has] = np.exp(lg.sum(axis=1) / (den[has] > 0).sum(axis=1))
                    n_pl = (~np.isnan(plain)).sum(axis=1)
                    pl_mean = np.where(n_pl > 0, np.nan_to_num(plain).sum(axis=1) / np.maximum(n_pl, 1), np.nan)
                    fill = np.where(has, geo, pl_mean)[:, None]
                    length = np.where(gap, fill, length)
        else:
            length = plain
        out["length"] = length
    if how != "counts":
        if "tpm" not in out or "counts" not in out:
            raise ValueError(f"how={how!r} needs counts and a tpm layer (salmon / kallisto / RSEM input)")
        lib = out["counts"].sum(axis=0)
        new = out["tpm"]
        if how == "length_scaled_tpm":
            if "length" not in out:
                raise ValueError("how='length_scaled_tpm' needs a length layer")
            with np.errstate(invalid="ignore"):
                new = new * np.nanmean(np.where(hit > 0, out["length"], np.nan), axis=1)[:, None]
        with np.errstate(invalid="ignore", divide="ignore"):
            scale = np.where(new.sum(axis=0) > 0, lib / new.sum(axis=0), 0.0)
        out["counts"] = new * scale
    none = hit == 0
    for k in out:
        out[k] = np.where(none, np.nan, out[k])
    return out


# ── the table ───────────────────────────────────────────────────────────────

class Expression(TableMixin):
    """A features x samples matrix aligned to the rows of a Genes or Loci table.

    ``E.counts`` (always), ``E.tpm``, ``E.length`` and ``E.layers[name]`` are
    float64 arrays of shape ``(len(E), E.n_samples)``; ``E.samples`` is a
    pandas DataFrame indexed by sample name; ``E.rows`` is the ``Genes`` or
    ``Loci`` the rows line up with (``None`` before :meth:`align`) and
    ``E.table`` that object's row table. As a table it is the wide frame of
    the default layer (``id``, ``gene_name``, one column per sample).
    """

    def __init__(self, layers: Dict[str, np.ndarray], samples, *, ids, names=None, rows=None,
                 level: Optional[str] = None, filename: Optional[str] = None, meta: Optional[dict] = None,
                 unmatched=None):
        import pandas as pd
        if not layers:
            raise ValueError("Expression: at least one layer (normally 'counts') is required")
        self.layers = {k: np.asarray(v, dtype=np.float64) for k, v in layers.items()}
        self._default = "counts" if "counts" in self.layers else next(iter(self.layers))
        if self.layers[self._default].ndim != 2:
            raise ValueError(f"Expression: layers are 2-D (features x samples), got shape "
                             f"{self.layers[self._default].shape}")
        n, m = self.layers[self._default].shape
        for k, v in self.layers.items():
            if v.shape != (n, m):
                raise ValueError(f"Expression: layer {k!r} has shape {v.shape}, {self._default} has {(n, m)}")
        self.ids = np.asarray(ids, dtype=object)
        if len(self.ids) != n:
            raise ValueError(f"Expression: {len(self.ids)} ids for {n} rows")
        self.names = None if names is None else np.asarray(names, dtype=object)
        if self.names is not None and len(self.names) != n:
            raise ValueError(f"Expression: {len(self.names)} names for {n} rows")
        if isinstance(samples, pd.DataFrame):
            s = samples.copy()
        else:
            s = pd.DataFrame(index=pd.Index([str(x) for x in samples], dtype=object))
        s.index = pd.Index([str(x) for x in s.index], dtype=object, name="sample")
        if len(s) != m:
            raise ValueError(f"Expression: {len(s)} samples for {m} columns")
        if s.index.has_duplicates:
            dup = s.index[s.index.duplicated()].tolist()[:3]
            raise ValueError(f"Expression: duplicated sample names {dup}; pass sample_names= to rename")
        self._samples = s
        if level is not None and level not in LEVELS:
            raise ValueError(f"level must be one of {LEVELS}, got {level!r}")
        self.level = level
        self.rows = rows
        if rows is not None:
            if len(self.table) != n:
                raise ValueError(f"Expression: {n} rows for a {self.level} table of {len(self.table)}; "
                                 f"build with Expression.make(..., genes=) or .align()")
        self.filename = filename
        self.meta = dict(meta or {})
        self.unmatched = np.zeros(0, dtype=object) if unmatched is None else np.asarray(unmatched, dtype=object)
        self._by_id = None

    # ── construction ──────────────────────────────────────────────────────
    @classmethod
    def make(cls, source, *, genes=None, loci=None, level: Optional[str] = None, how: str = "length_scaled_tpm",
             samples=None, key: str = "auto", sample_names=None, strand: str = "auto", strict: bool = False,
             verbose: bool = True, backend: Optional[str] = None) -> "Expression":
        """Read quantifications and (with ``genes=`` or ``loci=``) align them.

        ``source``: a glob or list of per-sample files (salmon / kallisto /
        RSEM / STAR / HTSeq; the sample is the directory or file name), a
        ``{sample: path}`` dict, or one or more multi-sample matrices
        (featureCounts, nf-core/rnaseq, any id + samples table; several
        files of one run fill several layers). ``level`` is the level of
        the result: a transcript-level input with ``genes=`` is summarised
        to ``"gene"`` by default (``how``, see :meth:`to_genes`);
        ``level="transcript"`` keeps it aligned to ``genes.transcripts``.
        ``samples`` is a DataFrame or TSV/CSV of sample metadata (index =
        sample names). ``key`` picks the matching rule (``"auto"``,
        ``"gene_id"``, ``"transcript_id"``, ``"gene_name"``; ``"name"`` /
        ``"uid"`` / ``"region"`` for a Loci). ``strand`` is for STAR tables.
        """
        paths, names = _paths(source)
        if sample_names is not None:
            if len(sample_names) != len(paths) and len(paths) != 1:
                raise ValueError(f"{len(sample_names)} sample_names for {len(paths)} files")
            if len(paths) == len(sample_names):
                names = [str(s) for s in sample_names]
        tables = [read_table(p, strand=strand, sample=s, backend=backend) for p, s in zip(paths, names)]
        t = _combine(tables, paths, verbose) if len(tables) > 1 else tables[0]
        if sample_names is not None and len(paths) == 1:
            if len(sample_names) != len(t.samples):
                raise ValueError(f"{len(sample_names)} sample_names for {len(t.samples)} samples")
            t.samples = [str(s) for s in sample_names]
        E = cls(t.layers, t.samples, ids=t.ids, names=t.names, level=t.level,
                filename=paths[0] if len(paths) == 1 else os.path.commonpath([os.path.abspath(p) for p in paths]),
                meta={"format": t.fmt, "sources": len(paths)})
        if verbose:
            what = f"{t.fmt} ({os.path.basename(paths[0])})" if len(paths) == 1 else f"{len(paths)} {t.fmt} files"
            print(f"[INFO] Expression: {len(E):,} {t.level or 'feature'}s × {E.n_samples} samples from {what}; "
                  f"layers {', '.join(E.layers)}")
            for note in t.notes:
                print(f"[INFO]   {note}")
        if samples is not None:
            E.samples = _samples_frame(samples, E.samples.index)
        if genes is not None or loci is not None:
            E = E.align(genes=genes, loci=loci, level=level, how=how, key=key, strict=strict, verbose=verbose)
        elif level is not None:
            if E.level is not None and E.level != level and not (E.level == "transcript" and level == "gene"):
                raise ValueError(f"the input is {E.level}-level; level={level!r} needs genes= to summarise")
            if level == "gene" and E.level == "transcript":
                raise ValueError("summarising transcripts to genes needs genes=")
            E.level = level
        return E

    @classmethod
    def from_frame(cls, df, *, ids=None, names=None, layer: str = "counts", genes=None, loci=None,
                   level: Optional[str] = None, how: str = "length_scaled_tpm", samples=None, key: str = "auto",
                   strict: bool = False, verbose: bool = True) -> "Expression":
        """From any table (pandas, polars, pyarrow, dict of columns, ...):
        ``ids`` names the feature-id column (default: the first), ``names`` a
        gene-name column (found by name otherwise), the numeric columns are
        the samples. ``layer`` names what the numbers are."""
        from .interop import frame
        f = frame(df)
        if f is None:
            raise TypeError(f"expected a table (pandas, polars, pyarrow, ...), got {type(df).__name__}")
        pdf = f.to_pandas()
        cols = [str(c) for c in pdf.columns]
        pdf.columns = cols
        d = {c: pdf[c].to_numpy() for c in cols}
        for c in cols:
            if d[c].dtype.kind not in "fiub":
                d[c] = d[c].astype(object)
        t = _matrix_table(cols, d, path="frame", level=level, ids=None if ids is None else str(ids),
                          name_col=None if names is None else str(names), layer=layer)
        E = cls(t.layers, t.samples, ids=t.ids, names=t.names, level=level, meta={"format": "frame"})
        if verbose:
            print(f"[INFO] Expression: {len(E):,} {level or 'feature'}s × {E.n_samples} samples from a frame; "
                  f"layers {', '.join(E.layers)}")
            for note in t.notes:
                print(f"[INFO]   {note}")
        if samples is not None:
            E.samples = _samples_frame(samples, E.samples.index)
        if genes is not None or loci is not None:
            E = E.align(genes=genes, loci=loci, level=level, how=how, key=key, strict=strict, verbose=verbose)
        return E

    @classmethod
    def from_anndata(cls, adata, *, layer: Optional[str] = None, genes=None, loci=None, level: Optional[str] = None,
                     how: str = "length_scaled_tpm", key: str = "auto", strict: bool = False,
                     verbose: bool = True) -> "Expression":
        """From an AnnData of samples (obs) x features (var): ``X`` (or
        ``layer``) becomes ``counts``, every other layer comes along, ``obs``
        becomes ``samples``, ``var['gene_name']`` the names when present."""
        X = adata.layers[layer] if layer else adata.X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        info = dict(adata.uns.get("genomeblocks", {})) if hasattr(adata, "uns") else {}
        others = [k for k in adata.layers.keys() if k is not None and k != layer]
        x_name = layer or info.get("X") or ("counts" if "counts" not in others else "X")
        layers = {x_name: X.T.astype(np.float64)}
        for k in others:
            v = adata.layers[k]
            v = v.toarray() if hasattr(v, "toarray") else np.asarray(v)
            layers[k] = v.T.astype(np.float64)
        if "counts" in layers and x_name != "counts":          # counts first: the default layer
            layers = {"counts": layers.pop("counts"), **layers}
        level = level or (info.get("level") or None)
        names = None
        for c in _NAME_COLUMNS:
            if c in adata.var.columns:
                names = adata.var[c].to_numpy(object)
                break
        E = cls(layers, adata.obs, ids=adata.var_names.to_numpy(object), names=names, level=level,
                meta={"format": "anndata"})
        if verbose:
            print(f"[INFO] Expression: {len(E):,} {level or 'feature'}s × {E.n_samples} samples from an AnnData; "
                  f"layers {', '.join(E.layers)}")
        if genes is not None or loci is not None:
            E = E.align(genes=genes, loci=loci, level=level, how=how, key=key, strict=strict, verbose=verbose)
        return E

    # ── alignment ─────────────────────────────────────────────────────────
    def align(self, genes=None, loci=None, *, level: Optional[str] = None, how: str = "length_scaled_tpm",
              key: str = "auto", strict: bool = False, verbose: bool = True) -> "Expression":
        """The same numbers on the rows of ``genes`` (``genes.genes`` or
        ``genes.transcripts``) or of ``loci``. Input features matching no
        row are dropped (``unmatched``), rows without input are NaN
        (``missing``); several inputs on one row are summed (lengths
        abundance-averaged). A transcript-level input aligned at
        ``level="gene"`` is summarised with :meth:`to_genes` (``how``)."""
        from .genes import Genes
        if (genes is None) == (loci is None):
            raise ValueError("pass exactly one of genes= or loci=")
        if genes is not None and not isinstance(genes, Genes):
            raise TypeError(f"genes= must be a Genes, got {type(genes).__name__}")
        if loci is not None:
            from .interop import as_loci
            loci = as_loci(loci)
        rows = genes if genes is not None else loci
        want = level
        if loci is not None:
            if want not in (None, "loci"):
                raise ValueError("level= must be 'loci' (or None) with loci=")
            match_level = want = "loci"
        else:
            own = self.level
            if want is None:
                want = "gene" if own in (None, "gene", "transcript") else own
            if want not in ("gene", "transcript"):
                raise ValueError(f"level must be 'gene' or 'transcript' with genes=, got {want!r}")
            if own == "gene" and want == "transcript":
                raise ValueError("a gene-level table cannot be aligned to transcripts")
            match_level = "transcript" if own == "transcript" else want
        idx, rule, n_hit = _match(self.ids, rows, match_level, key)
        if n_hit == 0:
            raise ValueError(f"none of the {len(self.ids)} input ids matches the {match_level} table (tried {rule!r}); "
                             f"first ids: {self.ids[:3].tolist()} — pass key= to name the id kind")
        table = _row_table(rows, match_level)
        unmatched = self.ids[idx < 0]
        if strict and len(unmatched):
            raise ValueError(f"{len(unmatched)} of {len(self.ids)} input ids match no {match_level} row "
                             f"(first: {unmatched[:3].tolist()}); strict=False keeps going")
        n_dup = int(n_hit - len(np.unique(idx[idx >= 0])))
        layers = _summarise(self.layers, idx, len(table), how="counts")
        E = Expression(layers, self._samples, ids=_row_ids(rows, match_level), names=_row_names(rows, match_level),
                       rows=rows, level=match_level, filename=self.filename,
                       meta={**self.meta, "key_rule": rule, "rows": _rows_label(rows, match_level)}, unmatched=unmatched)
        if verbose:
            miss = int(E.missing.sum())
            print(f"[INFO] Aligned to {_rows_label(rows, match_level)} by {rule}: {n_hit:,}/{len(self.ids):,} input ids matched"
                  + (f" ({len(unmatched):,} dropped)" if len(unmatched) else "")
                  + (f", {n_dup:,} shared a row and were summed" if n_dup else "")
                  + f"; {miss:,} of {len(table):,} rows not quantified")
        if match_level == "transcript" and want == "gene":
            E = E.to_genes(how=how, verbose=verbose)
        return E

    def to_genes(self, how: str = "length_scaled_tpm", *, verbose: bool = True) -> "Expression":
        """Transcripts to genes over ``transcripts['gene']``, as
        ``tximport::summarizeToGene``: counts and TPM sum per gene; the gene
        length is the TPM-weighted mean of its transcripts' effective
        lengths (plain mean where TPM is 0). ``how="counts"`` keeps the
        summed counts (use length offsets in a model); ``"scaled_tpm"``
        rescales gene TPM to the library size; ``"length_scaled_tpm"``
        (default) also multiplies by the gene's mean length first — the
        counts tximport recommends for DESeq2 / edgeR without offsets."""
        from .genes import Genes
        if how not in HOW:
            raise ValueError(f"how must be one of {HOW}, got {how!r}")
        if not isinstance(self.rows, Genes) or self.level != "transcript":
            raise ValueError("to_genes needs a transcript-level table aligned to a Genes (Expression.make(..., "
                             "genes=genes, level='transcript'))")
        G = self.rows
        target = G.transcripts.cols["gene"].astype(np.int64)
        layers = _summarise(self.layers, target, len(G.genes), how=how, fill="tximport")
        E = Expression(layers, self._samples, ids=G.genes.cols["gene_id"], names=G.genes.cols["gene_name"],
                       rows=G, level="gene", filename=self.filename,
                       meta={**self.meta, "counts_from": how, "rows": "genes.genes"}, unmatched=self.unmatched)
        if verbose:
            print(f"[INFO] Summarised {int((~self.missing).sum()):,} transcripts to {int((~E.missing).sum()):,} "
                  f"genes ({how})")
        return E

    # ── samples ───────────────────────────────────────────────────────────
    @property
    def samples(self):
        """Sample metadata: a DataFrame indexed by sample name, in column order."""
        return self._samples

    @samples.setter
    def samples(self, df):
        import pandas as pd
        if not isinstance(df, pd.DataFrame):
            raise TypeError("samples must be a pandas DataFrame indexed by sample name")
        if len(df) != self.n_samples:
            raise ValueError(f"{len(df)} rows for {self.n_samples} samples")
        df = df.copy()
        df.index = pd.Index([str(x) for x in df.index], dtype=object, name="sample")
        if set(df.index) == set(self._samples.index) and not df.index.equals(self._samples.index):
            df = df.loc[self._samples.index]
        self._samples = df

    @property
    def n_samples(self) -> int:
        return self.layers[self._default].shape[1]

    @property
    def sample_names(self) -> list:
        return self._samples.index.tolist()

    # ── layers ────────────────────────────────────────────────────────────
    @property
    def default_layer(self) -> str:
        """``'counts'`` when the table has counts, else its first layer."""
        return self._default

    @property
    def X(self) -> np.ndarray:
        """The default layer."""
        return self.layers[self._default]

    @property
    def counts(self) -> np.ndarray:
        if "counts" not in self.layers:
            raise AttributeError(f"this table has no counts layer (it holds {', '.join(self.layers)}); counts are "
                                 f"needed for library sizes, normalisation and differential expression")
        return self.layers["counts"]

    @property
    def length(self) -> np.ndarray:
        if "length" not in self.layers:
            raise AttributeError("this table has no length layer (featureCounts / salmon / kallisto / RSEM give one)")
        return self.layers["length"]

    @property
    def tpm(self) -> np.ndarray:
        """The TPM layer, or TPM computed from counts and effective lengths."""
        if "tpm" in self.layers:
            return self.layers["tpm"]
        if "length" not in self.layers:
            raise AttributeError("no tpm layer and no length layer to compute it from")
        with np.errstate(invalid="ignore", divide="ignore"):
            rate = self.layers["counts"] / self.layers["length"]
            return rate / np.nansum(rate, axis=0) * 1e6

    def layer(self, name: Optional[str] = None) -> np.ndarray:
        """A layer by name (``None``: the default layer)."""
        if name is None:
            return self.layers[self._default]
        if name == "counts":
            return self.counts
        if name not in self.layers:
            raise KeyError(f"no layer {name!r}; layers are {list(self.layers)}")
        return self.layers[name]

    @property
    def missing(self) -> np.ndarray:
        """Rows the input did not quantify (NaN in every layer)."""
        return np.isnan(self.layers[self._default]).all(axis=1) if self.n_samples else np.zeros(len(self), bool)

    @property
    def table(self):
        """The Loci the rows are aligned to (``genes.genes``, ``genes.transcripts`` or the Loci)."""
        return None if self.rows is None else _row_table(self.rows, self.level)

    @property
    def lib_size(self) -> np.ndarray:
        """Total counts per sample over the quantified rows."""
        from ._norm import lib_sizes
        return lib_sizes(self.counts)

    # ── normalisation ─────────────────────────────────────────────────────
    def cpm(self, *, log: bool = False, prior: float = 2.0, norm: Optional[str] = None) -> np.ndarray:
        """Counts per million (edgeR's ``cpm``); ``norm="tmm"`` / ``"size"``
        scales the library sizes by :meth:`tmm` / :meth:`size_factors`;
        ``log=True`` gives log2 with a library-scaled prior count."""
        from ._norm import cpm
        f = None
        if norm == "tmm":
            f = self.tmm()
        elif norm in ("size", "size_factors", "deseq2"):
            f = self.size_factors()
        elif norm not in (None, "none"):
            raise ValueError(f"norm must be None, 'tmm' or 'size', got {norm!r}")
        return cpm(self.counts, factors=f, log=log, prior=prior)

    def log_cpm(self, *, prior: float = 2.0, norm: Optional[str] = None) -> np.ndarray:
        return self.cpm(log=True, prior=prior, norm=norm)

    def tmm(self, **kw) -> np.ndarray:
        """edgeR's TMM factors, one per sample (geometric mean 1)."""
        from ._norm import tmm_factors
        return tmm_factors(self.counts, **kw)

    def size_factors(self) -> np.ndarray:
        """DESeq2's median-of-ratios size factors, one per sample."""
        from ._norm import size_factors
        return size_factors(self.counts)

    def filter_by_expr(self, group=None, **kw) -> np.ndarray:
        """edgeR's ``filterByExpr`` mask over the rows; ``group`` is a column
        of ``samples`` (by name) or an array."""
        from ._norm import filter_by_expr
        if isinstance(group, str):
            group = self._samples[group].to_numpy()
        return filter_by_expr(self.counts, group=group, **kw)

    # ── rows ──────────────────────────────────────────────────────────────
    def _key_index(self):
        if self._by_id is None:
            ix = _first_index(self.ids)
            if self.names is not None:
                for i, n in enumerate(self.names.tolist()):
                    if n not in (None, "") and n not in ix:
                        ix[n] = i
            self._by_id = ix
        return self._by_id

    def row(self, key) -> int:
        """The row of a gene / transcript id, a gene name or a Loci uid."""
        from .genes import Genes
        if isinstance(self.rows, Genes) and self.level == "gene":
            try:
                return int(self.rows._key_index()[key])
            except KeyError:
                pass
        if self.rows is not None and self.level == "loci":
            try:
                return int(self.rows.row(key))
            except KeyError:
                pass
        try:
            return int(self._key_index()[key])
        except KeyError:
            raise KeyError(f"no row {key!r} (ids and gene names are keys)") from None

    def get(self, key, layer: Optional[str] = None):
        """One feature's values over the samples as a pandas Series (``layer``: the default layer)."""
        import pandas as pd
        i = self.row(key)
        return pd.Series(self.layer(layer)[i], index=self._samples.index, name=str(key))

    def lookup(self, keys, values=None, *, layer: Optional[str] = None, agg: str = "mean", missing=np.nan) -> np.ndarray:
        """``values`` (one per row; default the ``agg`` of ``layer`` over the
        samples) picked for each of ``keys`` (gene names / ids / uids); an
        unknown or empty key gives ``missing``. For ``A.vp[...] =
        E.lookup(A.vp.gene, res.log2fc)``-style joins."""
        if values is None:
            red = {"mean": np.nanmean, "median": np.nanmedian, "sum": np.nansum, "max": np.nanmax}[agg]
            with np.errstate(invalid="ignore"):
                values = red(self.layer(layer), axis=1)
        v = np.asarray(values.to_numpy() if hasattr(values, "to_numpy") else values)
        if len(v) != len(self):
            raise ValueError(f"{len(v)} values for {len(self)} rows")
        ix = self._key_index()
        from .genes import Genes
        if isinstance(self.rows, Genes) and self.level == "gene":
            ix = {**ix, **self.rows._key_index()}
        get = ix.get
        keys = [keys] if isinstance(keys, str) else list(np.asarray(keys, dtype=object).tolist())
        rows = np.fromiter((get(k, -1) if k not in (None, "") else -1 for k in keys), np.int64, len(keys))
        out = np.where(rows >= 0, v[np.maximum(rows, 0)], missing) if v.dtype.kind == "f" else \
            np.array([v[r] if r >= 0 else missing for r in rows], dtype=object)
        return out

    def take(self, idx) -> "Expression":
        """Rows ``idx``: the alignment is dropped (the result no longer lines
        up with the whole gene table); ids and names stay."""
        if isinstance(idx, list):
            idx = np.asarray(idx)
        return Expression({k: v[idx] for k, v in self.layers.items()}, self._samples, ids=self.ids[idx],
                          names=None if self.names is None else self.names[idx], level=self.level,
                          filename=self.filename, meta={**self.meta, "rows": None})

    def take_samples(self, idx) -> "Expression":
        """Samples ``idx`` (names, mask, rows or slice); the row alignment is kept."""
        import pandas as pd
        if isinstance(idx, str):
            idx = [idx]
        if isinstance(idx, (list, np.ndarray, pd.Index)) and len(idx) and isinstance(np.asarray(idx).ravel()[0], str):
            pos = {s: i for i, s in enumerate(self.sample_names)}
            idx = np.array([pos[s] for s in idx], np.int64)
        if isinstance(idx, pd.Series):
            idx = idx.to_numpy()
        return Expression({k: v[:, idx] for k, v in self.layers.items()}, self._samples.iloc[idx], ids=self.ids,
                          names=self.names, rows=self.rows, level=self.level, filename=self.filename,
                          meta=dict(self.meta), unmatched=self.unmatched)

    def __getitem__(self, key):
        """``E['TP53']`` a Series over samples; ``E[['TP53', 'MYC']]`` a frame;
        ``E[i]`` row ``i``; ``E[mask]`` / ``E[rows]`` a row subset (alignment
        dropped); ``E[:, samples]`` a sample subset."""
        import pandas as pd
        if isinstance(key, tuple):
            if len(key) != 2:
                raise IndexError("use E[rows, samples]")
            r, s = key
            out = self if (isinstance(r, slice) and r == slice(None)) else self.take(r)
            return out if (isinstance(s, slice) and s == slice(None)) else out.take_samples(s)
        if isinstance(key, (int, np.integer)):
            n = len(self)
            if not -n <= key < n:
                raise IndexError(key)
            i = int(key) % n
            return pd.Series(self.X[i], index=self._samples.index, name=str(self.ids[i]))
        if isinstance(key, str):
            return self.get(key)
        if isinstance(key, (list, tuple, np.ndarray)) and len(key) and isinstance(np.asarray(key).ravel()[0], str):
            rows = [self.row(k) for k in key]
            return pd.DataFrame(self.X[rows], index=pd.Index(list(key), dtype=object), columns=self._samples.index)
        return self.take(key)

    def __contains__(self, key) -> bool:
        try:
            self.row(key)
            return True
        except KeyError:
            return False

    def __len__(self):
        return len(self.ids)

    def __iter__(self):
        return iter(self.ids.tolist())

    # ── like a DataFrame ─────────────────────────────────────────────────
    @property
    def shape(self):
        """(features, samples) — the matrix shape."""
        return self.layers[self._default].shape

    @property
    def columns(self) -> list:
        """The sample names (the columns of the matrix)."""
        return self.sample_names

    def head(self, n: int = 5):
        return self.to_pandas().head(n)

    def tail(self, n: int = 5):
        return self.to_pandas().tail(n)

    def describe(self):
        import pandas as pd
        lib = self.lib_size if "counts" in self.layers else np.zeros(0)
        q = (np.nanmin(lib), np.nanmedian(lib), np.nanmax(lib)) if len(lib) else (0, 0, 0)
        rows = [("rows", len(self)), ("level", self.level or "unknown"), ("samples", self.n_samples),
                ("layers", ", ".join(self.layers)), ("aligned to", self.meta.get("rows") or "nothing"),
                ("quantified rows", int((~self.missing).sum())), ("rows not quantified", int(self.missing.sum())),
                ("input ids dropped", len(self.unmatched)), ("id rule", self.meta.get("key_rule", "")),
                ("counts from", self.meta.get("counts_from", self.meta.get("format", ""))),
                ("library size min / median / max", f"{q[0]:,.0f} / {q[1]:,.0f} / {q[2]:,.0f}")]
        return pd.DataFrame(rows, columns=["", "value"]).set_index("").rename_axis(None)

    summary = describe

    # ── exits ─────────────────────────────────────────────────────────────
    def _id_column(self) -> str:
        return {"gene": "gene_id", "transcript": "transcript_id", "loci": "uid"}.get(self.level, "id") \
            if self.rows is not None else "id"

    def to_arrow(self, layer: Optional[str] = None):
        import pyarrow as pa
        arrs = {self._id_column(): pa.array(self.ids.tolist(), pa.string())}
        if self.names is not None:
            arrs["gene_name"] = pa.array([None if x is None else str(x) for x in self.names.tolist()], pa.string())
        X = self.layer(layer)
        for j, s in enumerate(self.sample_names):
            arrs[s] = pa.array(X[:, j])
        return pa.table(arrs)

    def to_pandas(self, layer: Optional[str] = None, *, long: bool = False):
        """The wide table (id, gene_name, one column per sample) of ``layer``
        (default: the default layer); ``long=True`` gives one row per
        (feature, sample) with the sample metadata joined on, for seaborn
        and friends."""
        import pandas as pd
        layer = layer or self._default
        X = self.layer(layer)
        d = {self._id_column(): self.ids}
        if self.names is not None:
            d["gene_name"] = self.names
        if not long:
            for j, s in enumerate(self.sample_names):
                d[s] = X[:, j]
            return pd.DataFrame(d)
        n, m = X.shape
        out = pd.DataFrame({k: np.repeat(v, m) for k, v in d.items()})
        out["sample"] = np.tile(np.asarray(self.sample_names, dtype=object), n)
        out[layer] = X.ravel()
        if self._samples.shape[1]:
            out = out.join(self._samples, on="sample")
        return out

    def to_polars(self, layer: Optional[str] = None):
        import polars as pl
        return pl.from_arrow(self.to_arrow(layer))

    def to_anndata(self, layer: Optional[str] = None):
        """AnnData of samples (obs) x features (var): ``X`` is ``layer`` (default
        counts), every other layer is kept, ``obs`` is ``samples``."""
        import anndata as ad
        import pandas as pd
        main = layer or self._default
        var = pd.DataFrame(index=pd.Index(self.ids.astype(str), dtype=object, name=self._id_column()))
        if self.names is not None:
            var["gene_name"] = self.names
        if self.rows is not None:
            var["missing"] = self.missing
        obs = self._samples.copy()
        obs.index = obs.index.astype(str)
        a = ad.AnnData(X=self.layer(main).T.copy(), obs=obs, var=var)
        for k, v in self.layers.items():                    # assigned after construction: anndata 0.13 adds a
            if k != main:                                   # spurious None layer when they go through the constructor
                a.layers[k] = v.T.copy()
        a.uns["genomeblocks"] = {"X": main, "level": self.level or "", "default_layer": self._default}
        return a

    def save(self, path: str):
        """``path/``: ``samples.parquet``, ``features.parquet`` (id, gene_name,
        missing), one ``<layer>.parquet`` per layer, ``meta.json``."""
        import json
        import pyarrow as pa
        import pyarrow.parquet as pq
        os.makedirs(path, exist_ok=True)
        s = self._samples.reset_index()
        pq.write_table(pa.Table.from_pandas(s, preserve_index=False), os.path.join(path, "samples.parquet"))
        feats = {"id": pa.array(self.ids.tolist(), pa.string())}
        if self.names is not None:
            feats["gene_name"] = pa.array([None if x is None else str(x) for x in self.names.tolist()], pa.string())
        feats["missing"] = pa.array(self.missing)
        pq.write_table(pa.table(feats), os.path.join(path, "features.parquet"))
        for k, v in self.layers.items():
            pq.write_table(pa.table({s: pa.array(v[:, j]) for j, s in enumerate(self.sample_names)}),
                           os.path.join(path, f"{k}.parquet"))
        meta = {**self.meta, "level": self.level, "layers": list(self.layers), "filename": self.filename,
                "unmatched": self.unmatched.tolist()[:10_000], "format_version": "genomeblocks.Expression/1"}
        with open(os.path.join(path, "meta.json"), "w") as f:
            json.dump(meta, f)

    @classmethod
    def load(cls, path: str, *, genes=None, loci=None, key: str = "auto", verbose: bool = False) -> "Expression":
        """Read a directory written by :meth:`save`; ``genes=`` / ``loci=``
        re-align the rows (by id) to the table you pass."""
        import json
        import pyarrow.parquet as pq
        with open(os.path.join(path, "meta.json")) as f:
            meta = json.load(f)
        s = pq.read_table(os.path.join(path, "samples.parquet")).to_pandas().set_index("sample")
        feats = pq.read_table(os.path.join(path, "features.parquet"))
        ids = np.asarray(feats.column("id").to_pylist(), dtype=object)
        names = (np.asarray(feats.column("gene_name").to_pylist(), dtype=object)
                 if "gene_name" in feats.schema.names else None)
        layers = {}
        for k in meta["layers"]:
            t = pq.read_table(os.path.join(path, f"{k}.parquet"))
            layers[k] = np.column_stack([t.column(c).to_numpy(zero_copy_only=False) for c in s.index.tolist()]) \
                if len(s) else np.zeros((len(ids), 0))
        level = meta.get("level")
        E = cls(layers, s, ids=ids, names=names, level=level, filename=meta.get("filename"),
                meta={k: v for k, v in meta.items() if k not in ("layers", "level", "filename", "unmatched",
                                                                  "format_version", "rows")},
                unmatched=meta.get("unmatched"))
        if genes is not None or loci is not None:
            E = E.align(genes=genes, loci=loci, level=level if level != "loci" else None, key=key, verbose=verbose)
        return E

    # ── display ───────────────────────────────────────────────────────────
    @property
    def _what(self) -> str:
        return {"gene": "genes", "transcript": "transcripts", "loci": "loci"}.get(self.level, "features")

    def __repr__(self):
        return (f"Expression({len(self):,} {self._what} × {self.n_samples} samples, "
                f"layers=[{', '.join(self.layers)}]"
                + (f", aligned to {self.meta.get('rows')}" if self.rows is not None else "") + ")")

    def _repr_html_(self):
        from ._display import kv_html
        mb = sum(v.nbytes for v in self.layers.values()) / 1e6
        return kv_html(f"Expression · {len(self):,} {self._what} × {self.n_samples} samples", [
            ("layers", ", ".join(f"{k} ({v.shape[0]:,} × {v.shape[1]})" for k, v in self.layers.items())),
            ("aligned to", f"{self.meta.get('rows')} by {self.meta.get('key_rule', '?')}" if self.rows is not None
             else "nothing (call .align(genes=...))"),
            ("quantified rows", f"{int((~self.missing).sum()):,} of {len(self):,}"
             + (f" · {len(self.unmatched):,} input ids dropped" if len(self.unmatched) else "")),
            ("samples", ", ".join(self.sample_names[:6]) + (" ..." if self.n_samples > 6 else "")),
            ("sample columns", ", ".join(map(str, self._samples.columns)) or "none"),
            ("memory", f"{mb:.1f} MB"),
        ])


# ── helpers ─────────────────────────────────────────────────────────────────

def _paths(source):
    """(paths, sample names or None per path) for the forms ``make`` takes."""
    if isinstance(source, dict):
        paths = [str(p) for p in source.values()]
        return paths, [str(k) for k in source.keys()]
    if isinstance(source, (str, os.PathLike)):
        s = str(source)
        if any(ch in s for ch in "*?[") and not os.path.exists(s):
            paths = sorted(_glob.glob(s))
            if not paths:
                raise FileNotFoundError(f"no file matches {s}")
        else:
            if not os.path.exists(s):
                raise FileNotFoundError(s)
            paths = [s]
        return paths, [None] * len(paths)
    paths = [str(p) for p in source]
    if not paths:
        raise ValueError("no files given")
    return paths, [None] * len(paths)


def _samples_frame(samples, index):
    import pandas as pd
    if isinstance(samples, (str, os.PathLike)):
        p = str(samples)
        sep = "," if p.lower().removesuffix(".gz").endswith(".csv") else "\t"
        samples = pd.read_csv(p, sep=sep, index_col=0)
    if not isinstance(samples, pd.DataFrame):
        raise TypeError("samples must be a DataFrame or a TSV/CSV path (first column = sample name)")
    samples = samples.copy()
    samples.index = samples.index.astype(str)
    missing = [s for s in index if s not in samples.index]
    if missing:
        raise KeyError(f"samples table lacks {missing[:5]}; its index is {samples.index.tolist()[:5]}...")
    return samples.loc[list(index)]


def _row_table(rows, level):
    from .genes import Genes
    if isinstance(rows, Genes):
        return rows.transcripts if level == "transcript" else rows.genes
    return rows


def _row_ids(rows, level):
    from .genes import Genes
    if isinstance(rows, Genes):
        return rows.transcripts.cols["transcript_id"] if level == "transcript" else rows.genes.cols["gene_id"]
    return rows.cols["name"] if "name" in rows.cols else rows.uid


def _row_names(rows, level):
    from .genes import Genes
    if isinstance(rows, Genes):
        if level == "transcript":
            return rows.genes.cols["gene_name"][rows.transcripts.cols["gene"]]
        return rows.genes.cols["gene_name"]
    return None


def _rows_label(rows, level=None) -> str:
    from .genes import Genes
    if isinstance(rows, Genes):
        return "genes.transcripts" if level == "transcript" else "genes.genes"
    return f"Loci(n={len(rows):,})"
