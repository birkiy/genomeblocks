"""Sequence-logo visualization for motifs / archetypes (logomaker).

Split out of ``motifs.py`` so motif scanning/clustering stays free of plotting
deps. matplotlib, logomaker and scipy are lazy-imported inside each function,
so importing this module is cheap. These are plain functions over PFM/PWM
arrays — nothing here is attached to Loci.
"""

def plot_archetype(pfm, ax=None, *, title=None, alphabet: str = 'ACGT',
                   show_xticks: bool = True, ylim=(0, 2)):
    """Plot a (4, W) PFM as an information-content sequence logo (bits)."""
    import logomaker
    import matplotlib.pyplot as plt
    import pandas as pd

    df = pd.DataFrame(pfm.T, columns=list(alphabet))
    ic_df = logomaker.transform_matrix(
        df, from_type='probability', to_type='information')

    if ax is None:
        _, ax = plt.subplots(figsize=(max(2.5, ic_df.shape[0] * 0.35), 1.4))

    logo = logomaker.Logo(ic_df, ax=ax, show_spines=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    ax.set_ylabel('bits', fontsize=8)
    ax.set_ylim(*ylim)
    if show_xticks:
        ax.set_xticks(range(ic_df.shape[0]))
        ax.set_xticklabels(range(1, ic_df.shape[0] + 1), fontsize=7)
    if title:
        ax.set_title(title, fontsize=9)
    return ax


def plot_archetypes(archetypes, *, ncols: int = 4, figsize_per=(2.8, 1.3),
                    members=None, sort_by_size: bool = True):
    """Grid of archetype logos. Pass ``members=`` to annotate cluster size."""
    import matplotlib.pyplot as plt
    import math

    items = list(archetypes.items())
    if sort_by_size and members is not None:
        items.sort(key=lambda kv: -len(members.get(kv[0], [])))
    n = len(items)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize_per[0] * ncols, figsize_per[1] * nrows),
        squeeze=False,
    )
    flat = axes.flatten()
    for ax, (name, pfm) in zip(flat, items):
        title = name
        if members is not None and name in members:
            title = f"{name}  (n={len(members[name])})"
        plot_archetype(pfm, ax=ax, title=title)
    for ax in flat[n:]:
        ax.axis('off')
    plt.tight_layout()
    return fig


def plot_cluster_members(archetype_pfm, member_pwms, member_names, *,
                         ncols: int = 3, figsize_per=(2.8, 1.3),
                         archetype_title: str = 'ARCHETYPE'):
    """Archetype on top + each member underneath for QC."""
    import matplotlib.pyplot as plt
    import math

    n = len(member_pwms) + 1
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize_per[0] * ncols, figsize_per[1] * nrows),
        squeeze=False,
    )
    flat = axes.flatten()
    plot_archetype(archetype_pfm, ax=flat[0], title=archetype_title)
    for ax, pfm, name in zip(flat[1:], member_pwms, member_names):
        plot_archetype(pfm, ax=ax, title=name)
    for ax in flat[n:]:
        ax.axis('off')
    plt.tight_layout()
    return fig


def plot_dendrogram(Z, names=None, *, cutoff=None, ax=None,
                    color_threshold=None, leaf_font_size: int = 6,
                    no_labels: bool = False):
    """Dendrogram from a linkage matrix Z (output of ``cluster_motifs``).
    Draws a horizontal line at ``cutoff`` if provided."""
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram
    if ax is None:
        _, ax = plt.subplots(figsize=(14, 4))
    color_threshold = color_threshold if color_threshold is not None else cutoff
    dendrogram(
        Z, labels=names, ax=ax,
        color_threshold=color_threshold,
        leaf_font_size=leaf_font_size,
        leaf_rotation=90,
        no_labels=no_labels,
    )
    if cutoff is not None:
        ax.axhline(cutoff, color='red', ls='--', lw=0.8, alpha=0.6,
                   label=f'cutoff={cutoff}')
        ax.legend(loc='upper right', fontsize=7)
    ax.set_ylabel('SW distance', fontsize=9)
    return ax
