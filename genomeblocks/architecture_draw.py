"""Region drawing for :class:`~genomeblocks.architecture.Architecture` graphs.

Kept separate from ``architecture.py`` so the core graph object stays lightweight
and free of matplotlib / drawing dependencies. Draw a genomic region with::

    from genomeblocks.draw import draw
    draw(arch, loci, region=("chr1", 1_000_000, 2_000_000))

``arch`` is an :class:`Architecture`, ``loci`` the same ``Loci`` used to build it.
"""
import graph_tool.all as gt


def _merge_nearby(arch, sub_loci, merge_distance, edge_key='w', label_prop='gene',
                  vertex_size_by=None, color_prop=None):
    """Merge loci within *merge_distance* bp into single nodes.

    Builds a fresh gt.Graph where each cluster of nearby loci becomes one
    vertex, and inter-cluster edges carry the summed weight of member edges.

    Args:
        sub_loci: Loci in the region that are also present in the graph.
        merge_distance: Max center-to-center bp for two consecutive sorted
            loci to be merged (uses Locus.distance_to).
        edge_key: Edge property to aggregate (default 'w').
        label_prop: Vertex property used for labels (default 'gene').
        vertex_size_by: Vertex property to sum for merged node sizes, or None.
        color_prop: Optional vertex PropertyMap on arch whose values should be
            aggregated per cluster into ``mg.vp.color_val``. String values take
            the first non-empty member; vector values average per-component;
            scalars average.

    Returns:
        (merged_graph, clusters, uid_to_cluster)
        merged_graph  – gt.Graph with vp.uid, vp.label, vp.size_val, ep.w
                        and (if color_prop given) vp.color_val
        clusters      – list[list[Locus]]
        uid_to_cluster – dict mapping original uid → cluster index
    """
    from collections import defaultdict
    import numpy as np

    graph_loci = sorted(
        [l for l in sub_loci if l.uid in arch.index],
        key=lambda l: (l.chrom, l.start),
    )

    # Greedy clustering on sorted loci
    clusters = [[graph_loci[0]]]
    for l in graph_loci[1:]:
        if clusters[-1][-1].distance_to(l) <= merge_distance:
            clusters[-1].append(l)
        else:
            clusters.append([l])

    uid_to_cluster = {}
    for ci, cluster in enumerate(clusters):
        for loc in cluster:
            uid_to_cluster[loc.uid] = ci

    # Build merged graph
    mg = gt.Graph(directed=False)
    mg.vp.uid   = mg.new_vertex_property("string")
    mg.vp.label = mg.new_vertex_property("string")
    mg.vp.size_val = mg.new_vertex_property("float")
    mg.ep.w     = mg.new_edge_property("float")

    color_vtype = color_prop.value_type() if color_prop is not None else None
    if color_prop is not None:
        mg.vp.color_val = mg.new_vertex_property(color_vtype)

    cluster_verts = []
    for ci, cluster in enumerate(clusters):
        v = mg.add_vertex()
        cluster_verts.append(v)

        first, last = cluster[0], cluster[-1]
        mg.vp.uid[v] = f"{first.chrom}:{first.start}-{last.end}"

        # Label
        if label_prop in arch.vp:
            labels = []
            for loc in cluster:
                if loc.uid in arch.index:
                    lbl = str(arch.vp[label_prop][arch.index[loc.uid]])
                    if lbl and lbl not in labels:
                        labels.append(lbl)
            mg.vp.label[v] = "/".join(labels) if labels else mg.vp.uid[v]
        else:
            mg.vp.label[v] = mg.vp.uid[v]

        # Size value
        if vertex_size_by and vertex_size_by in arch.vp:
            mg.vp.size_val[v] = sum(
                arch.vp[vertex_size_by][arch.index[loc.uid]]
                for loc in cluster if loc.uid in arch.index
            )
        else:
            mg.vp.size_val[v] = float(len(cluster))

        # Color value
        if color_prop is not None:
            vals = [color_prop[arch.index[loc.uid]]
                    for loc in cluster if loc.uid in arch.index]
            if vals:
                if color_vtype == "string":
                    mg.vp.color_val[v] = next((x for x in vals if x), vals[0])
                elif "vector" in color_vtype:
                    arr = np.array([list(x) for x in vals], dtype=float)
                    mg.vp.color_val[v] = arr.mean(axis=0).tolist()
                else:
                    mg.vp.color_val[v] = float(np.mean(vals))

    # Aggregate inter-cluster edges
    ep_key = edge_key if edge_key and edge_key in arch.ep else 'w'
    edge_weights = defaultdict(float)
    for ci, cluster in enumerate(clusters):
        for loc in cluster:
            if loc.uid not in arch.index:
                continue
            v_orig = arch.index[loc.uid]
            for nb in v_orig.all_neighbors():
                n_uid = arch.vp.uid[nb]
                if n_uid not in uid_to_cluster:
                    continue
                cj = uid_to_cluster[n_uid]
                if ci == cj:
                    continue
                pair = (min(ci, cj), max(ci, cj))
                edge_weights[pair] += arch.ep[ep_key][arch.edge(v_orig, nb)]

    for (ci, cj), w in edge_weights.items():
        e = mg.add_edge(cluster_verts[ci], cluster_verts[cj])
        mg.ep.w[e] = w

    return mg, clusters, uid_to_cluster


def draw(arch, loci, region, *,
         merge_distance=None,
         vertex_size_by=None,
         edge_width_by='w',
         vertex_size_range=(10, 50),
         edge_width_range=(1, 10),
         vertex_color=None,
         edge_color='#CCCCCC',
         figsize=(12, 8),
         layout='spring',
         show_labels=True,
         label_prop='gene',
         label_position='outside',
         vertex_font_size=10,
         ax=None,
         **kwargs):
    """Draw the subgraph for CREs overlapping a genomic region.

    Args:
        loci: Loci collection (same one used to build the graph)
        region: Tuple of (chrom, start, end) or Locus object
        merge_distance: If set, merge loci within this many bp into single
            nodes before drawing (center-to-center via Locus.distance_to).
        vertex_size_by: Vertex property name to scale node sizes (e.g., 'Agg_H1')
        edge_width_by: Edge property name to scale edge widths (default: 'w')
        vertex_size_range: (min_size, max_size) for vertices (default: (10, 50))
        edge_width_range: (min_width, max_width) for edges (default: (1, 10))
        vertex_color: Color for vertices - can be:
            - String (e.g., '#4A90E2') - single color for all
            - Vertex property name (e.g., 'crest_color') - color by values
            - PropertyMap (e.g., arch.vp.crest_color) - color by values
            - None - defaults to blue
            Under ``merge_distance``, per-vertex colors are aggregated across
            each cluster (string → first non-empty; vector → mean per channel;
            scalar → mean).
        edge_color: Color for edges (default: gray)
        figsize: Figure size as (width, height)
        layout: Layout algorithm ('spring', 'circular', 'kamada_kawai')
        show_labels: Whether to show node labels
        label_prop: Vertex property to use for labels (default: 'gene')
        label_position: Where to render labels — 'outside' (default; preserves
            vertex sizes by placing text radially) or 'inside' (may inflate
            vertices to fit text). Ignored when ``show_labels=False``.
        vertex_font_size: Font size for vertex labels (default: 10)
        ax: Matplotlib axis to plot on (creates new if None)
        **kwargs: Additional arguments passed to graph_tool's graph_draw

    Returns:
        Matplotlib axis object
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from graph_tool.draw import prop_to_size

    # Parse region
    if hasattr(region, 'chrom'):  # Locus object
        chrom, start, end = region.chrom, region.start, region.end
    else:
        chrom, start, end = region

    # Get overlapping loci
    sub_loci = loci.overlaps(chrom, start, end)
    if len(sub_loci) == 0:
        print(f"[WARNING] No loci found in region {chrom}:{start}-{end}")
        return ax

    # Get UIDs in region
    region_uids = {l.uid for l in sub_loci}
    region_uids_in_graph = {uid for uid in region_uids if uid in arch}

    if len(region_uids_in_graph) == 0:
        print(f"[WARNING] No graph vertices found in region {chrom}:{start}-{end}")
        return ax

    # ── Resolve vertex_color once, shared across paths ───────────────
    # color_prop: a PropertyMap on arch (needs to be remapped under merge)
    # color_literal: a scalar color value (string / tuple) or None
    color_prop = None
    color_literal = None
    if vertex_color is None:
        color_literal = '#4A90E2'
    elif isinstance(vertex_color, str):
        if vertex_color in arch.vp:
            color_prop = arch.vp[vertex_color]
        else:
            color_literal = vertex_color
    elif isinstance(vertex_color, gt.PropertyMap):
        if vertex_color.get_graph().base is not arch.base:
            raise ValueError(
                "vertex_color PropertyMap is not bound to this Architecture"
            )
        color_prop = vertex_color
    else:
        color_literal = vertex_color  # tuple/list RGBA, etc.

    # ── Merge nearby loci if requested ───────────────────────────────
    if merge_distance is not None:
        graph_loci = [l for l in sub_loci if l.uid in region_uids_in_graph]
        if len(graph_loci) == 0:
            print(f"[WARNING] No graph loci to merge in region {chrom}:{start}-{end}")
            return ax

        n_before = len(graph_loci)
        mg, clusters, _ = _merge_nearby(arch, 
            sub_loci, merge_distance,
            edge_key=edge_width_by, label_prop=label_prop,
            vertex_size_by=vertex_size_by,
            color_prop=color_prop,
        )
        subgraph = mg

        # Vertex sizes from aggregated size_val
        if vertex_size_by and vertex_size_by in arch.vp:
            vals = mg.vp.size_val.a
            if vals.max() > vals.min():
                normed = (vals - vals.min()) / (vals.max() - vals.min())
                mg.vp.size_val.a = vertex_size_range[0] + normed * (vertex_size_range[1] - vertex_size_range[0])
            else:
                mg.vp.size_val.a[:] = np.mean(vertex_size_range)
            vertex_sizes = mg.vp.size_val
        else:
            vertex_sizes = np.mean(vertex_size_range)

        # Edge widths from aggregated ep.w
        if mg.num_edges() > 0:
            vals = mg.ep.w.a
            if vals.max() > vals.min():
                normed = (vals - vals.min()) / (vals.max() - vals.min())
                mg.ep.w.a = edge_width_range[0] + normed * (edge_width_range[1] - edge_width_range[0])
            else:
                mg.ep.w.a[:] = np.mean(edge_width_range)
            edge_widths = mg.ep.w
        else:
            edge_widths = np.mean(edge_width_range)

        vertex_fill_color = mg.vp.color_val if color_prop is not None else color_literal
        labels = mg.vp.label if show_labels else None

        title_extra = (f"{subgraph.num_vertices()} nodes "
                       f"({n_before} loci merged at {merge_distance}bp), "
                       f"{subgraph.num_edges()} edges")
    else:
        # ── Standard (unmerged) path ─────────────────────────────────
        vfilt = arch.new_vertex_property("bool")
        for v in arch.vertices():
            vfilt[v] = arch.vp.uid[v] in region_uids_in_graph

        subgraph = gt.GraphView(arch, vfilt=vfilt)

        if subgraph.num_vertices() == 0:
            print(f"[WARNING] No vertices in subgraph for region {chrom}:{start}-{end}")
            return ax

        # Vertex sizes
        if vertex_size_by and vertex_size_by in arch.vp:
            vertex_sizes = prop_to_size(
                arch.vp[vertex_size_by],
                mi=vertex_size_range[0],
                ma=vertex_size_range[1],
                power=1.0,
            )
        else:
            vertex_sizes = np.mean(vertex_size_range)

        # Edge widths
        if edge_width_by and edge_width_by in arch.ep:
            edge_widths = prop_to_size(
                arch.ep[edge_width_by],
                mi=edge_width_range[0],
                ma=edge_width_range[1],
                power=1.0,
            )
        else:
            edge_widths = np.mean(edge_width_range)

        vertex_fill_color = color_prop if color_prop is not None else color_literal

        # Labels
        if show_labels and label_prop in arch.vp:
            labels = subgraph.new_vertex_property("string")
            for v in subgraph.vertices():
                labels[v] = str(arch.vp[label_prop][v])
        else:
            labels = None

        title_extra = (f"{subgraph.num_vertices()} nodes, "
                       f"{subgraph.num_edges()} edges")

    # ── Common drawing code ──────────────────────────────────────────
    if layout == 'spring':
        pos = gt.sfdp_layout(subgraph)
    elif layout == 'circular':
        pos = gt.circular_layout(subgraph)
    elif layout == 'kamada_kawai':
        pos = gt.kamada_kawai_layout(subgraph)
    else:
        pos = gt.sfdp_layout(subgraph)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Label placement: -1 draws inside the vertex and forces it to grow to fit
    # the text (breaks size scaling); any non-negative value renders text
    # radially outside, preserving vertex_size.
    if labels is not None:
        if label_position == 'outside':
            kwargs.setdefault('vertex_text_position', 0)
        elif label_position == 'inside':
            kwargs.setdefault('vertex_text_position', -1)
        else:
            raise ValueError(
                f"label_position must be 'inside' or 'outside', got {label_position!r}"
            )

    gt.graph_draw(
        subgraph,
        pos=pos,
        vertex_size=vertex_sizes,
        vertex_fill_color=vertex_fill_color,
        edge_pen_width=edge_widths,
        edge_color=edge_color,
        vertex_text=labels,
        vertex_font_size=vertex_font_size,
        output_size=tuple(int(x * 100) for x in figsize),
        mplfig=ax,
        **kwargs,
    )

    ax.set_title(f"{arch._name}: {chrom}:{start:,}-{end:,}\n{title_extra}",
                 fontsize=14, pad=10)
    ax.axis('off')

    return ax
