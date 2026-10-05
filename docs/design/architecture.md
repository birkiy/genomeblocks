---
title: Architecture
parent: Design
layout: default
nav_order: 7
---

# Architecture
{: .no_toc }

`Architecture` is a chromatin-contact graph: CREs are vertices, loops are
edges, and Hi-C contacts become edge weights. It subclasses graph-tool's
`Graph`, so every graph-tool algorithm runs on it directly.
{: .fs-5 .fw-300 }

1. TOC
{:toc}

## The pipeline
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/arch-pipeline.svg %}
</div><figcaption>
<strong>Six chained calls, each writing graph properties.</strong> Edge properties (<code>ep.w</code>, <code>ep.d</code>, <code>ep.n</code>) and vertex properties (<code>vp.annot</code>, <code>vp.gene</code>, <code>vp.strength</code>) are graph-tool property maps. Each step reads what the previous one wrote, so the steps run in this order.
</figcaption></figure>

Vertices are keyed by CRE uid: `vp.uid` holds it and `G.index` maps a uid to
its vertex. `G[uid]` returns that CRE's neighbours with their weights, and
`G[uid1, uid2]` returns the edge's properties.

## From loops to weighted edges
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/arch-weights.svg %}
</div><figcaption>
<strong>Loops make edges, Hi-C weighs them, distance normalises them.</strong> <code>make</code> takes every CRE within ± <code>r</code> of each loop anchor's midpoint and links every CRE under one anchor to every CRE under the other. <code>add_mcool</code> maps each CRE to its nearest cooler bin and gives each edge the pixel count of its bin pair, divided among the edges that share that pixel. <code>normalize</code> divides by the expected contact at that distance.
</figcaption></figure>

- **`make(loci, bedpe, r=2500, dmax=1e9)`** reads loops with
  [`read_bedpe`]({{ '/design/bedpe/' | relative_url }}), looks anchors up in the CREs' interval index, and
  adds each CRE pair once (self-pairs are skipped). Loops longer than `dmax`
  are ignored.
- **`add_mcool(loci, mcool, resolution=5000)`** joins every edge's
  `(bin1, bin2)` key against the cooler's sorted pixel table in one pass. A
  pixel's count is split evenly across the edges in it, so dense clusters of
  CREs do not multiply one contact. Edges without a pixel keep `w = 0`.
- **`normalize(loci)`** fits `w ≈ C·d^−α` on cis edges with positive distance and
  weight, and writes `ep.n = w / E(d)` and the distance `ep.d`. Trans edges have
  no distance: `ep.d = ∞` and their expectation is the mean trans weight.
  Co-located CREs (`d = 0`) get O/E 0; `prune()` removes them, once, at the
  end of a pipeline.

## Genes and hubs
{: .sec-purple }

<figure class="gb-fig"><div class="gb-fig-body">
{% include diagrams/arch-hubs.svg %}
</div><figcaption>
<strong>Contacts assign genes; strength finds hubs.</strong> <code>annotate</code> labels each CRE with its region class (see <a href="{{ '/design/genes/' | relative_url }}">Genes</a>). A promoter CRE takes its nearest gene. Every other CRE takes the gene of the promoter it contacts most strongly, by the edge property <code>key</code> (O/E by default). <code>strength</code> sums each vertex's incident O/E. <code>prime_hubs</code> ranks CREs by that strength, cuts at the point where the scaled curve's slope reaches 1, and returns the hubs' genes.
</figcaption></figure>

- Gene assignment is vectorised over the edge array: one sort by
  (vertex, −weight) and the first row per vertex wins. Use a different `name=`
  per weight key (`key="n_CM", name="gene_CM"`) to keep several assignments
  side by side.
- `elbow(key)` measures the slope on a smoothed copy of the curve, so the cut
  follows the curve's shape rather than single-step noise.

## Storing and combining graphs
{: .sec-purple }

- **Pickling** stores vertex uids, edges as uid pairs, and every property map
  with its exact graph-tool value type, so `bool` and `int` maps survive a
  round trip.
- **`G | H`, `G & H`** combine graphs by uid; **`subgraph(...)`** filters by uid,
  property value or a function; **`copy()`** is deep.
- **`architecture_draw.draw(arch, loci, region)`** lays out the subgraph of CREs
  in a region as a network, with vertex size and edge width taken from
  properties.

The [columnar Architecture]({{ '/design/columnar/' | relative_url }}) keeps the same pipeline but stores edges
as one sorted table, which makes per-chromosome views free and handles trans
loops throughout. Measured numbers are on the
[Benchmarks]({{ '/benchmarks/#architecture' | relative_url }}) page.
