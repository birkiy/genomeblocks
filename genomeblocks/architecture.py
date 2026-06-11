"""Architecture graph utilities (graph-tool wrapper)."""
from typing import Optional

import graph_tool.all as gt # ignore import


class Architecture(gt.Graph):

    def __init__(s, name: Optional[str]=None):
        super().__init__(directed=False)
        s._name = name or "Architecture"
        # --- Vertex Properties (for loci) ---
        s.vp.uid = s.new_vertex_property("string")
        s.index = {}

        # --- Edge Properties (for links) ---
        s.ep.w = s.new_edge_property("float")
        s.ep.n = s.new_edge_property("float")
        s.ep.d = s.new_edge_property("float")
        print("[INFO] Initialized an empty Architecture graph. 🏗️")

    def _add_vertex(s, uid):
        if uid not in s.index:
            v = s.add_vertex()
            s.vp.uid[v] = uid
            s.index[uid] = v
        return s.index[uid]

    @property
    def n_loci(s) -> int: return s.num_vertices()
    @property
    def n_links(s) -> int: return s.num_edges()

    def __str__(s) -> str:
        eprops = ", ".join(s.ep.keys())
        vprops = ", ".join(s.vp.keys())
        return f"Architecture(name='{s._name}', loci={s.n_loci}, links={s.n_links}, vertex_props=[{vprops}], edge_props=[{eprops}])"
    __repr__ = __str__

    def __len__(s) -> int: return s.num_vertices()
    def __contains__(s, uid: str) -> bool: return uid in s.index

    def __getitem__(s, key):
        if isinstance(key, str):
            if key not in s.index: raise KeyError(f"Locus UID '{key}' not found in the graph.")
            v = s.index[key]
            neighbors = {}
            for neighbor_v in v.all_neighbors():
                neighbor_uid = s.vp.uid[neighbor_v]
                edge = s.edge(v, neighbor_v)
                neighbors[neighbor_uid] = s.ep.w[edge]
            return neighbors
        if isinstance(key, tuple) and len(key) == 2:
            uid1, uid2 = key
            if uid1 not in s.index or uid2 not in s.index: raise KeyError(f"One or both UIDs ('{uid1}', '{uid2}') not found.")
            v1, v2 = s.index[uid1], s.index[uid2]
            edge = s.edge(v1, v2)
            if edge is None: return None
            return {name: prop[edge] for name, prop in s.ep.items()}
        raise TypeError("Key must be a Locus UID (string) or a tuple of two UIDs.")

    def copy(s) -> "Architecture":
        """Create a deep copy of the Architecture graph.
        
        Copies all vertices, edges, and ALL vertex and edge properties.
        
        Returns:
            New Architecture instance with all data copied
        """
        new_arch = s.__class__(name=s._name + "_copy")
        new_arch.set_directed(s.is_directed())
        
        # Create all vertex property maps in the new graph
        for vp_name, vp_map in s.vp.items():
            if vp_name == 'uid':
                continue  # uid is created by default
            # Infer property type from the first value
            sample_v = next(s.vertices(), None)
            if sample_v is not None:
                sample_val = vp_map[sample_v]
                if isinstance(sample_val, str):
                    new_arch.vp[vp_name] = new_arch.new_vertex_property("string")
                elif isinstance(sample_val, (int, float)):
                    new_arch.vp[vp_name] = new_arch.new_vertex_property("double")
                elif isinstance(sample_val, bool):
                    new_arch.vp[vp_name] = new_arch.new_vertex_property("bool")
                else:
                    new_arch.vp[vp_name] = new_arch.new_vertex_property("object")
        
        # Create all edge property maps in the new graph
        for ep_name, ep_map in s.ep.items():
            if ep_name in ['w', 'n', 'd']:
                continue  # These are created by default
            sample_e = next(s.edges(), None)
            if sample_e is not None:
                sample_val = ep_map[sample_e]
                if isinstance(sample_val, str):
                    new_arch.ep[ep_name] = new_arch.new_edge_property("string")
                elif isinstance(sample_val, (int, float)):
                    new_arch.ep[ep_name] = new_arch.new_edge_property("double")
                elif isinstance(sample_val, bool):
                    new_arch.ep[ep_name] = new_arch.new_edge_property("bool")
                else:
                    new_arch.ep[ep_name] = new_arch.new_edge_property("object")
        
        # Add vertices and copy all vertex properties
        for v_old in s.vertices():
            uid = s.vp.uid[v_old]
            v_new = new_arch._add_vertex(uid)
            
            # Copy all vertex properties
            for vp_name, vp_map in s.vp.items():
                if vp_name == 'uid':
                    continue  # Already set by _add_vertex
                new_arch.vp[vp_name][v_new] = vp_map[v_old]
        
        # Add edges and copy all edge properties
        for e_old in s.edges():
            uid1 = s.vp.uid[e_old.source()]
            uid2 = s.vp.uid[e_old.target()]
            v1_new = new_arch.index[uid1]
            v2_new = new_arch.index[uid2]
            e_new = new_arch.add_edge(v1_new, v2_new)
            
            # Copy all edge properties
            for ep_name, ep_map in s.ep.items():
                new_arch.ep[ep_name][e_new] = ep_map[e_old]
        
        print(f"[INFO] Created a copy of '{s._name}' with all properties. 🐑")
        return new_arch
    
    def subgraph(s, filter_func=None, vp_name=None, vp_values=None, uids=None, name=None) -> "Architecture":
        """Create a subgraph by filtering vertices.
        
        Multiple filtering options (use only one):
        - filter_func: Custom function that takes vertex and returns bool
        - vp_name + vp_values: Filter by vertex property values
        - uids: Filter by specific UIDs
        
        Args:
            filter_func: Function(vertex) -> bool to select vertices
            vp_name: Name of vertex property to filter by
            vp_values: Value or list of values to match for vp_name
            uids: Single UID or list of UIDs to include
            name: Name for the new subgraph (default: original_name_sub)
        
        Returns:
            New Architecture instance with filtered vertices and their edges
            
        Examples:
            >>> # Filter by UIDs
            >>> sub = G.subgraph(uids=['chr1:100000', 'chr2:200000'])
            
            >>> # Filter by vertex property
            >>> sub = G.subgraph(vp_name='gene', vp_values=['MYC', 'TP53'])
            
            >>> # Filter by custom function
            >>> sub = G.subgraph(filter_func=lambda v: G.vp.score[v] > 0.5)
        """
        if name is None:
            name = f"{s._name}_sub"
        
        new_arch = s.__class__(name=name)
        new_arch.set_directed(s.is_directed())
        
        # Determine which vertices to include
        include_vertices = set()
        
        if uids is not None:
            # Filter by UIDs
            if isinstance(uids, str):
                uids = [uids]
            for v in s.vertices():
                if s.vp.uid[v] in uids:
                    include_vertices.add(v)
        
        elif vp_name is not None and vp_values is not None:
            # Filter by vertex property
            if vp_name not in s.vp:
                raise ValueError(f"Vertex property '{vp_name}' not found in graph")
            
            if not isinstance(vp_values, (list, set, tuple)):
                vp_values = [vp_values]
            
            for v in s.vertices():
                if s.vp[vp_name][v] in vp_values:
                    include_vertices.add(v)
        
        elif filter_func is not None:
            # Filter by custom function
            for v in s.vertices():
                if filter_func(v):
                    include_vertices.add(v)
        
        else:
            raise ValueError("Must provide one of: filter_func, vp_name+vp_values, or uids")
        
        if len(include_vertices) == 0:
            print(f"[WARNING] No vertices matched filter criteria")
            return new_arch
        
        # Create all vertex property maps in the new graph
        for vp_name_iter, vp_map in s.vp.items():
            if vp_name_iter == 'uid':
                continue
            sample_v = next(iter(include_vertices))
            sample_val = vp_map[sample_v]
            if isinstance(sample_val, str):
                new_arch.vp[vp_name_iter] = new_arch.new_vertex_property("string")
            elif isinstance(sample_val, (int, float)):
                new_arch.vp[vp_name_iter] = new_arch.new_vertex_property("double")
            elif isinstance(sample_val, bool):
                new_arch.vp[vp_name_iter] = new_arch.new_vertex_property("bool")
            else:
                new_arch.vp[vp_name_iter] = new_arch.new_vertex_property("object")
        
        # Create all edge property maps
        for ep_name, ep_map in s.ep.items():
            if ep_name in ['w', 'n', 'd']:
                continue
            sample_e = next(s.edges(), None)
            if sample_e is not None:
                sample_val = ep_map[sample_e]
                if isinstance(sample_val, str):
                    new_arch.ep[ep_name] = new_arch.new_edge_property("string")
                elif isinstance(sample_val, (int, float)):
                    new_arch.ep[ep_name] = new_arch.new_edge_property("double")
                elif isinstance(sample_val, bool):
                    new_arch.ep[ep_name] = new_arch.new_edge_property("bool")
                else:
                    new_arch.ep[ep_name] = new_arch.new_edge_property("object")
        
        # Add filtered vertices with all properties
        for v_old in include_vertices:
            uid = s.vp.uid[v_old]
            v_new = new_arch._add_vertex(uid)
            
            # Copy all vertex properties
            for vp_name_iter, vp_map in s.vp.items():
                if vp_name_iter == 'uid':
                    continue
                new_arch.vp[vp_name_iter][v_new] = vp_map[v_old]
        
        # Add edges between included vertices
        edge_count = 0
        for e_old in s.edges():
            source_old = e_old.source()
            target_old = e_old.target()
            
            # Only add edge if both endpoints are in the subgraph
            if source_old in include_vertices and target_old in include_vertices:
                uid1 = s.vp.uid[source_old]
                uid2 = s.vp.uid[target_old]
                v1_new = new_arch.index[uid1]
                v2_new = new_arch.index[uid2]
                e_new = new_arch.add_edge(v1_new, v2_new)
                
                # Copy all edge properties
                for ep_name, ep_map in s.ep.items():
                    new_arch.ep[ep_name][e_new] = ep_map[e_old]
                
                edge_count += 1
        
        print(f"[INFO] Created subgraph '{name}': {new_arch.n_loci} vertices, {new_arch.n_links} edges "
              f"(from {s.n_loci} vertices, {s.n_links} edges)")
        return new_arch
    
    def __or__(s, other: "Architecture") -> "Architecture":
        """Union of two Architecture graphs (vertices and edges).
        
        Creates a new graph containing all vertices and edges from both graphs.
        Edge properties are taken from the first graph when edges exist in both.
        
        Args:
            other: Another Architecture instance
            
        Returns:
            New Architecture with union of vertices and edges
        """
        if not isinstance(other, Architecture):
            raise TypeError("Union requires another Architecture instance.")
        
        new_arch = s.__class__(name=f"{s._name}|{other._name}")
        
        # Add all vertices from both graphs
        all_uids = set()
        for v in s.vertices():
            uid = s.vp.uid[v]
            all_uids.add(uid)
            new_arch._add_vertex(uid)
        
        for v in other.vertices():
            uid = other.vp.uid[v]
            if uid not in all_uids:
                all_uids.add(uid)
                new_arch._add_vertex(uid)
        
        # Track edges to avoid duplicates
        edge_set = set()
        
        # Add edges from first graph
        for e in s.edges():
            uid1 = s.vp.uid[e.source()]
            uid2 = s.vp.uid[e.target()]
            edge_key = tuple(sorted([uid1, uid2]))
            
            if edge_key not in edge_set:
                v1 = new_arch.index[uid1]
                v2 = new_arch.index[uid2]
                e_new = new_arch.add_edge(v1, v2)
                
                # Copy edge properties from first graph
                for key, prop_map in s.ep.items():
                    if key not in new_arch.ep:
                        new_arch.ep[key] = new_arch.new_edge_property("float")
                    new_arch.ep[key][e_new] = prop_map[e]
                
                edge_set.add(edge_key)
        
        # Add edges from second graph (only if not already added)
        for e in other.edges():
            uid1 = other.vp.uid[e.source()]
            uid2 = other.vp.uid[e.target()]
            edge_key = tuple(sorted([uid1, uid2]))
            
            if edge_key not in edge_set:
                v1 = new_arch.index[uid1]
                v2 = new_arch.index[uid2]
                e_new = new_arch.add_edge(v1, v2)
                
                # Copy edge properties from second graph
                for key, prop_map in other.ep.items():
                    if key not in new_arch.ep:
                        new_arch.ep[key] = new_arch.new_edge_property("float")
                    new_arch.ep[key][e_new] = prop_map[e]
                
                edge_set.add(edge_key)
        
        print(f"[INFO] Union: {new_arch.n_loci} vertices, {new_arch.n_links} edges "
              f"(from {s.n_loci}+{other.n_loci} vertices, {s.n_links}+{other.n_links} edges)")
        return new_arch
    
    def __and__(s, other: "Architecture") -> "Architecture":
        """Intersection of two Architecture graphs (vertices and edges).
        
        Creates a new graph containing only vertices present in both graphs,
        and only edges that exist in both graphs.
        
        Args:
            other: Another Architecture instance
            
        Returns:
            New Architecture with intersection of vertices and edges
        """
        if not isinstance(other, Architecture):
            raise TypeError("Intersection requires another Architecture instance.")
        
        new_arch = s.__class__(name=f"{s._name}&{other._name}")
        
        # Find common vertices
        uids_s = {s.vp.uid[v] for v in s.vertices()}
        uids_other = {other.vp.uid[v] for v in other.vertices()}
        common_uids = uids_s & uids_other
        
        # Add common vertices
        for uid in common_uids:
            new_arch._add_vertex(uid)
        
        # Find edges in first graph with both endpoints in common vertices
        edges_s = set()
        for e in s.edges():
            uid1 = s.vp.uid[e.source()]
            uid2 = s.vp.uid[e.target()]
            if uid1 in common_uids and uid2 in common_uids:
                edges_s.add(tuple(sorted([uid1, uid2])))
        
        # Find edges in second graph with both endpoints in common vertices
        edges_other = set()
        edge_map_other = {}
        for e in other.edges():
            uid1 = other.vp.uid[e.source()]
            uid2 = other.vp.uid[e.target()]
            if uid1 in common_uids and uid2 in common_uids:
                edge_key = tuple(sorted([uid1, uid2]))
                edges_other.add(edge_key)
                edge_map_other[edge_key] = e
        
        # Find common edges and add them
        common_edges = edges_s & edges_other
        edge_map_s = {}
        for e in s.edges():
            uid1 = s.vp.uid[e.source()]
            uid2 = s.vp.uid[e.target()]
            edge_key = tuple(sorted([uid1, uid2]))
            if edge_key in common_edges:
                edge_map_s[edge_key] = e
        
        for edge_key in common_edges:
            uid1, uid2 = edge_key
            v1 = new_arch.index[uid1]
            v2 = new_arch.index[uid2]
            e_new = new_arch.add_edge(v1, v2)
            
            # Copy edge properties from first graph
            e_s = edge_map_s[edge_key]
            for key, prop_map in s.ep.items():
                if key not in new_arch.ep:
                    new_arch.ep[key] = new_arch.new_edge_property("float")
                new_arch.ep[key][e_new] = prop_map[e_s]
        
        print(f"[INFO] Intersection: {new_arch.n_loci} vertices, {new_arch.n_links} edges "
              f"(from {s.n_loci}∩{other.n_loci} vertices, {s.n_links}∩{other.n_links} edges)")
        return new_arch
    
    # lets add getstate and setstate for pickling
    def __getstate__(s):
        """Serialize Architecture to a pickleable dict.
        
        We can't pickle graph-tool's C++ internals directly, so we extract:
        - vertex UIDs in order
        - edges as (uid1, uid2) pairs
        - all vertex property maps as {name: [values...]}
        - all edge property maps as {name: [values...]}
        """
        verts = list(s.vertices())
        edges = list(s.edges())
        
        # Extract vertex UIDs
        vertex_uids = [s.vp.uid[v] for v in verts]
        
        # Extract edges as UID pairs
        edge_pairs = [(s.vp.uid[e.source()], s.vp.uid[e.target()]) for e in edges]
        
        # Extract vertex properties
        vprops = {}
        for name, prop in s.vp.items():
            vprops[name] = [prop[v] for v in verts]
        
        # Extract edge properties
        eprops = {}
        for name, prop in s.ep.items():
            eprops[name] = [prop[e] for e in edges]
        
        return {
            '_name': s._name,
            'directed': s.is_directed(),
            'vertex_uids': vertex_uids,
            'edge_pairs': edge_pairs,
            'vprops': vprops,
            'eprops': eprops,
        }

    def __setstate__(s, state):
        """Reconstruct Architecture from pickled state.
        
        We rebuild the graph from scratch:
        1. Initialize a new graph-tool Graph
        2. Add all vertices and recreate vertex properties
        3. Add all edges and recreate edge properties
        4. Rebuild the index mapping
        """
        # Initialize the graph-tool base class
        gt.Graph.__init__(s, directed=state.get('directed', False))
        
        # Set the name
        s._name = state.get('_name', 'Architecture')
        
        # Create vertex property maps first
        vprops_data = state.get('vprops', {})
        for name in vprops_data.keys():
            # Infer type from first non-None value
            values = vprops_data[name]
            sample_val = next((v for v in values if v is not None), None)
            if sample_val is None or isinstance(sample_val, str):
                prop_type = "string"
            elif isinstance(sample_val, (int, float)):
                prop_type = "float"
            else:
                prop_type = "string"
            s.vp[name] = s.new_vertex_property(prop_type)
        
        # Add vertices and set their properties
        vertex_uids = state.get('vertex_uids', [])
        vertex_map = []
        for i, uid in enumerate(vertex_uids):
            v = s.add_vertex()
            vertex_map.append(v)
            # Set all vertex properties for this vertex
            for name, values in vprops_data.items():
                s.vp[name][v] = values[i]
        
        # Build UID to vertex mapping for edge creation
        uid_to_vertex = {s.vp.uid[v]: v for v in vertex_map}
        
        # Create edge property maps
        eprops_data = state.get('eprops', {})
        for name in eprops_data.keys():
            values = eprops_data[name]
            sample_val = next((v for v in values if v is not None), None)
            if sample_val is None or isinstance(sample_val, str):
                prop_type = "string"
            elif isinstance(sample_val, (int, float)):
                prop_type = "float"
            else:
                prop_type = "string"
            s.ep[name] = s.new_edge_property(prop_type)
        
        # Add edges and set their properties
        edge_pairs = state.get('edge_pairs', [])
        for i, (uid1, uid2) in enumerate(edge_pairs):
            v1 = uid_to_vertex[uid1]
            v2 = uid_to_vertex[uid2]
            e = s.add_edge(v1, v2)
            # Set all edge properties for this edge
            for name, values in eprops_data.items():
                s.ep[name][e] = values[i]
        
        # Rebuild the index for fast UID lookup
        s.index = {s.vp.uid[v]: v for v in s.vertices()}


    def to_frame(self):
        import pandas as pd
        return pd.DataFrame({k: self.vp[k] for k in self.vp})


# ── Construction & contact overlay ──────────────────────────────────────────

@classmethod
def make(cls, loci, bedpe: str, *, name: str="Skeleton", r: int=2500, dmax=1e9, verbose: bool=True):
    from tqdm import tqdm
    G = cls(name=name)
    mapped_loops = 0
    total_loops = 0
    edge_set = set()
    with open(bedpe) as f:
        for line in tqdm(f, desc='[INFO] Building Architecture from loops'):
            if line.startswith('#') or not line.strip(): continue
            total_loops += 1
            fields = line.strip().split()
            if len(fields) < 6: continue
            try:
                chrom1, start1, end1 = fields[0], int(fields[1]), int(fields[2])
                chrom2, start2, end2 = fields[3], int(fields[4]), int(fields[5])
            except ValueError:
                continue
            mid1 = (start1 + end1) // 2; mid2 = (start2 + end2) // 2
            if abs(mid1 - mid2) > dmax: continue
            a1 = [loci[j] for *_, j in loci.cgr.overlap(chrom1, mid1 - r, mid1 + r)]
            a2 = [loci[j] for *_, j in loci.cgr.overlap(chrom2, mid2 - r, mid2 + r)]
            if not a1 or not a2: continue
            mapped_loops += 1
            for locus1 in a1:
                for locus2 in a2:
                    if locus1.uid == locus2.uid: continue
                    edge_key = tuple(sorted([locus1.uid, locus2.uid]))
                    if edge_key in edge_set: continue
                    edge_set.add(edge_key)
                    v1 = G._add_vertex(locus1.uid)
                    v2 = G._add_vertex(locus2.uid)
                    G.add_edge(v1, v2)
    if verbose:
        pct_mapped = 100 * mapped_loops / max(total_loops, 1)
        print(f"[INFO] {total_loops} loops | {mapped_loops} mapped ({pct_mapped:.1f}%) | "
              f"loci={G.n_loci}, links={G.n_links}")
    return G


def add_mcool(s: Architecture, loci, mcool: str, *, resolution: Optional[int]=None, name: str="w", verbose: bool=True) -> Architecture:
    import cooler
    from tqdm import tqdm
    import numpy as np
    from collections import defaultdict
    print(f"[INFO] Adding '{name}' weights to the graph. 🏗️")
    uri = f"{mcool}::resolutions/{resolution}" if resolution else mcool
    clr = cooler.Cooler(uri)
    bins = clr.bins()[:][['chrom', 'start', 'end']].reset_index()
    from .loci import Loci
    from .locus import Locus
    bins_l = Loci(Locus(row[1], row[2], row[3]) for row in bins.itertuples(index=False))
    near = loci.nearest(bins_l)
    uid_to_bin = dict(zip(near['Name'], near['Name_b'].map(bins_l.uids)))
    bin_pair_edge_counts = defaultdict(int)
    for edge in s.edges():
        v1, v2 = edge.source(), edge.target()
        uid1, uid2 = s.vp.uid[v1], s.vp.uid[v2]
        bin1 = uid_to_bin.get(uid1)
        bin2 = uid_to_bin.get(uid2)
        if bin1 is None or bin2 is None: continue
        if bin1 > bin2: bin1, bin2 = bin2, bin1
        bin_pair_edge_counts[(bin1, bin2)] += 1
    pixels = clr.pixels()[:].set_index(['bin1_id', 'bin2_id'])['count']
    edges_set = 0
    s.ep[name] = s.new_edge_property("float")
    for edge in tqdm(s.edges(), desc='[INFO] Assigning distributed weights to edges', total=s.n_links):
        v1, v2 = edge.source(), edge.target()
        uid1, uid2 = s.vp.uid[v1], s.vp.uid[v2]
        bin1 = uid_to_bin.get(uid1)
        bin2 = uid_to_bin.get(uid2)
        if bin1 is None or bin2 is None: continue
        bin1, bin2 = (min(bin1, bin2), max(bin1, bin2))
        try:
            total_count = float(pixels[bin1][bin2].sum())
            if total_count > 0:
                num_edges_in_pair = bin_pair_edge_counts.get((bin1, bin2), 1)
                distributed_weight = total_count / num_edges_in_pair
                s.ep[name][edge] += distributed_weight
                edges_set += 1
        except KeyError:
            s.ep[name][edge] += 0
        except TypeError:
            print(pixels[bin1][bin2])
            break
    if verbose:
        print(f"[INFO] Set distributed weights for {edges_set}/{s.n_links} edges from cooler. [{name}]")
    return s


# ── Distance-decay normalization ────────────────────────────────────────────

def _pl_model(x, C, alpha):
    import numpy as np
    x = np.asarray(x, float)
    out = np.full(x.shape, np.inf)          # x<=0 is singular → inf (→ O/E 0)
    pos = x > 0
    out[pos] = C * (x[pos] ** (-alpha))
    return out


def _pl_expect(x, y):
    import numpy as np
    from scipy.optimize import curve_fit
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    if not m.any():
        return np.full_like(y, np.nan), {"alpha": np.nan, "C": np.nan}
    try:
        b, a = np.polyfit(np.log10(x[m]), np.log10(y[m]), 1)
        guess = [10.0**a, -b]
        popt, _ = curve_fit(_pl_model, x[m], y[m], p0=guess)
        C_fit, alpha_fit = popt
        y_expected = _pl_model(x, C_fit, alpha_fit)
        return y_expected, {"alpha": float(alpha_fit), "C": float(C_fit)}
    except (RuntimeError, np.linalg.LinAlgError, ValueError):
        return np.full_like(y, np.nan), {"alpha": np.nan, "C": np.nan}


def normalize(s: Architecture, loci, *, source: str = "w", name: str = "n", verbose: bool = True) -> Architecture:
    """Normalize edge weights by their distance-decay (power-law) expectation.

    Fits a power law to the (distance, ``ep[source]``) cloud across all edges,
    then stores the observed/expected ratio in ``ep[name]``. Edge distances are
    computed from ``loci`` and cached in ``ep.d``.

    Zero-distance edges (overlapping/co-located loci) are a power-law
    singularity: ``_pl_expect`` fits on positive distances only, and
    ``_pl_model`` returns inf at distance 0 so their O/E falls to 0. They are
    not removed here — call :meth:`prune` once at the end to drop them (removing
    edges mid-pipeline desyncs the edge-index range across repeated calls).
    """
    import numpy as np
    from tqdm import tqdm
    if verbose:
        print(f"[INFO] Normalizing '{source}' by power-law expectation. Storing in '{name}'. 📏")
    if source not in s.ep:
        raise ValueError(f"Source edge property '{source}' not found in the graph.")
    if name not in s.ep:
        s.ep[name] = s.new_edge_property("float")
    edge_list = list(s.edges())
    distances = np.zeros(len(edge_list), dtype=float)
    raw_weights = np.zeros(len(edge_list), dtype=float)
    for i, edge in tqdm(enumerate(edge_list), desc='[INFO] Calculating distances for edges', total=len(edge_list)):
        v1, v2 = edge.source(), edge.target()
        uid1, uid2 = s.vp.uid[v1], s.vp.uid[v2]
        dist = loci[uid1].distance_to(loci[uid2])
        s.ep.d[edge] = dist
        raw_weights[i] = s.ep[source][edge]
        distances[i] = dist
    expected, fit_params = _pl_expect(distances, raw_weights)
    if verbose:
        print(f"[INFO] Power-law fit complete: alpha={fit_params['alpha']:.3f}, C={fit_params['C']:.3e}")
    norm_weights = raw_weights / np.maximum(expected, 1e-12)
    for i, e in enumerate(s.edges()):
        s.ep[name][e] = norm_weights[i]
    if verbose:
        print(f"[INFO] Set O/E weights for {s.n_links} intra-chromosomal edges to `ep.{name}`.")
    return s


def prune(s: Architecture, *, dist_prop: str = "d", verbose: bool = True) -> Architecture:
    """Remove zero-distance edges (overlapping/co-located loci).

    These share a center (``distance_to == 0``) and are a power-law
    singularity, so they carry no meaningful O/E weight. Distances are read
    from ``ep[dist_prop]`` (populated by ``normalize``), so run ``normalize``
    before ``prune``.

    Call this **once, at the very end** of a multi-cell pipeline. Removing
    edges mid-pipeline leaves the edge-index range uncompacted, which desyncs
    ``new_edge_property().a`` from ``list(edges())`` on the next ``normalize``.
    """
    import numpy as np
    if dist_prop not in s.ep:
        raise ValueError(
            f"Edge property '{dist_prop}' not found — run normalize() first "
            f"(it populates ep.{dist_prop} with edge distances).")
    # Derive the mask from ep.d's array so it matches new_edge_property().a
    # sizing (both indexed by edge-index range) — mixing a list-of-edges mask
    # with a property array is exactly what broke before.
    keep = np.asarray(s.ep[dist_prop].a) > 0
    n_zero = int((~keep).sum())
    if n_zero == 0:
        if verbose:
            print("[INFO] prune: no zero-distance edges to remove.")
        return s
    keep_ep = s.new_edge_property("bool")
    keep_ep.a = keep
    s.set_edge_filter(keep_ep)
    s.purge_edges()
    s.set_edge_filter(None)
    if verbose:
        print(f"[INFO] prune: removed {n_zero} zero-distance edges → {s.n_links} edges.")
    return s


# ── Gene annotation ─────────────────────────────────────────────────────────

def annotate(self, loci, genes, *, key='n', name='gene', verbose=True):
    """Add gene and genomic region annotations to vertices.

    Two-stage gene assignment:
      1. Region label (``vp.annot``) and promoter genes. A CRE labelled
         'Promoter-TSS' gets ``vp[name]`` = its nearest gene.
      2. Every non-promoter CRE is assigned (``vp[name]``) to the gene of its
         **highest-weight promoter neighbour**, scored by edge property
         ``key``. CREs with no promoter contact stay unassigned ('').

    Stores per vertex:
        vp.annot   : region label (Promoter-TSS / 5UTR / 3UTR / Exonic /
                     Intronic / Intergenic).
        vp[name]   : promoters → nearest gene; other CREs → gene of the
                     top-``key``-weight promoter they contact ('' if none).
                     Use a distinct ``name`` per ``key`` (e.g. key='n_CM',
                     name='gene_CM') to keep assignments side by side.

    Requires edge property ``key`` (run add_mcool/normalize first) so the
    interaction-based assignment in stage 2 has weights to rank by.
    """
    import numpy as np

    if key not in self.ep:
        raise ValueError(
            f"Edge property '{key}' not found — compute edge weights before "
            f"annotate (e.g. add_mcool/normalize) so promoter assignment can "
            f"rank by interaction strength.")

    cre_a_df = genes.annotations(loci)
    annot_map = cre_a_df.set_index('uid')['annotation'].to_dict()

    cre_g_df = genes.nearest_genes(loci)
    nearest_map = dict(zip(cre_g_df['Name'], cre_g_df['Name_b']))

    self.vp[name] = self.new_vp('string')
    self.vp.annot = self.new_vp('string')
    gene_pm = self.vp[name]
    annot_pm = self.vp.annot

    # ── Stage 1: bulk-extract uids, derive labels/genes, single write pass ──
    # PropertyMap iteration is C-level vertex-index order — much faster than
    # per-vertex indexing.
    uids = list(self.vp.uid)
    annot_vals = [annot_map.get(u, '') for u in uids]
    is_prom = np.fromiter(('Promoter' in a for a in annot_vals),
                          dtype=bool, count=len(annot_vals))

    n_prom = 0
    for v, a, ip, u in zip(self.vertices(), annot_vals, is_prom, uids):
        annot_pm[v] = a
        if ip:
            gene_pm[v] = nearest_map.get(u, '')
            n_prom += 1

    # ── Stage 2: vectorized highest-weight promoter neighbour per non-prom ──
    # Old per-vertex/per-edge loop did ~10M PropertyMap lookups for a
    # 1.9M-edge graph. Replaced with numpy ops on the edge array.
    ew = self.get_edges(eprops=[self.ep[key]])     # (E, 3): src, tgt, weight
    src = ew[:, 0].astype(np.int64, copy=False)
    tgt = ew[:, 1].astype(np.int64, copy=False)
    w   = ew[:, 2].astype(float, copy=False)

    gene_arr = np.asarray(list(gene_pm), dtype=object)
    eligible = is_prom & (gene_arr != '')          # promoter w/ a gene

    # Each edge contributes at most one (non-prom-vertex, prom-gene, weight)
    # row; self-loops and prom-prom / non-prom-non-prom edges drop out.
    case1 = ~is_prom[src] & eligible[tgt]
    case2 = ~is_prom[tgt] & eligible[src]

    np_idx = np.concatenate([src[case1], tgt[case2]])
    pgene  = np.concatenate([gene_arr[tgt[case1]], gene_arr[src[case2]]])
    we     = np.concatenate([w[case1], w[case2]])

    n_assigned = 0
    if np_idx.size:
        # Sort by (non-prom-vertex asc, weight desc) — first row per vertex
        # group is the max-weight promoter contact.
        order = np.lexsort((-we, np_idx))
        np_sorted = np_idx[order]
        pgene_sorted = pgene[order]
        _, first = np.unique(np_sorted, return_index=True)
        for vi, gn in zip(np_sorted[first], pgene_sorted[first]):
            if gn:
                gene_pm[self.vertex(int(vi))] = gn
                n_assigned += 1

    if verbose:
        n_others = int((~is_prom).sum())
        print(f"[INFO] Annotated {self.n_loci} loci: {n_prom} promoter CREs | "
              f"{n_assigned}/{n_others} non-promoter CREs assigned to a "
              f"top-'{key}' promoter gene → vp.{name}.")
    return self


# ── Hubs (node strength → slope-1 knee) ─────────────────────────────────────

def strength(s: Architecture, key: str = "n", name: str = "strength", *, verbose: bool = True) -> Architecture:
    """Compute node strength (sum of incident edge weights) per vertex.

    For each vertex, sums the edge property ``key`` across all incident edges
    and stores the raw sum in ``vp[name]``. No normalization is applied — if you
    want a normalized version, divide ``vp[name]`` by its total yourself.

    Args:
        key:  Edge property to sum (default ``"n"`` = O/E weights).
        name: Vertex property name for the result (default ``"strength"``).

    Returns:
        self (for chaining).
    """
    if key not in s.ep:
        raise ValueError(f"Edge property '{key}' not found.")

    s.vp[name] = s.new_vp("float")
    for e in s.edges():
        w = s.ep[key][e]
        s.vp[name][e.source()] += w
        s.vp[name][e.target()] += w

    if verbose:
        print(f"[INFO] Summed ep.{key} → vp.{name} (node strength).")
    return s


def elbow(s: Architecture, key: str, *, verbose: bool = True):
    """Find a cutoff on a vertex property's sorted curve via the slope-1 knee.

    Sorts vertices by ``vp[key]`` descending (dropping zeros), then on the
    ascending curve with both axes normalized to ``[0, 1]`` finds the point
    where the local tangent slope equals ``1`` (the 45° knee). The slope is
    measured on a smoothed curve so it reflects the curve's overall shape
    rather than single-step jitter. Hubs are everything above the cutoff.

    Args:
        key: Vertex property name to analyze (e.g. ``"strength"``).

    Returns:
        ``(cutoff_index, sorted_uids)`` where ``cutoff_index`` is the number
        of elements above the cutoff and ``sorted_uids`` lists UIDs descending.
    """
    import numpy as np
    from scipy.ndimage import uniform_filter1d

    if key not in s.vp:
        raise ValueError(f"Vertex property '{key}' not found.")

    vals = [(s.vp.uid[v], float(s.vp[key][v])) for v in s.vertices() if s.vp[key][v] > 0]
    vals.sort(key=lambda x: -x[1])
    uids = [v[0] for v in vals]
    y = np.array([v[1] for v in vals])

    if len(y) < 3:
        return len(y), uids

    # Tangent slope = 1 on the [0, 1]-normalized ascending curve (rank vs node
    # strength). The slope is computed on a smoothed curve because with many
    # CREs the normalized x-step is tiny (~1/m), so raw np.gradient is dominated
    # by single-step jitter and crosses 1 at the very bottom of the curve.
    ys = y[::-1]                                   # ascending (y is descending)
    m = len(ys)
    xn = np.arange(m) / (m - 1)
    yn = (ys - ys[0]) / (ys[-1] - ys[0])
    w = max(11, m // 200)
    yn_s = uniform_filter1d(yn, size=w, mode="nearest")
    slope = np.gradient(yn_s, xn)
    crossed = np.where(slope >= 1.0)[0]
    i = int(crossed[0]) if len(crossed) else m
    cutoff = m - i                                 # count of hubs above contact
    if verbose:
        val = ys[i] if i < m else ys[-1]
        print(f"[INFO] Slope-1 on vp.{key} (value≈{val:.3g} at cutoff): "
              f"cutoff at {cutoff}/{m} ({100 * cutoff / m:.1f}%)")
    return cutoff, uids


def prime_hubs(s: Architecture, key: str = "n", gene: str = 'gene', *, verbose: bool = True):
    """Find prime genes: the genes of hub CREs (slope-1 knee on node strength).

    Pipeline:
        1. :func:`strength` — node strength (sum of incident ``ep[key]``), if
           not already present on the graph.
        2. :func:`elbow` — slope-1 knee on the strength curve → hub CREs.
        3. Each hub's gene (``vp[gene]``, from :func:`annotate`) is collected,
           split by whether the hub is a promoter or not.

    Requires ``vp.annot`` and ``vp[gene]`` (from :func:`annotate`).

    Returns:
        dict with ``'prime_genes'``, ``'promoter_genes'``, ``'enhancer_genes'``,
        ``'hub_uids'``, ``'cutoff'``, ``'promoter_uids'``, ``'enhancer_uids'``.
    """
    for req in ("annot", gene):
        if req not in s.vp:
            raise ValueError(f"vp.{req} missing — run annotate() first.")

    strength_name = f"strength_{key.replace('n_', '')}" if key != 'n' else 'strength'
    if strength_name not in s.vp:
        s.strength(key=key, name=strength_name, verbose=verbose)

    cutoff, sorted_uids = s.elbow(strength_name, verbose=verbose)
    hub_uids = sorted_uids[:cutoff]

    p_genes, e_genes = set(), set()
    prom_hubs, enh_hubs = [], []
    for uid in hub_uids:
        v = s.index[uid]
        g = s.vp[gene][v]
        if "Promoter" in s.vp.annot[v]:
            prom_hubs.append(uid)
            if g: p_genes.add(g)
        else:
            enh_hubs.append(uid)
            if g: e_genes.add(g)

    prime_genes = p_genes | e_genes

    if verbose:
        print(f"[INFO] Prime hubs: {len(hub_uids)} hubs "
              f"({len(prom_hubs)} promoters, {len(enh_hubs)} enhancers)")
        print(f"[INFO] Prime genes: {len(prime_genes)} = "
              f"{len(p_genes)} promoter + {len(e_genes)} enhancer "
              f"(overlap: {len(p_genes & e_genes)})")

    return {
        "prime_genes": prime_genes,
        "promoter_genes": p_genes,
        "enhancer_genes": e_genes,
        "hub_uids": hub_uids,
        "cutoff": cutoff,
        "promoter_uids": prom_hubs,
        "enhancer_uids": enh_hubs,
    }


# attach these utilities to the Architecture class to preserve the method API
Architecture.make = make
Architecture.add_mcool = add_mcool
Architecture.normalize = normalize
Architecture.prune = prune
Architecture.annotate = annotate
Architecture.strength = strength
Architecture.elbow = elbow
Architecture.prime_hubs = prime_hubs
