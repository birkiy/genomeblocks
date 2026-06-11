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


# attach methods defined elsewhere (keeps same top-level API in core shim)
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


@classmethod
def make_clique(cls, loci, *, name: str="Clique", verbose: bool=True):
    """Build a fully connected (clique) graph from a Loci collection.
    
    Each locus becomes a vertex, and all pairs are connected with an edge.
    Edge distances are automatically calculated and stored in ep.d.
    
    Args:
        loci: Loci collection to build the clique from
        name: Name for the Architecture graph
        verbose: Whether to print info messages
    
    Returns:
        Architecture instance with all loci fully connected
    """
    from itertools import combinations
    from tqdm import tqdm
    
    G = cls(name=name)
    
    # Add all vertices
    vertex_map = []
    for l in loci:
        v = G._add_vertex(l.uid)
        vertex_map.append(v)
    
    # Add edges for every pair (clique)
    n_pairs = len(vertex_map) * (len(vertex_map) - 1) // 2
    desc = '[INFO] Building clique edges' if verbose else None
    
    for i, (v1, v2) in enumerate(tqdm(combinations(vertex_map, 2), total=n_pairs, desc=desc, disable=not verbose)):
        e = G.add_edge(v1, v2)
        # Calculate distance between the two loci
        locus1 = loci[int(v1)]
        locus2 = loci[int(v2)]
        dist = locus1.distance_to(locus2)
        G.ep.d[e] = dist if dist is not NotImplemented else 0
        # Initialize other edge properties
        G.ep.w[e] = 1.0
        G.ep.n[e] = 0.0
    
    if verbose:
        print(f"[INFO] Created clique graph with {G.n_loci} nodes and {G.n_links} edges. ✓")
    
    return G


@classmethod
def make_spread(cls, source_loci, bedpe: str, *, 
                name: str="Spread", 
                r: int=2500, 
                dmax=1e9, 
                hops: int=1,
                directed: bool=False,
                verbose: bool=True):
    """Build an Architecture graph by spreading from source loci through loops.
    
    This method discovers and builds a network starting from source loci (e.g., promoters).
    Unlike `make()`, you don't need a pre-defined set of all CREs - the network is built
    entirely from what's in the BEDPE file.
    
    The algorithm works by spreading through loops:
    - Hop 1: Find all loop anchors that connect to source loci
    - Hop 2: Find all loop anchors that connect to Hop 1 anchors
    - Hop N: Continue spreading...
    
    Each loop anchor is represented by its CENTER point, which becomes a vertex in the graph.
    
    Args:
        source_loci: Loci collection of starting points (e.g., promoters)
        bedpe: Path to BEDPE file containing chromatin loops
        name: Name for the Architecture graph
        r: Radius around loop anchor centers to search for source loci (default: 2500bp)
        dmax: Maximum distance between loop anchors (default: 1e9, effectively no limit)
        hops: Number of hops to explore (1 = direct connections only, 2 = second-order, etc.)
        directed: If True, creates directed edges (source → target); if False, undirected
        verbose: Whether to print progress information
    
    Returns:
        Architecture instance with the spreading network
        
    Example:
        >>> promoters = genes.promoters(upstream=2000, downstream=500)
        >>> G = Architecture.make_spread(
        ...     source_loci=promoters,
        ...     bedpe="loops.bedpe",
        ...     hops=2,
        ...     directed=True
        ... )
        >>> # Now G contains: promoters + their loop partners + partners of partners
    
    Note:
        - Vertices are loop anchor centers (chr:pos format as UID)
        - You don't need a pre-defined Loci collection of all CREs
        - The network grows organically from the BEDPE file
        - Use hops=1 for direct connections, hops=2 for neighborhood exploration
    """
    from tqdm import tqdm
    from collections import defaultdict
    from .loci import Loci
    from .locus import Locus
    
    G = cls(name=name)
    G.set_directed(directed)
    
    if len(source_loci) == 0:
        if verbose:
            print(f"[WARNING] No source loci provided!")
        return G
    
    if verbose:
        print(f"[INFO] Starting with {len(source_loci)} source loci")
    
    # Build a temporary Loci collection from active anchors for fast overlap queries
    # Start with source loci
    active_loci = source_loci.copy()
    
    # Track all discovered anchors to avoid re-processing
    all_discovered_uids = {l.uid for l in source_loci}
    
    # Add source vertices to graph (using their centers as UIDs)
    for l in source_loci:
        #center = (loc.start + loc.end) // 2
        #uid = f"{loc.chrom}:{center}"
        G._add_vertex(l.uid)
    
    # Statistics
    total_loops = 0
    mapped_loops_per_hop = defaultdict(int)
    edge_set = set()

    # Perform spreading for each hop
    for hop in range(hops):
        if len(active_loci) == 0:
            if verbose:
                print(f"[INFO] No active loci at hop {hop+1}, stopping early")
            break

        new_loci = []

        # Scan BEDPE file
        with open(bedpe) as f:
            desc = f'[INFO] Hop {hop+1}/{hops}: Scanning loops'
            for line in tqdm(f, desc=desc, disable=not verbose):
                if line.startswith('#') or not line.strip():
                    continue

                if hop == 0:
                    total_loops += 1

                fields = line.strip().split()
                if len(fields) < 6:
                    continue

                try:
                    chrom1, start1, end1 = fields[0], int(fields[1]), int(fields[2])
                    chrom2, start2, end2 = fields[3], int(fields[4]), int(fields[5])
                except ValueError:
                    continue

                # Calculate centers
                mid1 = (start1 + end1) // 2
                mid2 = (start2 + end2) // 2

                # Check distance constraint
                if chrom1 == chrom2 and abs(mid1 - mid2) > dmax:
                    continue

                # Use cgr.overlap() to efficiently check if anchors overlap active loci
                # Check anchor1 against active loci
                a1_overlaps = [active_loci[j] for *_, j in active_loci.cgr.overlap(chrom1, mid1 - r, mid1 + r)]
                # Check anchor2 against active loci
                a2_overlaps = [active_loci[j] for *_, j in active_loci.cgr.overlap(chrom2, mid2 - r, mid2 + r)]

                # Create UIDs for anchors (centers)
                #uid1 = f"{chrom1}:{mid1}"
                uid1 = f"{chrom1}:{mid1-r}-{mid1+r}(.)"
                uid2 = f"{chrom2}:{mid2-r}-{mid2+r}(.)"

                # If anchor1 overlaps active loci, spread to anchor2
                if a1_overlaps and uid1 != uid2:
                    mapped_loops_per_hop[hop] += 1

                    # Add vertices
                    G._add_vertex(uid1)
                    G._add_vertex(uid2)

                    # Add edge (only once per uid pair)
                    edge_key = (uid1, uid2) if directed else tuple(sorted([uid1, uid2]))
                    if edge_key not in edge_set:
                        edge_set.add(edge_key)
                        v1 = G.index[uid1]
                        v2 = G.index[uid2]
                        G.add_edge(v1, v2)

                    # Mark anchor2 for next hop if not already discovered
                    if uid2 not in all_discovered_uids:
                        # Create a Locus object for this anchor center
                        new_loc = Locus(chrom2, mid2, mid2 + 1)  # Single bp at center
                        new_loci.append(new_loc)
                        all_discovered_uids.add(uid2)

                # If anchor2 overlaps active loci, spread to anchor1
                if a2_overlaps and uid1 != uid2:
                    if not a1_overlaps:  # Don't double-count
                        mapped_loops_per_hop[hop] += 1

                    # Add vertices
                    G._add_vertex(uid1)
                    G._add_vertex(uid2)

                    # Add edge (only once per uid pair)
                    if directed:
                        edge_key = (uid2, uid1)  # v2 (active) → v1 (target)
                    else:
                        edge_key = tuple(sorted([uid1, uid2]))
                    if edge_key not in edge_set:
                        edge_set.add(edge_key)
                        v1 = G.index[uid1]
                        v2 = G.index[uid2]
                        if directed:
                            G.add_edge(v2, v1)
                        else:
                            G.add_edge(v1, v2)

                    # Mark anchor1 for next hop if not already discovered
                    if uid1 not in all_discovered_uids:
                        new_loc = Locus(chrom1, mid1, mid1 + 1)
                        new_loci.append(new_loc)
                        all_discovered_uids.add(uid1)
        
        # Update for next hop
        if verbose:
            print(f"[INFO] Hop {hop+1}: {mapped_loops_per_hop[hop]} loops mapped, "
                  f"{len(new_loci)} new anchors discovered")
        
        # Next hop uses newly discovered loci
        if new_loci:
            active_loci = Loci(new_loci)
        else:
            active_loci = Loci([])
    
    if verbose:
        total_mapped = sum(mapped_loops_per_hop.values())
        pct_mapped = 100 * total_mapped / max(total_loops, 1) if total_loops > 0 else 0
        print(f"[INFO] Complete: {total_loops} loops scanned | {total_mapped} mapped ({pct_mapped:.1f}%) | "
              f"loci={G.n_loci}, links={G.n_links}")
        print(f"[INFO] Graph is {'directed' if directed else 'undirected'} ({'→' if directed else '↔'})")
    
    return G


def surrounds(s: Architecture, positions, cre_loci, window: int, mcool: str, *,
              resolution: Optional[int]=None, source: str="w",
              name: str="s", verbose: bool=True) -> Architecture:
    """Sum normalized contact frequency between each TSS-anchored CRE and all
    CREs within ±``window`` bp, stored as ``s.vp[name]``.

    Positions are mapped to ``cre_loci`` via overlap to pick anchor CREs. For
    each anchor, raw bin counts to neighbouring CREs are pulled in bulk from
    ``mcool`` (one CSR fetch per chromosome), divided by the power-law expected
    from ``s._pl_fits[source]`` (cached by :meth:`normalize`), and summed.
    Anchors not present in the graph are added; non-anchor vertices stay NaN.
    Run :meth:`normalize` first.
    """
    import cooler
    import numpy as np
    from tqdm import tqdm
    from collections import defaultdict

    if not hasattr(s, "_pl_fits") or source not in getattr(s, "_pl_fits", {}):
        raise RuntimeError(
            f"surrounds: no cached pl_fit for source='{source}'. "
            f"Call normalize(source='{source}') first."
        )
    alpha = float(s._pl_fits[source]["alpha"])
    C     = float(s._pl_fits[source]["C"])
    if not (np.isfinite(alpha) and np.isfinite(C)):
        raise RuntimeError(f"surrounds: pl_fit['{source}'] has non-finite params.")

    uri = f"{mcool}::resolutions/{resolution}" if resolution else mcool
    clr = cooler.Cooler(uri)
    res = clr.binsize
    if verbose:
        print(f"[INFO] surrounds: ±{window}bp on {res}bp cool, "
              f"pl_fit['{source}'] alpha={alpha:.3f} C={C:.3e}")

    uid_to_bin, chrom_of_bin, chrom_offsets = _build_uid_to_bin(cre_loci, clr)

    tss_to_anchor = positions.map(cre_loci)
    n_tss = len(tss_to_anchor)
    skipped_tss = sum(1 for vs in tss_to_anchor.values() if not vs)
    anchor_uids = {u for vs in tss_to_anchor.values() for u in vs}

    if name not in s.vp:
        s.vp[name] = s.new_vertex_property("float")
    nan_val = float("nan")
    for v in s.vertices():
        s.vp[name][v] = nan_val

    anchors_by_chrom = defaultdict(list)
    skipped_no_bin = 0
    for anchor_uid in anchor_uids:
        anchor = cre_loci[anchor_uid]
        bin_id = uid_to_bin.get(anchor_uid)
        if bin_id is None:
            skipped_no_bin += 1
            continue
        anchors_by_chrom[anchor.chrom].append((anchor_uid, anchor, bin_id))

    set_count = 0
    iterator = anchors_by_chrom.items()
    if verbose:
        iterator = tqdm(iterator, desc='[INFO] surrounds: per-chrom contacts',
                        total=len(anchors_by_chrom))
    for chrom, anchors in iterator:
        try:
            mat = clr.matrix(balance=False, sparse=True).fetch(chrom).tocsr()
        except Exception as e:
            if verbose: print(f"[WARN] surrounds: skipping chrom {chrom}: {e}")
            continue
        coff = chrom_offsets[chrom]
        n_bins = mat.shape[0]

        for anchor_uid, anchor, anchor_bin_global in anchors:
            anchor_bin = int(anchor_bin_global - coff)
            if anchor_bin < 0 or anchor_bin >= n_bins:
                continue
            ctr = anchor.center
            hits = cre_loci.cgr.overlap(chrom, max(0, ctr - window), ctr + window)
            neighbor_bins = []
            neighbor_dists = []
            for *_, j in hits:
                nb = cre_loci[j]
                if nb.uid == anchor_uid: continue
                gbin = uid_to_bin.get(nb.uid)
                if gbin is None: continue
                lbin = int(gbin - coff)
                if lbin < 0 or lbin >= n_bins: continue
                neighbor_bins.append(lbin)
                neighbor_dists.append(abs(anchor.center - nb.center))
            if not neighbor_bins:
                continue
            row = mat[anchor_bin].toarray().ravel()
            raw = row[np.asarray(neighbor_bins, dtype=np.int64)].astype(np.float64)
            # Floor at one bin width — power-law diverges as d → 0; mirrors add_patches.
            d = np.maximum(np.asarray(neighbor_dists, dtype=np.float64), float(res))
            expected = C * d ** (-alpha)
            norm = raw / np.maximum(expected, 1e-12)
            v = s._add_vertex(anchor_uid)
            s.vp[name][v] = float(norm.sum())
            set_count += 1

    if verbose:
        print(f"[INFO] surrounds: {n_tss} positions ({skipped_tss} no-overlap) | "
              f"{len(anchor_uids)} anchor CREs ({skipped_no_bin} without bin) | "
              f"set vp.{name} on {set_count} vertices.")
    return s


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
    for i, edge in tqdm(enumerate(edge_list),desc='[INFO] Calculating distances for edges', total=len(edge_list)):
        v1, v2 = edge.source(), edge.target()
        uid1, uid2 = s.vp.uid[v1], s.vp.uid[v2]
        dist = loci[uid1].distance_to(loci[uid2])
        s.ep.d[edge] = dist
        raw_weights[i] = s.ep[source][edge]
        distances[i] = dist
    # Zero-distance edges (overlapping/co-located loci) are a power-law
    # singularity: _pl_expect fits on positive distances only, and _pl_model
    # returns inf at distance 0 so their O/E falls to 0. They are not removed
    # here — call Architecture.prune() once at the end to drop them (removing
    # edges mid-pipeline desyncs the edge-index range across repeated calls).
    expected, fit_params = _pl_expect(distances, raw_weights)
    if verbose:
        print(f"[INFO] Power-law fit complete: alpha={fit_params['alpha']:.3f}, C={fit_params['C']:.3e}")
    # Cache the fit on the Architecture so `add_patches` (and any other
    # downstream step) can reuse it without refitting. Keyed by source
    # property name so multiple normalizations don't clobber each other.
    if not hasattr(s, "_pl_fits"):
        s._pl_fits = {}
    s._pl_fits[source] = dict(fit_params)
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


# ── Patch extraction (Akita / C.Origami style, HiChIP-correct) ─────────────
#
# What this is:
#   For each edge in the Architecture, extract a K×K matrix of contact
#   counts around the edge's (bin_a, bin_b) cell from the underlying
#   .mcool, then optionally normalize it by the SAME power-law fit
#   that `normalize` uses (NOT ICE — ICE's uniform-contact-probability
#   assumption is wrong for HiChIP, which is enriched at protein-bound
#   regions).
#
# Why this exists separately from add_mcool/normalize:
#   add_mcool stores ONE scalar count per edge (and distributes the
#   bin-pair count among co-binned edges). add_patches stores a K×K
#   neighborhood per edge — the raw bin-pair counts straight from the
#   cooler, which gives a richer, lower-variance target for downstream
#   models (Akita-style smoothing without distorting per-edge attribution).
#
# Speed:
#   - One sparse-CSR fetch per chromosome (cached) — 22-23 mcool calls
#     for hg38, not n_edges calls.
#   - Vectorized expected-grid computation per edge.
#   - Power-law expected reuses the fit from `normalize` if already done.

def _build_uid_to_bin(loci, clr):
    """Map each anchor uid to its nearest cooler bin id. Returns
    (uid_to_bin, chrom_of_bin, chrom_offsets) — same convention as
    add_mcool uses internally."""
    import numpy as np
    bins = clr.bins()[:][['chrom', 'start', 'end']].reset_index()
    from .loci import Loci
    from .locus import Locus
    bins_l = Loci(Locus(row[1], row[2], row[3]) for row in bins.itertuples(index=False))
    near = loci.nearest(bins_l)
    uid_to_bin = dict(zip(near['Name'], near['Name_b'].map(bins_l.uids)))
    chrom_of_bin = dict(zip(bins.index, bins['chrom']))
    # Per-chromosome bin offset (the bin id of the chromosome's first bin)
    chrom_offsets = bins.groupby('chrom', observed=True)['index'].min().astype(int).to_dict()
    return uid_to_bin, chrom_of_bin, chrom_offsets


def add_patches(s: Architecture, loci, mcool: str, *,
                resolution: "Optional[int]" = None,
                K: int = 3,
                source: str = "w",
                normalize: bool = True,
                clip: float = 0.0,
                smooth_sigma: float = 0.0,
                pseudocount: float = 1e-3,
                verbose: bool = True):
    """Extract K×K contact patches around each edge, normalized by the
    same power-law expectation as `Architecture.normalize`.

    Parameters
    ----------
    loci, mcool, resolution
        Same as add_mcool — the .mcool whose contacts you're sampling.
    K : int (odd, ≥1)
        Side length of the patch. K=3 → 3×3 (9 outputs/edge); K=5 → 5×5.
    source : str
        Edge property to use when fitting the genome-wide power-law
        decay (default 'w', i.e. raw counts already added by add_mcool).
        If `s._pl_fits[source]` is already cached (because `normalize`
        was called with the same source), the fit is reused.
    normalize : bool
        If True (default), divide each patch cell by its distance-
        expected value and take log2 → log2(O/E). If False, return raw
        counts.
    clip : float
        Clip log2(O/E) values to (-clip, +clip) (Akita uses 2.0). 0 disables.
    smooth_sigma : float
        Convolve each patch with a 5-wide 2D Gaussian of this sigma
        (Akita uses 1.0). 0 disables.
    pseudocount : float
        Added to numerator+denominator when computing O/E to avoid
        log(0). Same default as generate_log2n_dataset.py's --n-pseudocount.

    Returns
    -------
    np.ndarray of shape (n_edges, K, K), dtype float32
        Aligned with `s.edges()` iteration order. Edges that map to
        unknown bins (e.g. inter-chromosomal or off-end-of-chrom) get
        all-zero patches; the caller can detect these via the
        `valid` mask returned alongside.
    np.ndarray of shape (n_edges,), dtype bool
        Per-edge validity mask (True = patch is real, False = zero-filled).
    """
    import cooler
    import numpy as np
    from tqdm import tqdm
    from collections import defaultdict
    from scipy.signal import convolve2d

    if K < 1 or K % 2 == 0:
        raise ValueError(f"K must be a positive odd integer (got {K})")
    half = K // 2

    uri = f"{mcool}::resolutions/{resolution}" if resolution else mcool
    clr = cooler.Cooler(uri)
    res = clr.binsize
    if verbose:
        print(f"[INFO] add_patches: {K}x{K} patches at {res} bp resolution. 🪟")

    # ── 1. Bin mapping (same as add_mcool) ────────────────────────────────
    uid_to_bin, chrom_of_bin, chrom_offsets = _build_uid_to_bin(loci, clr)

    # ── 2. Per-edge metadata (chrom, bins, anchor distance) ───────────────
    edge_list = list(s.edges())
    n_edges = len(edge_list)
    edge_chrom = [None] * n_edges
    edge_bin_a = np.zeros(n_edges, dtype=np.int64)
    edge_bin_b = np.zeros(n_edges, dtype=np.int64)
    edge_dist  = np.zeros(n_edges, dtype=np.float64)   # bp distance, same as ep.d
    valid = np.zeros(n_edges, dtype=bool)
    for i, edge in enumerate(edge_list):
        v1, v2 = edge.source(), edge.target()
        uid1, uid2 = s.vp.uid[v1], s.vp.uid[v2]
        bin1 = uid_to_bin.get(uid1)
        bin2 = uid_to_bin.get(uid2)
        if bin1 is None or bin2 is None:
            continue
        ch1 = chrom_of_bin[bin1]
        ch2 = chrom_of_bin[bin2]
        if ch1 != ch2:
            continue
        edge_chrom[i] = ch1
        edge_bin_a[i] = bin1
        edge_bin_b[i] = bin2
        edge_dist[i]  = loci[uid1].distance_to(loci[uid2])
        valid[i] = True

    # ── 3. Power-law fit (reuse cached, else fit from edge raw counts) ────
    fit = None
    if hasattr(s, "_pl_fits") and source in getattr(s, "_pl_fits", {}):
        fit = s._pl_fits[source]
        if verbose:
            print(f"[INFO] add_patches: reusing cached pl_fit for '{source}': "
                  f"alpha={fit['alpha']:.3f}, C={fit['C']:.3e}")
    elif normalize:
        if source not in s.ep:
            raise ValueError(
                f"add_patches: source edge property '{source}' not found "
                f"and no cached pl_fit available. Call add_mcool + normalize "
                f"first, or pass normalize=False to get raw patches."
            )
        raw_w = np.array([s.ep[source][e] for e in edge_list], dtype=float)
        _, fit = _pl_expect(edge_dist, raw_w)
        if verbose:
            print(f"[INFO] add_patches: power-law fit on '{source}': "
                  f"alpha={fit['alpha']:.3f}, C={fit['C']:.3e}")
        if not hasattr(s, "_pl_fits"):
            s._pl_fits = {}
        s._pl_fits[source] = dict(fit)

    # ── 4. Group edges by chromosome ──────────────────────────────────────
    by_chrom = defaultdict(list)
    for i, ch in enumerate(edge_chrom):
        if ch is not None:
            by_chrom[ch].append(i)

    # ── 5. Pre-build offset grids (used in distance calc) ─────────────────
    offs = np.arange(-half, half + 1)
    da, db = np.meshgrid(offs, offs, indexing="ij")          # (K, K)
    rel_offset = (da - db).astype(np.float64) * res          # bp offset along diagonal

    # ── 6. Per-chromosome extraction ──────────────────────────────────────
    raw_patches = np.zeros((n_edges, K, K), dtype=np.float32)

    for chrom, idxs in tqdm(by_chrom.items(),
                            desc="[INFO] add_patches: per-chrom extraction",
                            total=len(by_chrom)):
        try:
            mat = clr.matrix(balance=False, sparse=True).fetch(chrom).tocsr()
        except Exception as e:
            print(f"[WARN] add_patches: skipping chrom {chrom}: {e}")
            for i in idxs:
                valid[i] = False
            continue
        coff = chrom_offsets[chrom]
        n_bins = mat.shape[0]

        for i in idxs:
            a = int(edge_bin_a[i] - coff)
            b = int(edge_bin_b[i] - coff)
            # Clamp + remember the slot offsets within the K×K patch
            a_lo = max(0, a - half); a_hi = min(n_bins, a + half + 1)
            b_lo = max(0, b - half); b_hi = min(n_bins, b + half + 1)
            if a_hi <= a_lo or b_hi <= b_lo:
                valid[i] = False
                continue
            sub = mat[a_lo:a_hi, b_lo:b_hi].toarray().astype(np.float32)
            # Slot offsets — handles patches that fall partially off the
            # chromosome edge.
            rs = a_lo - (a - half)
            cs = b_lo - (b - half)
            raw_patches[i, rs:rs+sub.shape[0], cs:cs+sub.shape[1]] = sub

    if verbose:
        n_valid = int(valid.sum())
        n_zero  = int((raw_patches.sum(axis=(1, 2)) == 0).sum())
        print(f"[INFO] add_patches: {n_valid}/{n_edges} edges have a patch; "
              f"{n_zero} of those have all-zero raw counts.")

    if not normalize:
        return raw_patches, valid

    # ── 7. Distance-correct (power-law) → log2(O/E) ───────────────────────
    if fit is None or not (np.isfinite(fit.get('alpha', np.nan))
                           and np.isfinite(fit.get('C', np.nan))):
        raise RuntimeError(
            "add_patches: cannot normalize — no valid power-law fit "
            "available. Either run normalize() first, or pass normalize=False."
        )
    alpha, C = float(fit['alpha']), float(fit['C'])

    # Per-edge expected grid: shape (n_edges, K, K)
    # d_cell_bp[i, da, db] = max(|edge_dist[i] + (da - db) * res|, res)
    # Floor at `res` because the power-law diverges as d → 0.
    if verbose:
        print(f"[INFO] add_patches: computing log2(O/E) "
              f"(alpha={alpha:.3f}, C={C:.3e}, clip={clip}, σ={smooth_sigma})")
    d_centered = edge_dist[:, None, None] + rel_offset[None, :, :]
    d_cell_bp = np.maximum(np.abs(d_centered), float(res))
    expected = C * d_cell_bp ** (-alpha)
    oe = (np.maximum(raw_patches, 0.0) + pseudocount) / np.maximum(expected, 1e-12)
    log2oe = np.log2(np.maximum(oe, 1e-12)).astype(np.float32)
    log2oe[~valid] = 0.0

    if clip > 0:
        log2oe = np.clip(log2oe, -clip, clip)

    if smooth_sigma > 0:
        # Per-edge 2D Gaussian smooth. For K ≤ 5 the convolution is cheap
        # enough that we just loop; for larger K replace with FFT-batched.
        r = 2  # 5-wide kernel as in Akita
        xs = np.arange(-r, r + 1)
        g = np.exp(-(xs[:, None] ** 2 + xs[None, :] ** 2) / (2 * smooth_sigma ** 2))
        kernel = (g / g.sum()).astype(np.float32)
        for i in range(n_edges):
            if valid[i]:
                log2oe[i] = convolve2d(log2oe[i], kernel,
                                       mode="same", boundary="symm").astype(np.float32)

    if verbose:
        v = log2oe[valid]
        if v.size:
            c = K // 2
            ctr = log2oe[valid, c, c]
            print(f"[INFO] add_patches: log2(O/E) range [{v.min():.2f}, {v.max():.2f}], "
                  f"mean {v.mean():.2f}; center cell mean {ctr.mean():.2f}.")

    return log2oe, valid


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


def _merge_nearby(self, sub_loci, merge_distance, edge_key='w', label_prop='gene',
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
        color_prop: Optional vertex PropertyMap on self whose values should be
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
        [l for l in sub_loci if l.uid in self.index],
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
        if label_prop in self.vp:
            labels = []
            for loc in cluster:
                if loc.uid in self.index:
                    lbl = str(self.vp[label_prop][self.index[loc.uid]])
                    if lbl and lbl not in labels:
                        labels.append(lbl)
            mg.vp.label[v] = "/".join(labels) if labels else mg.vp.uid[v]
        else:
            mg.vp.label[v] = mg.vp.uid[v]

        # Size value
        if vertex_size_by and vertex_size_by in self.vp:
            mg.vp.size_val[v] = sum(
                self.vp[vertex_size_by][self.index[loc.uid]]
                for loc in cluster if loc.uid in self.index
            )
        else:
            mg.vp.size_val[v] = float(len(cluster))

        # Color value
        if color_prop is not None:
            vals = [color_prop[self.index[loc.uid]]
                    for loc in cluster if loc.uid in self.index]
            if vals:
                if color_vtype == "string":
                    mg.vp.color_val[v] = next((x for x in vals if x), vals[0])
                elif "vector" in color_vtype:
                    arr = np.array([list(x) for x in vals], dtype=float)
                    mg.vp.color_val[v] = arr.mean(axis=0).tolist()
                else:
                    mg.vp.color_val[v] = float(np.mean(vals))

    # Aggregate inter-cluster edges
    ep_key = edge_key if edge_key and edge_key in self.ep else 'w'
    edge_weights = defaultdict(float)
    for ci, cluster in enumerate(clusters):
        for loc in cluster:
            if loc.uid not in self.index:
                continue
            v_orig = self.index[loc.uid]
            for nb in v_orig.all_neighbors():
                n_uid = self.vp.uid[nb]
                if n_uid not in uid_to_cluster:
                    continue
                cj = uid_to_cluster[n_uid]
                if ci == cj:
                    continue
                pair = (min(ci, cj), max(ci, cj))
                edge_weights[pair] += self.ep[ep_key][self.edge(v_orig, nb)]

    for (ci, cj), w in edge_weights.items():
        e = mg.add_edge(cluster_verts[ci], cluster_verts[cj])
        mg.ep.w[e] = w

    return mg, clusters, uid_to_cluster


def draw(self, loci, region, *,
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
    region_uids_in_graph = {uid for uid in region_uids if uid in self}

    if len(region_uids_in_graph) == 0:
        print(f"[WARNING] No graph vertices found in region {chrom}:{start}-{end}")
        return ax

    # ── Resolve vertex_color once, shared across paths ───────────────
    # color_prop: a PropertyMap on self (needs to be remapped under merge)
    # color_literal: a scalar color value (string / tuple) or None
    color_prop = None
    color_literal = None
    if vertex_color is None:
        color_literal = '#4A90E2'
    elif isinstance(vertex_color, str):
        if vertex_color in self.vp:
            color_prop = self.vp[vertex_color]
        else:
            color_literal = vertex_color
    elif isinstance(vertex_color, gt.PropertyMap):
        if vertex_color.get_graph().base is not self.base:
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
        mg, clusters, _ = self._merge_nearby(
            sub_loci, merge_distance,
            edge_key=edge_width_by, label_prop=label_prop,
            vertex_size_by=vertex_size_by,
            color_prop=color_prop,
        )
        subgraph = mg

        # Vertex sizes from aggregated size_val
        if vertex_size_by and vertex_size_by in self.vp:
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
        vfilt = self.new_vertex_property("bool")
        for v in self.vertices():
            vfilt[v] = self.vp.uid[v] in region_uids_in_graph

        subgraph = gt.GraphView(self, vfilt=vfilt)

        if subgraph.num_vertices() == 0:
            print(f"[WARNING] No vertices in subgraph for region {chrom}:{start}-{end}")
            return ax

        # Vertex sizes
        if vertex_size_by and vertex_size_by in self.vp:
            vertex_sizes = prop_to_size(
                self.vp[vertex_size_by],
                mi=vertex_size_range[0],
                ma=vertex_size_range[1],
                power=1.0,
            )
        else:
            vertex_sizes = np.mean(vertex_size_range)

        # Edge widths
        if edge_width_by and edge_width_by in self.ep:
            edge_widths = prop_to_size(
                self.ep[edge_width_by],
                mi=edge_width_range[0],
                ma=edge_width_range[1],
                power=1.0,
            )
        else:
            edge_widths = np.mean(edge_width_range)

        vertex_fill_color = color_prop if color_prop is not None else color_literal

        # Labels
        if show_labels and label_prop in self.vp:
            labels = subgraph.new_vertex_property("string")
            for v in subgraph.vertices():
                labels[v] = str(self.vp[label_prop][v])
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

    ax.set_title(f"{self._name}: {chrom}:{start:,}-{end:,}\n{title_extra}",
                 fontsize=14, pad=10)
    ax.axis('off')

    return ax


# ── Network metrics ──────────────────────────────────────────────────────────
def aggregate(s: Architecture, key: str = "n", name: str = "agg", *, verbose: bool = True) -> Architecture:
    """Compute node strength (sum of incident edge weights) per vertex.

    For each vertex, sums the edge property ``key`` across all incident edges.
    Stores the raw sum in ``vp[name]`` and the total-normalized version in
    ``vp['n' + name]``.

    Args:
        key:  Edge property to aggregate (default ``"n"`` = O/E weights).
        name: Vertex property name for the result (default ``"agg"``).

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

    nname = f"n{name}"
    s.vp[nname] = s.new_vp("float")
    total = s.vp[name].a.sum()
    if total > 0:
        s.vp[nname].a = s.vp[name].a / total

    if verbose:
        print(f"[INFO] Aggregated ep.{key} → vp.{name} (node strength) + vp.{nname} (normalized)")
    return s


def weighted_activity(s: Architecture, features, *, key: str = "n", name: str = "abc",
                      eps: float = 1e-12, verbose: bool = True) -> Architecture:
    """Per-vertex sum of ``ep[key] * GM_features(neighbour)`` over incident edges.

    For each vertex ``u``, sums over neighbours ``v``:
        ``ep[key][edge(u,v)] * (∏_f vp[f][v])^(1/n_features)``
    The vertex's own feature values are NOT in the geometric mean — only the
    neighbour's. Single feature collapses to ``Σ ep[key] * vp[f][v]`` over
    neighbours, the classical ABC numerator with activity sourced from the
    neighbour side.

    Asymmetric on each edge: edge ``(u, v)`` contributes ``w * gm[v]`` to ``u``
    and ``w * gm[u]`` to ``v``. ``vp[f]`` must be non-negative; zeros are
    clipped to ``eps`` to keep the log-space computation finite.
    """
    import numpy as np
    if isinstance(features, str):
        features = [features]
    if not features:
        raise ValueError("At least one vertex feature is required.")
    for f in features:
        if f not in s.vp:
            raise ValueError(f"Vertex property '{f}' not found.")
    if key not in s.ep:
        raise ValueError(f"Edge property '{key}' not found.")

    n_feat = len(features)
    log_vp = np.zeros(s.num_vertices(), dtype=float)
    for f in features:
        vp_arr = np.maximum(np.asarray(s.vp[f].a, dtype=float), eps)
        log_vp += np.log(vp_arr)
    gm = np.exp(log_vp / n_feat)

    edges = s.get_edges()
    src, tgt = edges[:, 0], edges[:, 1]
    w_arr = np.asarray(s.ep[key].a, dtype=float)

    out = np.zeros(s.num_vertices(), dtype=float)
    np.add.at(out, src, w_arr * gm[tgt])
    np.add.at(out, tgt, w_arr * gm[src])

    s.vp[name] = s.new_vertex_property("float")
    s.vp[name].a = out

    if verbose:
        nz = int((out > 0).sum())
        print(f"[INFO] weighted_activity(features={list(features)}, key='{key}') "
              f"→ vp.{name} ({nz}/{s.n_loci} non-zero, neighbour-only GM)")
    return s


def cluster(s: Architecture, key: str = "n", name: str = "lc", *, verbose: bool = True) -> Architecture:
    """Compute weighted local clustering coefficient per vertex.

    Uses graph-tool's ``local_clustering`` with optional edge weights.
    Stores the result in ``vp[name]``.

    Args:
        key:  Edge property for weighting (default ``"n"``).
              Use ``None`` for unweighted clustering.
        name: Vertex property name (default ``"lc"``).

    Returns:
        self (for chaining).
    """
    import graph_tool.clustering as gt_clust

    weight = s.ep[key] if key and key in s.ep else None
    s.vp[name] = gt_clust.local_clustering(s, weight=weight)

    if verbose:
        nz = sum(1 for v in s.vertices() if s.vp[name][v] > 0)
        print(f"[INFO] Local clustering → vp.{name}  ({nz}/{s.n_loci} non-zero)")
    return s


def elbow(s: Architecture, key: str, *, transform: str = "double_exp",
          method: str = "kneedle", verbose: bool = True):
    """Find a cutoff on a vertex property's sorted-descending curve.

    Sorts vertices by ``vp[key]`` descending, then selects a cutoff index by
    one of two methods:

    - ``"kneedle"`` (default): optionally applies a double-exponential
      transform (``exp(exp(y_norm))``) to amplify curvature, then returns the
      index of maximum perpendicular distance from the chord.
    - ``"slope1"``: on the ascending curve with both axes normalized to
      ``[0, 1]``, finds the point where the local tangent slope equals ``1``
      (the 45° knee). The slope is measured on a smoothed curve so it reflects
      the curve's shape rather than single-step jitter. The ``transform``
      argument is ignored for this method.

    Args:
        key:       Vertex property name to analyze (e.g. ``"agg"``).
        transform: ``"double_exp"`` (default) or ``"none"`` (kneedle only).
        method:    ``"kneedle"`` (default) or ``"slope1"``.

    Returns:
        ``(cutoff_index, sorted_uids)`` where ``cutoff_index`` is the number
        of elements above the cutoff and ``sorted_uids`` lists UIDs descending.
    """
    import numpy as np

    if key not in s.vp:
        raise ValueError(f"Vertex property '{key}' not found.")

    # Collect (uid, value) pairs, drop zeros
    vals = [(s.vp.uid[v], float(s.vp[key][v])) for v in s.vertices() if s.vp[key][v] > 0]
    vals.sort(key=lambda x: -x[1])
    uids = [v[0] for v in vals]
    y = np.array([v[1] for v in vals])

    if len(y) < 3:
        return len(y), uids

    # Normalize to [0, 1]
    y_norm = (y - y[-1]) / (y[0] - y[-1])

    if method == "slope1":
        # Tangent slope = 1 on the curve with both axes normalized to [0, 1].
        # Work on the ascending curve (rank vs node strength). The slope is
        # computed on a *smoothed* curve because with many CREs the normalized
        # x-step is tiny (~1/m), so raw np.gradient is dominated by single-step
        # jitter and crosses 1 at the very bottom of the curve. Smoothing makes
        # the slope reflect the actual shape, so the slope-1 point lands on the
        # real knee. Hubs are everything above the contact point.
        from scipy.ndimage import uniform_filter1d
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

    # Transform
    y_t = np.exp(np.exp(y_norm)) if transform == "double_exp" else y_norm

    # Kneedle: max perpendicular distance from chord
    x = np.arange(len(y_t), dtype=float)
    line_vec = np.array([x[-1] - x[0], y_t[-1] - y_t[0]])
    pts = np.column_stack([x - x[0], y_t - y_t[0]])
    distances = np.abs(np.cross(line_vec, pts)) / np.linalg.norm(line_vec)
    cutoff = int(np.argmax(distances))

    if verbose:
        print(f"[INFO] Elbow on vp.{key} ({transform}): cutoff at {cutoff}/{len(y)} "
              f"({100 * cutoff / len(y):.1f}%)")
    return cutoff, uids


def focus(s: Architecture, sources, key: str = "n", *, verbose: bool = True):
    """Find focus genes for source CREs via network-weighted candidate voting.

    For each source CRE we walk its promoter-annotated neighbors and tally edge
    weight ``ep[key]`` per *candidate* gene (from the neighbor's ``vp.genes``
    list — falls back to single ``vp.gene`` if ``annotate`` is pre-multi).
    This lets bidirectional promoters (TP53/WRAP53) and dense TSS clusters be
    disambiguated by contact strength rather than arbitrary pick.

    A per-source double-exponential Kneedle elbow on the sorted gene-weight
    vector selects the "focus" subset for that source.

    Args:
        sources: Iterable of UIDs (or single UID string).
        key:     Edge property used to weight the candidate vote (default ``"n"``).

    Returns:
        dict with keys:
            ``'genes'``   – union of focus genes across sources.
            ``'records'`` – per-source dict including ``focus_genes``,
                           ``gene_weights`` (full dict), ``n_candidates``,
                           ``n_focus``.
    """
    import numpy as np
    from collections import defaultdict

    if isinstance(sources, str):
        sources = [sources]
    if key not in s.ep:
        raise ValueError(f"Edge property '{key}' not found.")
    if "annot" not in s.vp:
        raise ValueError("vp.annot missing — run annotate() first.")
    use_multi = "genes" in s.vp
    has_single = "gene" in s.vp

    all_focus_genes = set()
    records = []

    for uid in sources:
        if uid not in s.index: continue
        v = s.index[uid]

        gene_weight: dict = defaultdict(float)
        for nb in v.all_neighbors():
            if "Promoter" not in s.vp.annot[nb]: continue
            e = s.edge(v, nb)
            w = float(s.ep[key][e])
            cands = []
            if use_multi and s.vp.genes[nb]:
                cands = [g for g in s.vp.genes[nb].split(";") if g]
            elif has_single and s.vp.gene[nb]:
                cands = [s.vp.gene[nb]]
            for gn in cands:
                gene_weight[gn] += w

        n = len(gene_weight)
        record = {"source": uid, "n_candidates": n, "gene_weights": dict(gene_weight)}

        if n == 0:
            record["focus_genes"] = set()
            record["n_focus"] = 0
            records.append(record)
            continue

        items = sorted(gene_weight.items(), key=lambda kv: -kv[1])
        w_sorted = np.array([kv[1] for kv in items], dtype=float)
        if n < 3:
            cutoff = n
        else:
            wmin, wmax = w_sorted[-1], w_sorted[0]
            if wmax == wmin:
                cutoff = n
            else:
                a = (w_sorted - wmin) / (wmax - wmin)
                y = np.exp(np.exp(a))
                x = np.arange(n, dtype=float)
                lv = np.array([x[-1] - x[0], y[-1] - y[0]])
                pts = np.column_stack([x - x[0], y - y[0]])
                dists = np.abs(np.cross(lv, pts)) / np.linalg.norm(lv)
                cutoff = int(np.argmax(dists)) + 1

        focus_genes = {items[i][0] for i in range(cutoff)} - {""}
        record["focus_genes"] = focus_genes
        record["n_focus"] = len(focus_genes)
        all_focus_genes |= focus_genes
        records.append(record)

    if verbose:
        print(f"[INFO] Focus (network-weighted): {len(sources)} sources → "
              f"{len(all_focus_genes)} focus genes.")
    return {"genes": all_focus_genes, "records": records}


def focus_genes(s: Architecture, sources, mapping=None, *, key: str="n", verbose: bool=True):
    """Focus genes via a gene→representative-CRE mapping (pseudo-supernodes).

    For each source uid, walks its neighbors and accumulates ``ep[key]`` into
    ``gene_weight[gene]`` whenever the neighbor is a representative CRE for
    that gene in ``mapping`` (or ``s.supernodes`` if ``mapping`` is None).
    A representative CRE can vote for multiple genes if its uid appears under
    several keys. Genes for which the source itself is a representative are
    excluded from that source's score (no self-vote). Per-source Kneedle
    elbow on the sorted weight vector picks the focus subset.

    Unlike :func:`focus`, this does not consult ``vp.annot``/``vp.genes`` —
    gene assignment comes from ``mapping``, typically built via
    :meth:`Genes.cre_supernodes`. Return shape mirrors :func:`focus`.
    """
    import numpy as np
    from collections import defaultdict

    if mapping is None:
        mapping = getattr(s, "supernodes", None)
        if mapping is None:
            raise RuntimeError(
                "focus_genes: no mapping given and s.supernodes not set. "
                "Pass mapping=... or assign s.supernodes = genes.cre_supernodes(...)."
            )
    if key not in s.ep:
        raise ValueError(f"Edge property '{key}' not found.")
    if isinstance(sources, str):
        sources = [sources]
    else:
        sources = list(sources)

    cre_to_genes: dict = defaultdict(list)
    for gene, cres in mapping.items():
        for cre in cres:
            cre_to_genes[cre].append(gene)

    all_focus_genes = set()
    records = []
    for uid in sources:
        if uid not in s.index: continue
        v = s.index[uid]
        self_genes = set(cre_to_genes.get(uid, []))

        gene_weight: dict = defaultdict(float)
        for nb in v.all_neighbors():
            nb_uid = s.vp.uid[nb]
            if nb_uid not in cre_to_genes: continue
            e = s.edge(v, nb)
            w = float(s.ep[key][e])
            for gn in cre_to_genes[nb_uid]:
                if gn in self_genes: continue
                gene_weight[gn] += w

        n = len(gene_weight)
        record = {"source": uid, "n_candidates": n, "gene_weights": dict(gene_weight)}
        if n == 0:
            record["focus_genes"] = set()
            record["n_focus"] = 0
            records.append(record)
            continue

        items = sorted(gene_weight.items(), key=lambda kv: -kv[1])
        w_sorted = np.array([kv[1] for kv in items], dtype=float)
        if n < 3:
            cutoff = n
        else:
            wmin, wmax = w_sorted[-1], w_sorted[0]
            if wmax == wmin:
                cutoff = n
            else:
                a = (w_sorted - wmin) / (wmax - wmin)
                y = np.exp(np.exp(a))
                x = np.arange(n, dtype=float)
                lv = np.array([x[-1] - x[0], y[-1] - y[0]])
                pts = np.column_stack([x - x[0], y - y[0]])
                dists = np.abs(np.cross(lv, pts)) / np.linalg.norm(lv)
                cutoff = int(np.argmax(dists)) + 1

        focus_set = {items[i][0] for i in range(cutoff)}
        record["focus_genes"] = focus_set
        record["n_focus"] = len(focus_set)
        all_focus_genes |= focus_set
        records.append(record)

    if verbose:
        print(f"[INFO] focus_genes (supernode mapping): {len(sources)} sources → "
              f"{len(all_focus_genes)} focus genes.")
    return {"genes": all_focus_genes, "records": records}


def supernode_metrics(s: Architecture, mapping=None, *, ep_keys=(), vp_keys=(),
                      vp_reduce: str = "sum", verbose: bool = True):
    """Per-supernode aggregations of edge and vertex properties.

    For each ``(gene, [cre_uid, ...])`` entry in ``mapping`` (or
    ``s.supernodes`` if ``mapping`` is None):
        - For each ``ep_keys`` entry: sums ``ep[k]`` across the union of edges
          incident to any CRE in the supernode (internal edges between two
          supernode members are counted once, not twice).
        - For each ``vp_keys`` entry: reduces ``vp[k]`` across supernode CREs
          via ``vp_reduce`` ∈ {"sum", "mean", "max"}.

    Returns a ``pandas.DataFrame`` indexed by gene name with columns
    ``[n_cres, *ep_keys, *vp_keys]``. CRE uids not present in the graph
    (e.g. window mapped to a CRE that no loop touched) are skipped and
    reported in the [INFO] line.
    """
    import pandas as pd
    if mapping is None:
        mapping = getattr(s, "supernodes", None)
        if mapping is None:
            raise RuntimeError(
                "supernode_metrics: no mapping given and s.supernodes not set."
            )
    if isinstance(ep_keys, str): ep_keys = (ep_keys,)
    if isinstance(vp_keys, str): vp_keys = (vp_keys,)
    ep_keys = tuple(ep_keys); vp_keys = tuple(vp_keys)
    for k in ep_keys:
        if k not in s.ep:
            raise ValueError(f"Edge property '{k}' not found.")
    for k in vp_keys:
        if k not in s.vp:
            raise ValueError(f"Vertex property '{k}' not found.")
    if vp_reduce not in ("sum", "mean", "max"):
        raise ValueError(f"vp_reduce must be one of sum/mean/max; got {vp_reduce!r}")

    rows = []
    n_missing_cres = 0
    for gene, cre_uids in mapping.items():
        valid = [u for u in cre_uids if u in s.index]
        n_missing_cres += len(cre_uids) - len(valid)
        row = {"gene": gene, "n_cres": len(valid), "degree": 0}

        if not valid:
            for k in ep_keys: row[k] = 0.0
            for k in vp_keys: row[k] = float("nan") if vp_reduce == "mean" else 0.0
            rows.append(row); continue

        # Edge dedup by sorted endpoint-idx pair, so an internal edge counts once.
        seen = set()
        unique_edges = []
        for u in valid:
            v = s.index[u]
            for e in v.all_edges():
                si, ti = int(e.source()), int(e.target())
                k = (si, ti) if si < ti else (ti, si)
                if k in seen: continue
                seen.add(k)
                unique_edges.append(e)
        row["degree"] = len(unique_edges)
        for k in ep_keys:
            row[k] = float(sum(s.ep[k][e] for e in unique_edges))

        for k in vp_keys:
            vals = [float(s.vp[k][s.index[u]]) for u in valid]
            if vp_reduce == "sum":   row[k] = sum(vals)
            elif vp_reduce == "mean": row[k] = sum(vals) / len(vals)
            else:                    row[k] = max(vals)
        rows.append(row)

    df = pd.DataFrame(rows).set_index("gene")
    if verbose:
        print(f"[INFO] supernode_metrics: {len(df)} genes | ep={list(ep_keys)} | "
              f"vp={list(vp_keys)} ({vp_reduce}) | {n_missing_cres} CRE uids not in graph")
    return df


def supernode_abc(s: Architecture, mapping=None, features=None, *, key: str = "n",
                  eps: float = 1e-12, verbose: bool = True):
    """ABC numerator per supernode, internal edges excluded.

    For each gene's supernode ``S``, sums over external edges ``(u, v)`` with
    ``u ∈ S, v ∉ S``:
        ``ep[key][edge] * (∏_f sqrt(vp[f][u] * vp[f][v]))^(1/|features|)``
    i.e. edge weight × geometric-mean-across-features of the endpoint-pair
    geometric mean. Both the supernode-side mark and the external CRE's mark
    appear. Edges with both endpoints inside ``S`` are skipped (no
    self-supernode contribution).

    Returns a ``pandas.DataFrame`` indexed by gene with columns
    ``[n_cres, abc]``.
    """
    import numpy as np
    import pandas as pd
    if mapping is None:
        mapping = getattr(s, "supernodes", None)
        if mapping is None:
            raise RuntimeError("supernode_abc: no mapping given and s.supernodes not set.")
    if features is None:
        raise ValueError("supernode_abc: 'features' is required.")
    if isinstance(features, str):
        features = [features]
    features = list(features)
    if not features:
        raise ValueError("supernode_abc: at least one feature required.")
    for f in features:
        if f not in s.vp:
            raise ValueError(f"Vertex property '{f}' not found.")
    if key not in s.ep:
        raise ValueError(f"Edge property '{key}' not found.")

    n_feat = len(features)
    log_vp = np.zeros(s.num_vertices(), dtype=float)
    for f in features:
        arr = np.maximum(np.asarray(s.vp[f].a, dtype=float), eps)
        log_vp += np.log(arr)
    # log_vp[v] = Σ_f log(vp[f][v]); for an edge (u, v) the combined endpoint-pair
    # GM-across-features = exp((log_vp[u] + log_vp[v]) / (2 * n_feat)).
    half_norm = 1.0 / (2.0 * n_feat)

    rows = []
    n_missing = 0
    for gene, cre_uids in mapping.items():
        valid = [u for u in cre_uids if u in s.index]
        n_missing += len(cre_uids) - len(valid)
        valid_set = {int(s.index[u]) for u in valid}
        score = 0.0
        for u in valid:
            v = s.index[u]
            u_idx = int(v)
            for e in v.all_edges():
                si, ti = int(e.source()), int(e.target())
                other = ti if si == u_idx else si
                if other in valid_set:
                    continue
                combined = np.exp((log_vp[u_idx] + log_vp[other]) * half_norm)
                score += float(s.ep[key][e]) * float(combined)
        rows.append({"gene": gene, "n_cres": len(valid), "abc": score})

    df = pd.DataFrame(rows).set_index("gene")
    if verbose:
        print(f"[INFO] supernode_abc: {len(df)} genes | features={features}, key='{key}' | "
              f"{n_missing} CRE uids not in graph")
    return df


def prime_hubs(s: Architecture, key: str = "n", gene: str = 'gene', *,
               method: str = "kneedle", verbose: bool = True):
    """Find prime genes: hub promoter genes ∪ focus genes of hub enhancers.

    Full pipeline:
        1. ``aggregate(key)`` — compute node strength if not yet present.
        2. ``elbow("agg", method=...)`` — find hub CREs above the cutoff.
        3. Split hubs into promoters (→ direct genes) and enhancers.
        4. ``focus(enhancer_hubs, key)`` — find focus neighbor genes.
        5. Union → prime genes.

    Requires ``vp.annot`` and ``vp.gene`` (from ``annotate()``).

    Args:
        key:    Edge property for weighting (default ``"n"``).
        method: Cutoff selection on the node-strength curve — ``"kneedle"``
                (Kneedle elbow, default) or ``"slope1"`` (tangent slope = 1,
                the 45° knee on the [0,1]-normalized ascending curve).

    Returns:
        dict with ``'prime_genes'``, ``'promoter_genes'``, ``'focus_genes'``,
        ``'hub_uids'``, ``'cutoff'``.
    """
    for req in ("annot", gene):
        if req not in s.vp:
            raise ValueError(f"vp.{req} missing — run annotate() first.")

    agg_name=f"agg_{key.replace('n_', '')}" if key != 'n' else 'agg'
    # 1. aggregate
    if agg_name not in s.vp:
        s.aggregate(key=key, name=agg_name, verbose=verbose)

    # 2. cutoff on node strength
    cutoff, sorted_uids = s.elbow(agg_name, method=method, verbose=verbose)
    hub_uids = sorted_uids[:cutoff]

    # 3. split by annotation — promoter hubs contribute *all* candidate genes
    p_genes = set()
    e_genes = set()
    prom_hubs, enh_hubs = [], []
    for uid in hub_uids:
        v = s.index[uid]
        
        annot_val = s.vp.annot[v]
        if "Promoter" in annot_val:
            prom_hubs.append(uid)
            p_genes.add(s.vp[gene][v])
        else:
            enh_hubs.append(uid)
            e_genes.add(s.vp[gene][v])

    # 5. union
    prime_genes = p_genes | e_genes

    if verbose:
        print(f"[INFO] Prime hubs: {len(hub_uids)} hubs "
              f"({len(prom_hubs)} promoters, {len(enh_hubs)} enhancers)")
        print(f"[INFO] Prime genes: {len(prime_genes)} = "
              f"{len(p_genes)} promoter + {len(e_genes)} focus "
              f"(overlap: {len(p_genes & e_genes)})")

    return {
        "prime_genes"   : prime_genes,
        "promoter_genes"   : p_genes,
        "enhancer_genes"   : e_genes,
        "hub_uids"      : hub_uids,
        "cutoff"        : cutoff,
        "promoter_uids" : prom_hubs,
        "enhancer_uids" : enh_hubs
    }


def mutual(s: Architecture, key: str = "n", *, transform: str = "double_exp", verbose: bool = True):
    """Find mutual-focus (reciprocal best hits) pairings between CREs.

    Applies the same double-exp + Kneedle machinery as :meth:`elbow` and
    :meth:`focus`, but *per row* of the contact matrix: for each vertex ``i``
    the sorted neighbor weights ``ep[key]`` define a focused tail F_i. The
    mutual mask is

        M = F & F^T

    i.e. edge (i, j) is mutual iff ``j ∈ F_i`` AND ``i ∈ F_j``. Equivalent to
    BLAST reciprocal best hits / mutual k-NN, but with a per-node data-driven
    cutoff rather than fixed k.

    Args:
        key:       Edge property used to rank neighbors (default ``"n"``).
        transform: ``"double_exp"`` (default) or ``"none"``.

    Returns:
        ``pandas.DataFrame`` with one row per mutual pair (``uid_a < uid_b``):
        ``uid_a, uid_b, weight, annot_a, annot_b, gene_a, gene_b`` (``genes_*``
        columns appear when ``vp.genes`` is present).
    """
    import numpy as np
    import pandas as pd
    from scipy import sparse

    if key not in s.ep:
        raise ValueError(f"Edge property '{key}' not found.")

    n = s.num_vertices()
    uids = np.array([s.vp.uid[v] for v in s.vertices()])

    # Build symmetric CSR of neighbor weights (each undirected edge → both directions)
    rows, cols, data = [], [], []
    for e in s.edges():
        w = float(s.ep[key][e])
        if w <= 0: continue
        i, j = int(e.source()), int(e.target())
        rows.append(i); cols.append(j); data.append(w)
        rows.append(j); cols.append(i); data.append(w)
    C = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))

    # Per-row Kneedle elbow on sorted neighbor weights → directed focus mask F
    F_rows, F_cols = [], []
    indptr, indices, values = C.indptr, C.indices, C.data
    for i in range(n):
        a, b = indptr[i], indptr[i + 1]
        if a == b: continue
        v_row = values[a:b]
        c_row = indices[a:b]
        order = np.argsort(-v_row)
        v_sorted = v_row[order]
        c_sorted = c_row[order]
        m = len(v_sorted)
        if m < 3:
            cutoff = m
        else:
            vmin, vmax = v_sorted[-1], v_sorted[0]
            if vmax == vmin:
                cutoff = m
            else:
                yn = (v_sorted - vmin) / (vmax - vmin)
                y = np.exp(np.exp(yn)) if transform == "double_exp" else yn
                x = np.arange(m, dtype=float)
                lv = np.array([x[-1] - x[0], y[-1] - y[0]])
                pts = np.column_stack([x - x[0], y - y[0]])
                dists = np.abs(np.cross(lv, pts)) / np.linalg.norm(lv)
                cutoff = int(np.argmax(dists)) + 1
        F_rows.extend([i] * cutoff)
        F_cols.extend(c_sorted[:cutoff].tolist())

    if not F_rows:
        if verbose:
            print(f"[INFO] Mutual focus on ep.{key}: 0 reciprocal pairs (empty focus).")
        return pd.DataFrame(columns=["uid_a", "uid_b", "weight",
                                     "annot_a", "annot_b", "gene_a", "gene_b"])

    F = sparse.csr_matrix(
        (np.ones(len(F_rows), dtype=np.int8), (F_rows, F_cols)),
        shape=(n, n),
    )
    M = sparse.triu(F.multiply(F.T), k=1).tocoo()

    has_annot = "annot" in s.vp
    has_multi = "genes" in s.vp
    has_gene  = "gene" in s.vp

    # Look up edge weights via the original CSR (already has both directions)
    Clil = C.tolil()
    records = []
    for i, j in zip(M.row.tolist(), M.col.tolist()):
        vi, vj = s.vertex(i), s.vertex(j)
        rec = {
            "uid_a":  uids[i],
            "uid_b":  uids[j],
            "weight": float(Clil[i, j]),
        }
        if has_annot:
            rec["annot_a"] = s.vp.annot[vi]
            rec["annot_b"] = s.vp.annot[vj]
        if has_gene:
            rec["gene_a"] = s.vp.gene[vi]
            rec["gene_b"] = s.vp.gene[vj]
        if has_multi:
            rec["genes_a"] = s.vp.genes[vi]
            rec["genes_b"] = s.vp.genes[vj]
        records.append(rec)

    df = pd.DataFrame(records).sort_values("weight", ascending=False).reset_index(drop=True)

    if verbose:
        print(f"[INFO] Mutual focus on ep.{key} ({transform}): "
              f"{F.nnz} directed focus edges → {len(df)} reciprocal pairs.")
    return df



# attach these utilities to the Architecture class to preserve the old API
Architecture.make = make
Architecture.make_clique = make_clique
Architecture.make_spread = make_spread
Architecture.surrounds = surrounds
Architecture.add_mcool = add_mcool
Architecture.normalize = normalize
Architecture.prune = prune
Architecture.add_patches = add_patches
Architecture.annotate = annotate
Architecture._merge_nearby = _merge_nearby
Architecture.draw = draw
Architecture.aggregate = aggregate
Architecture.weighted_activity = weighted_activity
Architecture.cluster = cluster
Architecture.elbow = elbow
Architecture.focus = focus
Architecture.focus_genes = focus_genes
Architecture.supernode_metrics = supernode_metrics
Architecture.supernode_abc = supernode_abc
Architecture.prime_hubs = prime_hubs
Architecture.mutual = mutual