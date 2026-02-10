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
                    
                    # Add edges (one for each active locus in anchor1)
                    for _ in a1_overlaps:
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
                    
                    # Add edges (one for each active locus in anchor2)
                    for _ in a2_overlaps:
                        v1 = G.index[uid1]
                        v2 = G.index[uid2]
                        if directed:
                            G.add_edge(v2, v1)  # v2 (active) → v1 (target)
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
    return C * (x**(-alpha))


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
    expected, fit_params = _pl_expect(distances, raw_weights)
    if verbose:
        print(f"[INFO] Power-law fit complete: alpha={fit_params['alpha']:.3f}, C={fit_params['C']:.3e}")
    norm_weights = raw_weights / np.maximum(expected, 1e-12)
    for i, e in enumerate(s.edges()):
        s.ep[name][e] = norm_weights[i]
    if verbose:
        print(f"[INFO] Set O/E weights for {s.n_links} intra-chromosomal edges to `ep.{name}`.")
    return s


def annotate(self, loci, genes, *, verbose=True):
    """Add gene and genomic region annotations to vertices.
    
    Args:
        loci: Loci collection (same one used to build the graph)
        genes: Genes object for annotation
    """
    # Add annotation mapping
    cre_a_df = genes.annotations(loci)
    annot_map = cre_a_df.set_index('uid')['annotation'].to_dict()
    
    # Add nearest gene mapping
    cre_g_df = genes.nearest_genes(loci)
    nearest_map = dict(zip(cre_g_df['Name'], cre_g_df['Name_b']))
    
    # Create vertex properties
    self.vp.gene = self.new_vp('string')
    self.vp.annot = self.new_vp('string')
    
    for v in self.vertices():
        uid = self.vp.uid[v]
        self.vp.gene[v] = nearest_map.get(uid, '')
        self.vp.annot[v] = annot_map.get(uid, '')
    
    if verbose:
        print(f"[INFO] Annotated {self.n_loci} loci with genes and regions.")
    
    return self


def draw(self, loci, region, *, 
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
         ax=None,
         **kwargs):
    """Draw the subgraph for CREs overlapping a genomic region.
    
    Args:
        loci: Loci collection (same one used to build the graph)
        region: Tuple of (chrom, start, end) or Locus object
        vertex_size_by: Vertex property name to scale node sizes (e.g., 'Agg_H1')
        edge_width_by: Edge property name to scale edge widths (default: 'w')
        vertex_size_range: (min_size, max_size) for vertices (default: (10, 50))
        edge_width_range: (min_width, max_width) for edges (default: (1, 10))
        vertex_color: Color for vertices - can be:
            - String (e.g., '#4A90E2') - single color for all
            - Vertex property name (e.g., 'Agg_H1') - color by values
            - None - defaults to blue
        edge_color: Color for edges (default: gray)
        figsize: Figure size as (width, height)
        layout: Layout algorithm ('spring', 'circular', 'kamada_kawai')
        show_labels: Whether to show node labels
        label_prop: Vertex property to use for labels (default: 'gene')
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
    
    # Create vertex filter for subgraph
    vfilt = self.new_vertex_property("bool")
    for v in self.vertices():
        vfilt[v] = self.vp.uid[v] in region_uids_in_graph
    
    # Create subgraph
    subgraph = gt.GraphView(self, vfilt=vfilt)
    
    if subgraph.num_vertices() == 0:
        print(f"[WARNING] No vertices in subgraph for region {chrom}:{start}-{end}")
        return ax
    
    # Prepare vertex sizes using prop_to_size
    if vertex_size_by and vertex_size_by in self.vp:
        vertex_sizes = prop_to_size(
            self.vp[vertex_size_by], 
            mi=vertex_size_range[0], 
            ma=vertex_size_range[1],
            power=1.0
        )
    else:
        vertex_sizes = np.mean(vertex_size_range)
    
    # Prepare edge widths using prop_to_size
    if edge_width_by and edge_width_by in self.ep:
        edge_widths = prop_to_size(
            self.ep[edge_width_by],
            mi=edge_width_range[0],
            ma=edge_width_range[1],
            power=1.0
        )
    else:
        edge_widths = np.mean(edge_width_range)
    
    # Prepare vertex colors
    if vertex_color is None:
        vertex_fill_color = '#4A90E2'
    elif isinstance(vertex_color, str):
        # Check if it's a property name or a color string
        if vertex_color in self.vp:
            # It's a property name - use it for coloring
            vertex_fill_color = self.vp[vertex_color]
        else:
            # It's a color string
            vertex_fill_color = vertex_color
    else:
        vertex_fill_color = vertex_color
    
    # Prepare labels
    if show_labels and label_prop in self.vp:
        labels = subgraph.new_vertex_property("string")
        for v in subgraph.vertices():
            labels[v] = str(self.vp[label_prop][v])
    else:
        labels = None
    
    # Choose layout
    if layout == 'spring':
        pos = gt.sfdp_layout(subgraph)
    elif layout == 'circular':
        pos = gt.circular_layout(subgraph)
    elif layout == 'kamada_kawai':
        pos = gt.kamada_kawai_layout(subgraph)
    else:
        pos = gt.sfdp_layout(subgraph)
    
    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    
    # Draw using graph-tool
    gt.graph_draw(
        subgraph,
        pos=pos,
        vertex_size=vertex_sizes,
        vertex_fill_color=vertex_fill_color,
        edge_pen_width=edge_widths,
        edge_color=edge_color,
        vertex_text=labels,
        vertex_font_size=10,
        output_size=tuple(int(x * 100) for x in figsize),
        mplfig=ax,
        **kwargs
    )
    
    # Set title
    ax.set_title(f"{self._name}: {chrom}:{start:,}-{end:,}\n"
                 f"{subgraph.num_vertices()} nodes, {subgraph.num_edges()} edges",
                 fontsize=14, pad=10)
    ax.axis('off')
    
    return ax


# attach these utilities to the Architecture class to preserve the old API
Architecture.make = make
Architecture.make_clique = make_clique
Architecture.make_spread = make_spread
Architecture.add_mcool = add_mcool
Architecture.normalize = normalize
Architecture.annotate = annotate
Architecture.draw = draw