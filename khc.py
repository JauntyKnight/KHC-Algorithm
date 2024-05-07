import itertools
import networkx as nx

from collections import defaultdict

import sage.all as sageall
from sage.graphs import connectivity

SPECIAL_BASE = 10**5

R_CLOSE = SPECIAL_BASE + 10
R_OPEN = SPECIAL_BASE + 9
S_CLOSE = SPECIAL_BASE + 8
S_OPEN = SPECIAL_BASE + 7
P_CLOSE = SPECIAL_BASE + 6
P_OPEN = SPECIAL_BASE + 5
B_CLOSE = SPECIAL_BASE + 4
B_OPEN = SPECIAL_BASE + 3
A_CLOSE = SPECIAL_BASE + 2
A_OPEN = SPECIAL_BASE + 1
STAR = SPECIAL_BASE


def get_skeleton(node):
    """
    Returns the skeleton of a node in a SPQR tree,
    i.e. the subgraph represented by this node in the tree.

    Params:
        - `node`: A node of a SPQR tree
    Returns:
        The skeleton of the node
    """
    return node[1].to_directed()


def is_virtual(edge):
    """
    Checks if an edge is a virtual edge in a SPQR tree.

    Params:
        - `edge`: the edge to check
    Returns:
        `True` if the edge is a virtual edge, `False` otherwise
    """
    return edge[-1] is not None


def count_virtual_edges(skeleton):
    """
    Counts the number of virtual edges in a skeleton of a node in a SPQR tree.

    Params:
        - `skeleton`: the skeleton of a node in a SPQR tree
    Returns:
        The number of virtual edges in the skeleton
    """
    num_virtual_edges = 0

    for edge in skeleton.edge_iterator():
        if is_virtual(edge):
            num_virtual_edges += 1

    return num_virtual_edges


def next_edge_tour(edge, skeleton):
    """
    Returns the next edge in a tour of the skeleton of a node in a SPQR tree.

    Params:
        - `edge`: the edge to start the tour from
        - `skeleton`: the skeleton of the node
    Returns:
        The next edge in the tour
    """
    u, v, _ = edge
    return next(
        filter(
            lambda x: x[1] != u and x[0] == v,
            skeleton.edges_incident(v),
        )
    )


def get_virtual_edges_codes(edge_iterator, node, tree, A):
    """
    Computes the codes of the virtual edges in a node of a SPQR tree.

    Params:
        - `edge_iterator`: the virtual edges to compute the codes for
        - `node`: the node containing the edges
        - `tree`: the SPQR tree
        - `A`: the dictionary of articulation points and their codes
    Returns:
        A dictionary of the virtual edges and their codes
    """
    virtual_edges_codes = {}

    for edge in edge_iterator:
        twin_edge, neigh = get_twin_edge(edge, tree, node)
        if twin_edge is None:
            continue

        virtual_edges_codes[edge] = find_code(twin_edge, neigh, A, tree)

    return virtual_edges_codes


def articulation_to_iterable(articulation):
    """
    Converts an articulation point to an iterable.
    """

    return (item for component in articulation for item in component)


def get_twin_edge(edge, tree, node):
    """
    Finds the twin edge, together with the node containing it, of a virtual edge in a SPQR tree.

    Params:
        - `edge`: the edge to find the twin of
        - `tree`: the SPQR tree to search for the twin in
        - `node`: the node in the SPQR tree where the edge is contained

    Returns:
    Tuple of:
        - the twin edge
        - the node containing the twin edge
    Returns `(None, None)` in case `node` has no descendants in `T`.

    Let `T` be the SPQR tree, `v` be a vertex in `T`, and `e` be a virtual edge in `v`.
    The twin edge of `e` is the edge inside another vertex of `T`, namely `u`, such that
    `u` is a descendant of `v` in `T`, and both `e` and its twin connect the same two vertices
    of the original graph `G`.
    """

    u, v, _ = edge
    for neigh in tree.neighbor_out_iterator(node):
        neigh_skeleton = set(neigh[1].vertex_iterator())
        if u in neigh_skeleton and v in neigh_skeleton:
            # get the edge from u to v
            for edge_ in get_skeleton(neigh).edge_iterator():
                if edge == edge_:
                    return edge, neigh

    return None, None


def find_code(edge, neigh, A, tree):
    """
    Entry point for finding the code of a virtual edge in a SPQR tree.

    Params:
        - `edge`: the virtual edge to find the code of
        - `neigh`: the node containing the virtual edge
        - `A`: the dictionary of articulation points and their codes
        - `tree`: the SPQR tree
    """
    type, _ = neigh

    if type == "R":
        return find_code_R_non_root(edge, neigh, A, tree)
    elif type == "S":
        return find_code_S_non_root(edge, neigh, A, tree)
    elif type == "P":
        return find_code_P_non_root(edge, neigh, A, tree)
    else:
        raise ValueError(f"Unknown node type: {type}")


def find_code_S_helper(e_in, stop_edge, virtual_edges_codes, node, A):
    """
    Helper function to compute the code of an edge in a node of type S in a SPQR tree.

    Params:
        - `e_in`: starting edge of the tour
        - `stop_edge`: the edge to stop the tour at
        - `virtual_edges_codes`: the dictionary of virtual edges and their codes
        - `node`: the node containing the edge
        - `A`: the dictionary of articulation points and their codes
    """

    skeleton = get_skeleton(node)
    code = []
    code.append(S_OPEN)
    code.append(skeleton.size())

    e = e_in
    tour_counter = 1
    tour_edges_codes = []

    while True:
        if e in virtual_edges_codes:
            code.append(tour_counter)
            tour_edges_codes.extend(virtual_edges_codes[e])

        if e[1] in A:
            code.append(tour_counter)
            code.append(STAR)
            code.extend(articulation_to_iterable(A[e[1]]))

        e = next_edge_tour(e, skeleton)

        if e == stop_edge:
            break

        tour_counter += 1

    code.extend(tour_edges_codes)
    code.append(S_CLOSE)

    return code


def code_of_S_root_node(root, A, tree):
    """
    Computes the code of a root node of type S in a SPQR tree.

    Params:
        - `root`: the root node
        - `A`: the dictionary of articulation points and their codes
        - `tree`: the SPQR tree

    Returns:
        The code of the root node (as a list of integers).
    """
    skeleton = get_skeleton(root)

    virtual_edges_codes = get_virtual_edges_codes(
        filter(is_virtual, skeleton.edge_iterator()), root, tree, A
    )

    # if there are no virtua edges, and no articulation points
    if not virtual_edges_codes and not (set(skeleton.vertex_iterator()) & set(A)):
        code = []
        code.append(S_OPEN)
        code.append(skeleton.size())
        code.append(S_CLOSE)
        return code

    return min(
        (
            find_code_S_helper(
                next_edge_tour(edge, skeleton),
                next_edge_tour(edge, skeleton),
                virtual_edges_codes,
                root,
                A,
            )
        )
        for edge in skeleton.edge_iterator()
    )


def find_code_S_non_root(e_in, node, A, tree):
    """
    Computes the code of a non-root node of type S in a SPQR tree.

    Params:
        - `e_in`: the edge to compute the code for
        - `node`: the node containing the edge
        - `A`: the dictionary of articulation points and their codes
        - `tree`: the SPQR tree

    Returns:
        The code of the edge (as a list of integers).
    """

    skeleton = get_skeleton(node)

    virtual_edges_codes = get_virtual_edges_codes(
        filter(is_virtual, skeleton.edge_iterator()), node, tree, A
    )

    return find_code_S_helper(
        next_edge_tour(e_in, skeleton), e_in, virtual_edges_codes, node, A
    )


def code_of_P_root_node(root, A, tree):
    """
    Computes the code of a root node of type P in a SPQR tree.

    Params:
        - `root`: the root node
        - `A`: the dictionary of articulation points and their codes
        - `tree`: the SPQR tree

    Returns:
        The code of the root node (as a list of integers).
    """

    skeleton = get_skeleton(root)

    virtual_edges_codes = get_virtual_edges_codes(
        filter(is_virtual, skeleton.edge_iterator()), root, tree, A
    )

    def get_code_for_vertex(v):
        """
        Helper function to get the code of one of the two vertices of the P node.
        """
        code = []
        code.append(P_OPEN)
        code.append(skeleton.size())
        code.append(count_virtual_edges(skeleton))

        for edge, edge_code in sorted(
            filter(lambda x: x[0][0] == v, virtual_edges_codes.items()),
            key=lambda x: x[1],
        ):
            code.extend(edge_code)

        if v in A:
            code.append(STAR)
            code.extend([item for articulation in A[v] for item in articulation])
            # del A[v]

        other_vertex = next(filter(lambda x: x != v, skeleton.vertex_iterator()))

        if other_vertex in A:
            code.append(STAR)
            code.extend(
                [item for articulation in A[other_vertex] for item in articulation]
            )
            # del A[other_vertex]

        code.append(P_CLOSE)
        return code

    return min((get_code_for_vertex(v) for v in skeleton.vertex_iterator()))


def weinberg(G, G_simple, start_edge, direction):
    """
    Computes Weinberg's code of a triconnected planar graph.

    DOI: 10.1109/TCT.1966.1082573

    Params:
        - `G`: a triconnected planar graph
        - `G_simple`: the simple (mainly undirected) version of `G`
        - `start_edge`: the edge to start the tour from
        - `direction`: the direction to follow the tour in, either "left" or "right"

    Returns:
        The Weinberg code of the graph corresponding to `start_edge` and `direction`.
    """

    if not hasattr(G_simple, "_embedding") or G_simple._embedding is None:
        G_simple.is_planar(set_embedding=True)

    embedding = G_simple._embedding

    visited_vertices = set()
    visited_edges = set()

    def get_next_edge(edge, direction):
        """
        Returns the next unused edge in the direction
        """
        u, v = edge
        index = embedding[v].index(u)

        while True:
            if direction == "left":
                index = (index - 1) % len(embedding[v])
            elif direction == "right":
                index = (index + 1) % len(embedding[v])

            if (v, embedding[v][index]) not in visited_edges:
                return (v, embedding[v][index])

    edge = (start_edge[0], start_edge[1])

    code = []

    while True:
        visited_vertices.add(edge[0])
        visited_edges.add(edge)
        code.append((edge[0], edge[1], G.edge_label(*edge)[0]))

        if len(code) == G.size():
            return code

        next_vertex = edge[1]

        if next_vertex not in visited_vertices:
            # if arrived at a new vertex, exit via the first edge
            # from the embedding in the corresponding direction
            edge = get_next_edge(edge, direction)
        elif next_vertex in visited_vertices and edge[::-1] not in visited_edges:
            # if arrived at a visited vertex,
            # for which the reverse edge was not used,
            # exit via the reverse edge
            edge = edge[::-1]
        else:
            # if arrived at a visited vertex,
            # and the reverse edge was used,
            # exit via the next edge in the embedding
            edge = get_next_edge(edge, direction)


def get_code_for_edge_R_node(e_in, direction, skeleton, virtual_edges_codes, A):
    """
    Computes the code of an edge in a node of type R in a SPQR tree.

    Params:
        - `e_in`: the edge to compute the code for
        - `direction`: the direction to follow the tour in, either "left" or "right"
        - `skeleton`: the skeleton of the node
        - `virtual_edges_codes`: the dictionary of virtual edges and their codes
        - `A`: the dictionary of articulation points and their codes

    Returns:
        The code of the edge (as a list of integers).
    """
    code = []
    code.append(R_OPEN)

    # necessary because sage cannot compute embeddings for directed graphs
    skeleton_simple = skeleton.to_simple()

    # compute Weinberg's code
    tour = weinberg(skeleton, skeleton_simple, e_in, direction)
    visited_nodes = {}

    for edge in tour:
        if edge[0] not in visited_nodes:
            visited_nodes[edge[0]] = len(visited_nodes) + 1
        code.append(visited_nodes[edge[0]])
        if is_virtual(edge) and edge in virtual_edges_codes:
            code.extend(virtual_edges_codes[edge])
        if edge[0] in A:
            code.append(STAR)
            code.extend([item for articulation in A[edge[0]] for item in articulation])
            # del A[edge[1]]

    # account for the last vertex
    last_edge = tour[-1]
    if last_edge[1] not in visited_nodes:
        visited_nodes[last_edge[1]] = len(visited_nodes) + 1
    code.append(visited_nodes[last_edge[1]])

    if last_edge[1] in A:
        code.append(STAR)
        code.extend([item for articulation in A[last_edge[1]] for item in articulation])
        # del A[last_edge[1]]

    # if tour[-1][0] in A:
    #     code.extend([item for articulation in A[tour[-1][0]] for item in articulation])
    # del A[tour[-1][0]]

    code.append(R_CLOSE)
    return code


def code_of_R_root_node(root, A, tree):
    """
    Computes the code of a root node of type R in a SPQR tree.

    Params:
        - `root`: the root node
        - `A`: the dictionary of articulation points and their codes
        - `tree`: the SPQR tree

    Returns:
        The code of the root node (as a list of integers).
    """
    virtual_edges_codes = {}
    skeleton = get_skeleton(root)

    # compute the codes of the virtual edges
    for edge in filter(is_virtual, skeleton.edge_iterator()):
        twin_edge, neigh = get_twin_edge(edge, tree, root)
        virtual_edges_codes[edge] = find_code(twin_edge, neigh, A, tree)

    return min(
        get_code_for_edge_R_node(edge, direction, skeleton, virtual_edges_codes, A)
        for edge in skeleton.edge_iterator()
        for direction in ["left", "right"]
    )


def find_code_R_non_root(e_in, node, A, tree):
    """
    Computes the code of a non-root node of type R in a SPQR tree.

    Params:
        - `e_in`: the edge to compute the code for
        - `node`: the node containing the edge
        - `A`: the dictionary of articulation points and their codes
        - `tree`: the SPQR tree

    Returns:
        The code of the edge (as a list of integers).
    """

    virtual_edges_codes = {}
    skeleton = get_skeleton(node)

    # compute the codes of the virtual edges
    for edge in filter(is_virtual, skeleton.edge_iterator()):
        twin_edge, neigh = get_twin_edge(edge, tree, node)
        if twin_edge is None:
            continue

        virtual_edges_codes[edge] = find_code(twin_edge, neigh, A, tree)

    return min(
        get_code_for_edge_R_node(e_in, direction, skeleton, virtual_edges_codes, A)
        for direction in ["left", "right"]
    )


def find_code_P_non_root(e_in, node, A, tree):
    """
    Computes the code of a non-root node of type P in a SPQR tree.

    Params:
        - `e_in`: the edge to compute the code for
        - `node`: the node containing the edge
        - `A`: the dictionary of articulation points and their codes
        - `tree`: the SPQR tree

    Returns:
        The code of the edge (as a list of integers).
    """

    skeleton = get_skeleton(node)
    virtual_edges_codes = {}

    # compute the codes of the virtual edges in the same direction as e_in
    for edge in filter(
        lambda x: x[0] == e_in[0],
        filter(is_virtual, skeleton.edge_iterator()),
    ):
        twin_edge, neigh = get_twin_edge(edge, tree, node)
        if twin_edge is None:
            continue

        virtual_edges_codes[edge] = find_code(twin_edge, neigh, A, tree)

    code = []
    code.append(P_OPEN)
    code.append(skeleton.size())
    code.append(count_virtual_edges(skeleton))

    for edge_code in sorted(virtual_edges_codes.values()):
        code.extend(edge_code)

    if e_in[0] in A:
        code.append(STAR)
        code.extend([item for articulation in A[e_in[0]] for item in articulation])
        # del A[e_in[0]]

    if e_in[1] in A:
        code.append(STAR)
        code.extend([item for articulation in A[e_in[1]] for item in articulation])
        # del A[e_in[1]]

    code.append(P_CLOSE)

    return code


def find_biconnected_codes_from_root(root, A, tree):
    """
    Entry point for finding the code of a SPQR tree from root.

    Params:
        - `root`: the root node of the SPQR tree
        - `A`: the dictionary of articulation points and their codes
        - `tree`: the SPQR tree

    Returns:
        The code of the SPQR tree (as a list of integers).
    """

    code = []
    code.append(B_OPEN)

    node_type = root[0]

    if node_type == "S":
        code.extend(code_of_S_root_node(root, A, tree))
    elif node_type == "P" or node_type == "Q":
        code.extend(code_of_P_root_node(root, A, tree))
    elif node_type == "R":
        code.extend(code_of_R_root_node(root, A, tree))
    else:
        raise ValueError(f"Unknown node type: {node_type}")

    code.append(B_CLOSE)

    return code


def make_directed_tree(tree, center):
    """
    Returns a directed version of a tree with a given center.

    Params:
        - `tree`: the tree to make directed
        - `center`: the center of the tree
    Returns:
        The directed version of the tree
    """
    tree_directed = sageall.DiGraph()
    tree_directed.add_vertices(tree.vertex_iterator())
    tree_directed.add_edges(tree.breadth_first_search(center, edges=True))
    return tree_directed


def find_biconnected_code(G, A):
    """
    Computes the code of a biconnected component in a planar graph.

    Params:
        - `G`: the biconnected component
        - `A`: the dictionary of articulation points and their codes

    Returns:
        The code of the biconnected component (as a list of integers).
    """

    tree = connectivity.TriconnectivitySPQR(G).get_spqr_tree()

    min_code, _ = min(
        (
            (
                find_biconnected_codes_from_root(
                    center, A, make_directed_tree(tree, center)
                ),
                center,
            )
            for center in tree.center()
        ),
        key=lambda x: x[0],
    )

    return min_code


def find_planar_code(G):
    """
    Computes the code of a planar graph.

    See DOI:10.7155/jgaa.00094
    Algorithm and Experiments in Testing Planar Graphs for Isomorphism, by Kuklu et al
    """
    assert G.is_planar()

    if G.order() == 1:
        return [1]

    A = {}
    G = sageall.Graph(G)
    C = {}
    T = connectivity.blocks_and_cuts_tree(G)

    for v in T.vertex_iterator():
        if v[0] == "C":
            A[v[1]] = []

    # while there is more than one vertex in the tree
    while T.order() > 1:
        leaves = [v for v in T if T.degree(v) == 1]

        for leaf in leaves:
            C[leaf] = find_biconnected_code(G.subgraph(leaf[1]), A)

        # a leaf can only be adjacent to a cut node anyway
        # articulation_points = set(itertools.chain([T[leaf] for leaf in leaves]))
        articulation_points = set()
        for leaf in leaves:
            for neigh in T[leaf]:
                articulation_points.add(neigh)

        for articulation_point in articulation_points:
            biconnected_neighs = []
            for leaf in [v for v in T[articulation_point] if T.degree(v) == 1]:
                biconnected_neighs.append(C[leaf])

            previous_code = (
                A[articulation_point[1]][1:-1] if A[articulation_point[1]] else []
            )
            A[articulation_point[1]] = [[A_OPEN]]

            A[articulation_point[1]].extend(sorted(previous_code + biconnected_neighs))

            A[articulation_point[1]].append([A_CLOSE])

        # delete from T all leaves
        T.delete_vertices(leaves)

        # delete from T all articulation points with degree 1
        T.delete_vertices([v for v in articulation_points if T.degree(v) == 1])

    # the last node in the blocks and cuts tree
    v = next(T.vertex_iterator())

    if v[0] == "C":  # articulation point
        return [item for articulation in A[v[1]] for item in articulation]

    # else the last node is a biconnected component
    return find_biconnected_code(G.subgraph(v[1]), A)


def code_to_string(code):
    special_chars_map = {
        R_CLOSE: ")R",
        R_OPEN: "(R",
        S_CLOSE: ")S",
        S_OPEN: "(S",
        P_CLOSE: ")P",
        P_OPEN: "(P",
        B_CLOSE: ")B",
        B_OPEN: "(B",
        A_CLOSE: ")A",
        A_OPEN: "(A",
        STAR: "*",
    }

    return " ".join(
        special_chars_map[char] if char in special_chars_map else str(char)
        for char in code
    )
