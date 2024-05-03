import itertools
import networkx as nx

from collections import defaultdict

import sage.all as sageall
from sage.graphs import connectivity

SPECIAL_BASE = 9999999

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


def get_twin_edge(edge, tree, node):
    """
    Let T be SPQR tree, and let e be a virtual edge
    """

    u, v, label = edge
    for neigh in tree.neighbor_out_iterator(node):
        neigh_skeleton = set(neigh[1].vertex_iterator())
        if u in neigh_skeleton and v in neigh_skeleton:
            # get the edge from u to v
            for edge_ in get_skeleton(neigh).edge_iterator():
                if edge == edge_:
                    return edge, neigh

    return None, None


def get_skeleton(node):
    return node[1].to_directed()


def find_code(edge, neigh, A, tree):
    type, _ = neigh

    if type == "T":
        return find_code_R_non_root(edge, neigh, A, tree)
    elif type == "S":
        return find_code_S_non_root(edge, neigh, A, tree)
    elif type == "P":
        return find_code_P_non_root(edge, neigh, A, tree)
    else:
        raise ValueError(f"Unknown node type: {type}")


def is_virtual(edge):
    return edge[-1] is not None


def weinberg(G, G_simple, start_edge, direction):
    if not hasattr(G_simple, "_embedding") or G_simple._embedding is None:
        G_simple.is_planar(set_embedding=True)

    embedding = G_simple._embedding
    # print("Embedding:", embedding)

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
    # print("Start edge:", edge)
    # print("Direction", direction)

    while True:
        visited_vertices.add(edge[0])
        visited_edges.add(edge)
        code.append((edge[0], edge[1], G.edge_label(*edge)[0]))

        if len(code) == G.size():
            return code

        next_vertex = edge[1]

        if next_vertex not in visited_vertices:
            edge = get_next_edge(edge, direction)
        elif next_vertex in visited_vertices and edge[::-1] not in visited_edges:
            edge = edge[::-1]
        else:
            edge = get_next_edge(edge, direction)


def code_of_S_root_node(root, A, tree):
    # print("entering S root">
    CV = {}
    skeleton = get_skeleton(root)
    for edge in filter(is_virtual, skeleton.edge_iterator()):
        u, v, label = edge
        twin_edge, neigh = get_twin_edge((u, v, label), tree, root)
        if twin_edge is None:
            continue

        CV[edge] = find_code(twin_edge, neigh, A, tree)

    # if len(CV) == 0:
    #     return f"<S {skeleton.size()} )S "

    def get_code_for_edge(edge):
        u, v, label = edge
        result = [S_OPEN]
        result.append(skeleton.size())

        e = next(filter(lambda x: x[1] != u and x[0] == v, skeleton.edges_incident(v)))
        # print("Start edge:", edge)
        # print("First edge:", e)
        e_in = e
        tour_counter = 1

        virtual_edge_codes = []

        while True:
            if is_virtual(e):
                result.append(tour_counter)
                virtual_edge_codes.extend(CV[e])
            # print("E:", e)
            if e[1] in A:
                result.append(tour_counter)
                result.append(STAR)
                result.extend(A[e[1]])
                # del A[e[1]]

            # print("Counter:", tour_counter, e, result)
            e = next(
                filter(
                    lambda x: x[1] != e[0] and x[0] == e[1],
                    skeleton.edges_incident(e[1]),
                )
            )
            tour_counter += 1

            if e == e_in:
                break

        result.extend(virtual_edge_codes)
        result.append(S_CLOSE)
        # print("Result:", result)
        return result

    return min((get_code_for_edge(edge) for edge in skeleton.edge_iterator()))


def code_of_P_root_node(root, A, tree):
    # print("entering P root">
    CV = {}
    skeleton = get_skeleton(root)

    for edge in filter(is_virtual, skeleton.edge_iterator()):
        u, v, label = edge
        twin_edge, neigh = get_twin_edge((u, v, label), tree, root)
        if twin_edge is None:
            continue

        CV[edge] = find_code(twin_edge, neigh, A, tree)

    def get_code_for_vertex(v):
        result = [P_OPEN]
        result.append(skeleton.size())
        result.append(len([1 for edge in skeleton.edge_iterator() if is_virtual(edge)]))

        for code in sorted(CV.values()):
            result.extend(code)

        if v in A:
            result.extend(A[v])
            result.append(STAR)
            # del A[v]

        other_vertex = next(filter(lambda x: x != v, skeleton.vertex_iterator()))

        if other_vertex in A:
            result.extend(A[other_vertex])
            result.append(STAR)
            # del A[other_vertex]

        result.append(P_CLOSE)
        return result

    return min(get_code_for_vertex(v) for v in skeleton.vertex_iterator())


def get_code_for_edge(skeleton, CV, A, e_in, direction):
    code = [R_OPEN]
    # print("CV:", CV)

    # print("Skeleton T">
    # print(list(skeleton.vertex_iterator()))
    # print(skeleton.is_directed())
    # print(list(skeleton.edge_iterator()))
    skeleton_simple = skeleton.to_simple()

    tour = weinberg(skeleton, skeleton_simple, e_in, direction)
    visited_nodes = {}

    for edge in tour:
        # print("Edge:", edge)
        if edge[0] not in visited_nodes:
            visited_nodes[edge[0]] = len(visited_nodes) + 1
        code.append(visited_nodes[edge[0]])
        if e_in != edge and is_virtual(edge) and edge in CV:
            code.extend(CV[edge])
        if edge[0] in A:
            # print(edge[0], A)
            code.append(STAR)
            # del A[edge[1]]

    # account for the last vertex
    last_edge = tour[-1]
    if last_edge[1] not in visited_nodes:
        visited_nodes[last_edge[1]] = len(visited_nodes) + 1
    code.append(visited_nodes[last_edge[1]])

    if last_edge[1] in A:
        code.append(STAR)
        # code += str(A[last_edge[1]]) + " "
        # del A[last_edge[1]]

    if tour[-1][0] in A:
        code.extend(A[tour[-1][0]])
        # del A[tour[-1][0]]

    code.append(R_CLOSE)
    return code


def code_of_R_root_node(root, A, tree):
    CV = {}
    skeleton = get_skeleton(root)

    for edge in filter(is_virtual, skeleton.edge_iterator()):
        u, v, label = edge
        twin_edge, neigh = get_twin_edge((u, v, label), tree, root)
        if twin_edge is None:
            continue

        CV[edge] = find_code(twin_edge, neigh, A, tree)

    if len(CV) == 0:
        return min(
            get_code_for_edge(skeleton, CV, A, edge, direction)
            for edge in skeleton.edge_iterator()
            for direction in ["left", "right"]
        )

    return min(
        get_code_for_edge(skeleton, CV, A, edge, direction)
        for edge in filter(is_virtual, skeleton.edge_iterator())
        for direction in ["left", "right"]
    )


def find_code_R_non_root(e_in, node, A, tree):
    # print("entering T non root">
    CV = {}
    skeleton = get_skeleton(node)

    for edge in filter(is_virtual, skeleton.edge_iterator()):
        u, v, label = edge
        twin_edge, neigh = get_twin_edge((u, v, label), tree, node)
        if twin_edge is None:
            continue

        CV[edge] = find_code(twin_edge, neigh, A, tree)

    return min(
        get_code_for_edge(skeleton, CV, A, e_in, direction)
        for direction in ["left", "right"]
    )


def find_code_S_non_root(e_in, node, A, tree):
    # print("entering S non root">
    skeleton = get_skeleton(node)
    code = [S_OPEN]
    code.append(skeleton.size())
    u, v, label = e_in

    e = next(filter(lambda x: x[1] != u and x[0] == v, skeleton.edges_incident(v)))
    tour_counter = 1
    while e != e_in:
        if is_virtual(e):
            code.append(tour_counter)
            twin_edge, neigh = get_twin_edge(e, tree, node)
            if twin_edge is None:
                continue
            code.extend(find_code(twin_edge, neigh, A, tree))
        if e[1] in A:
            code.append(tour_counter)
            code.append(STAR)
            code.extend(A[e[1]])
            # del A[e[1]]

        e = next(
            filter(
                lambda x: x[1] != e[0] and x[0] == e[1],
                skeleton.edges_incident(e[1]),
            )
        )
        tour_counter += 1

    code.append(S_CLOSE)
    return code


def find_code_P_non_root(e_in, node, A, tree):
    # print("entering P non root">
    skeleton = get_skeleton(node)

    CV = {}
    # print("P edges:", list(skeleton.edge_iterator()))
    for edge in filter(
        lambda x: x[0] == e_in[0] and x != e_in,
        filter(is_virtual, skeleton.edge_iterator()),
    ):
        u, v, label = edge
        twin_edge, neigh = get_twin_edge((u, v, label), tree, node)
        if twin_edge is None:
            continue

        CV[edge] = find_code(twin_edge, neigh, A, tree)

    # print("P CV:", CV)

    result = [P_OPEN]
    result.append(skeleton.size())
    result.append(len([1 for edge in skeleton.edge_iterator() if is_virtual(edge)]))

    for code in sorted(CV.values()):
        result.extend(code)

    if e_in[1] in A:
        result.append(STAR)
        result.extend(A[e_in[1]])
        # del A[e_in[1]]

    result.append(P_CLOSE)

    return result


def find_biconnected_codes_from_root(root, A, tree):
    # print("Is directed:", tree.is_directed())
    # print("Skeleton:")
    # print(list(tree.edge_iterator()))
    # for node in tree.vertex_iterator():
    #     print(
    #         "Node:",
    #         node,
    #         list(get_skeleton(node).vertex_iterator()),
    #         list(get_skeleton(node).edge_iterator()),
    #     )
    #     print(
    #         "Edges:",
    #         [list(get_skeleton(node).vertex_iterator()) for node in tree[node]],
    #     )

    # print(
    #     "Centers:",
    #     [list(get_skeleton(center).vertex_iterator()) for center in tree.center()],
    # )

    code = [B_OPEN]

    node_type = root[0]

    if node_type == "S":
        code.extend(code_of_S_root_node(root, A, tree))
    elif node_type == "P" or node_type == "Q":
        code.extend(code_of_P_root_node(root, A, tree))
    elif node_type == "T":
        code.extend(code_of_R_root_node(root, A, tree))
    else:
        raise ValueError(f"Unknown node type: {node_type}")

    code.append(B_CLOSE)

    return code


def biconnected_components_to_tree(bcs):
    T = nx.Graph()

    for component in bcs:
        component_tuple = tuple(component)
        T.add_node(component_tuple)

        for v in T.nodes:
            if v == component_tuple:
                continue

            for u in v:
                if u in component:
                    T.add_edge(component_tuple, v)

    return T


def make_directed_tree(tree, center):
    ### List to store directed edges and avoid the gb
    # directed_edges = list(tree.breadth_first_search(center, edges=True))
    # print("Directed edges:", list(directed_edges))
    tree_directed = sageall.DiGraph()
    tree_directed.add_vertices(tree.vertex_iterator())
    tree_directed.add_edges(tree.breadth_first_search(center, edges=True))

    # for edge in tree.breadth_first_search(center, edges=True):
    #     u, v = edge
    #     tree_directed.add_edge(u, v)

    return tree_directed


def find_biconnected_code(G, A):
    tree = connectivity.TriconnectivitySPQR(G).get_spqr_tree()

    min_code, min_center = min(
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
    assert nx.is_planar(G)
    A = defaultdict(str)

    G = sageall.Graph(G)
    C = defaultdict(str)
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
            codes = []
            for leaf in [v for v in T[articulation_point] if T.degree(v) == 1]:
                codes.append(C[leaf])
                # print("Leaf:", leaf, C[leaf])

            A[articulation_point[1]] = [A_OPEN]
            for code in sorted(codes):
                A[articulation_point[1]].extend(code)

            A[articulation_point[1]].append(A_CLOSE)

        # print("A:", A)
        # print("Finished a round")

        # delete from T all leaves
        T.delete_vertices(leaves)

        # delete from T all articulation points with degree 1
        T.delete_vertices([v for v in articulation_points if T.degree(v) == 1])

    # the last node in the blocks and cuts tree
    v = next(T.vertex_iterator())

    if v[0] == "C":  # articulation point
        return A[v[1]]

    # else the last node is a biconnected component
    # print(v[1])
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


G = nx.Graph()
for i in range(19):
    G.add_node(i)

G.add_edges_from(
    [
        (0, 4),
        (1, 5),
        (2, 5),
        (3, 5),
        (2, 3),
        (5, 6),
        (4, 5),
        (5, 8),
        (8, 9),
        (7, 9),
        (4, 7),
        (9, 10),
        (10, 11),
        (11, 12),
        (12, 9),
        (9, 13),
        (13, 14),
        (14, 17),
        (17, 15),
        (15, 13),
        (15, 16),
        (17, 18),
    ]
)

print(code_to_string(find_planar_code(G)))
