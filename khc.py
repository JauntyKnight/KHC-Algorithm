import itertools
import networkx as nx

from collections import defaultdict

import sage.all as sageall
from sage.graphs import connectivity

import warnings

# warnings.filterwarnings("ignore")


def get_twin_edge(edge, tree, node):
    """
    Let T be SPQR tree, and let e be a virtual edge
    """

    u, v, _ = edge
    for neigh in tree.neighbor_out_iterator(node):
        neigh_skeleton = set(neigh[1].vertex_iterator())
        if u in neigh_skeleton and v in neigh_skeleton:
            # get the edge from u to v
            for u_, v_, label in get_skeleton(neigh).edge_iterator():
                if u_ == u and v_ == v:
                    return (u, v, label), neigh

    return None, None


def get_skeleton(node):
    return node[1].to_directed()


def find_code(edge, neigh, A, tree):
    type, _ = neigh

    if type == "R":
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
    print("Embedding:", embedding)

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
    print("Start edge:", edge)
    print("Direction", direction)

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
    print("entering S root")
    CV = {}
    skeleton = get_skeleton(root)
    for edge in filter(is_virtual, skeleton.edge_iterator()):
        u, v, label = edge
        twin_edge, neigh = get_twin_edge((u, v, label), tree, root)
        if twin_edge is None:
            continue

        CV[edge] = find_code(twin_edge, neigh, A, tree)

    def get_code_for_edge(edge):
        u, v, label = edge
        result = "(S "
        result += str(skeleton.size()) + " "

        e = next(filter(lambda x: x[1] != u and x[0] == v, skeleton.edges_incident(v)))
        e_in = e
        tour_counter = 1
        while True:
            if is_virtual(e):
                result += str(tour_counter) + " "
                result += CV[e]
            if e[1] in A:
                result += str(tour_counter) + " "
                result += "* "
                result += str(A[e[1]]) + " "
                del A[e[1]]

            e = next(
                filter(
                    lambda x: x[1] != e[0] and x[0] == e[1],
                    skeleton.edges_incident(e[1]),
                )
            )
            tour_counter += 1

            if e == e_in:
                break

        result += ")S "
        print("Result:", result)
        return result

    return min(
        (
            get_code_for_edge(edge)
            for edge in filter(is_virtual, skeleton.edge_iterator())
        )
    )


def code_of_P_root_node(root, A, tree):
    print("entering P root")
    CV = {}
    skeleton = get_skeleton(root)

    for edge in filter(is_virtual, skeleton.edge_iterator()):
        u, v, label = edge
        twin_edge, neigh = get_twin_edge((u, v, label), tree, root)
        if twin_edge is None:
            continue

        CV[edge] = find_code(twin_edge, neigh, A, tree)

    def get_code_for_vertex(v):
        result = "(P "
        result += str(skeleton.size()) + " "
        result += str(len(filter(is_virtual, skeleton.edge_iterator()))) + " "

        for code in sorted(CV.values()):
            result += code

        if v in A:
            result += "* "
            result += str(A[v]) + " "
            del A[v]

        result += ")P "
        return result

    return min(get_code_for_vertex(v) for v in skeleton.vertex_iterator())


def get_code_for_edge(skeleton, CV, A, edge, direction):
    code = "(R "

    print("Skeleton R")
    print(list(skeleton.vertex_iterator()))
    print(skeleton.is_directed())
    print(list(skeleton.edge_iterator()))
    skeleton_simple = skeleton.to_simple()

    tour = weinberg(skeleton, skeleton_simple, edge, direction)

    for edge in tour:
        print("Edge:", edge)
        if is_virtual(edge):
            code += CV[edge]
        if edge[0] in A:
            code += str(A[edge[1]]) + " "
            del A[edge[1]]

    if tour[-1][0] in A:
        code += str(A[tour[-1][0]]) + " "
        del A[tour[-1][0]]

    code += ")R "
    return code


def code_of_R_root_node(root, A, tree):
    print("entering R root")
    CV = {}
    skeleton = get_skeleton(root)

    for edge in filter(is_virtual, skeleton.edge_iterator()):
        u, v, label = edge
        twin_edge, neigh = get_twin_edge((u, v, label), tree, root)
        if twin_edge is None:
            continue

        CV[edge] = find_code(twin_edge, neigh, A, tree)

    return min(
        get_code_for_edge(skeleton, CV, A, edge, direction)
        for edge in filter(is_virtual, skeleton.edge_iterator())
        for direction in ["left", "right"]
    )


def find_code_R_non_root(e_in, node, A, tree):
    print("entering R non root")
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
    print("entering S non root")
    skeleton = get_skeleton(node)
    code = "(S "
    code += str(skeleton.size()) + " "
    u, v, label = e_in

    e = next(filter(lambda x: x[1] != u and x[0] == v, skeleton.edges_incident(v)))
    tour_counter = 1
    while e != e_in:
        if is_virtual(e):
            code += str(tour_counter) + " "
            twin_edge, neigh = get_twin_edge(e, tree, node)
            if twin_edge is None:
                continue
            code += find_code(twin_edge, neigh, A, tree)
        if e[1] in A:
            code += str(tour_counter) + " "
            code += "* "
            code += str(A[e[1]]) + " "
            del A[e[1]]

        e = next(
            filter(
                lambda x: x[1] != e[0] and x[0] == e[1],
                skeleton.edges_incident(e[1]),
            )
        )
        tour_counter += 1

    code += ")S "
    return code


def find_code_P_non_root(e_in, node, A, tree):
    print("entering P non root")
    skeleton = get_skeleton(node)

    CV = {}
    for edge in filter(
        lambda x: x[0] == e_in[0] and x != e_in,
        filter(is_virtual, skeleton.edge_iterator()),
    ):
        u, v, label = edge
        twin_edge, neigh = get_twin_edge((u, v, label), tree, node)
        if twin_edge is None:
            continue

        CV[edge] = find_code(twin_edge, neigh, A, tree)

    result = "(P "
    result += str(skeleton.size()) + " "
    result += (
        str(len([1 for edge in skeleton.edge_iterator() if is_virtual(edge)])) + " "
    )

    for code in sorted(CV.values()):
        result += code

    if e_in[1] in A:
        result += "* "
        result += str(A[e_in[1]]) + " "
        del A[e_in[1]]

    result += ")P "

    return result


def find_biconnected_codes_from_root(root, A, tree):
    print("Is directed:", tree.is_directed())
    print(list(tree.edge_iterator()))
    for node in tree.vertex_iterator():
        print("Node:", node, list(get_skeleton(node).vertex_iterator()))
        print(
            "Edges:",
            [list(get_skeleton(node).vertex_iterator()) for node in tree[node]],
        )

    print(
        "Centers:",
        [list(get_skeleton(center).vertex_iterator()) for center in tree.center()],
    )

    code = "(B "

    node_type = root[0]

    if node_type == "S" or node_type == "Q":
        code += code_of_S_root_node(root, A, tree)
    elif node_type == "P":
        code += code_of_P_root_node(root, A, tree)
    elif node_type == "R":
        code += code_of_R_root_node(root, A, tree)
    else:
        raise ValueError(f"Unknown node type: {node_type}")

    code += ")B "

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
    directed_edges = list(tree.breadth_first_search(center, edges=True))
    print("Directed edges:", list(directed_edges))
    tree_directed = sageall.DiGraph()
    tree_directed.add_vertices(tree.vertex_iterator())
    tree_directed.add_edges(tree.breadth_first_search(center, edges=True))

    # for edge in tree.breadth_first_search(center, edges=True):
    #     u, v = edge
    #     tree_directed.add_edge(u, v)

    return tree_directed


def find_biconnected_code(G, A, B):
    print("Subgraph:", G.nodes, G.edges)
    G_sage = sageall.Graph(G)
    tree = connectivity.TriconnectivitySPQR(G_sage).get_spqr_tree()

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
    bcs = list(nx.biconnected_components(G))
    articulation_points = set(nx.articulation_points(G))
    print("Biconnected components:", bcs)
    A = defaultdict(str)
    B = defaultdict(str)
    C = defaultdict(str)

    T = biconnected_components_to_tree(bcs)
    print("Tree:", T.nodes, T.edges)

    while len(T.nodes) > 1:
        leafs = [v for v in T.nodes if T.degree(v) == 1]

        for v in leafs[::-1]:
            print("Leaf:", v)
            C[v] = find_biconnected_code(G.subgraph(v), A, B[v])
            print("C[v]:", C[v])
            return

        # for all articulation points adjacent to leafs
        # leafs_neighbors = set(itertools.chain(*[T.neighbors(v) for v in leafs]))
        # for articulation_points in leafs_neighbors.intersection(articulation_points):
        #     # B[articulation_points] = find_bionnected_code(A, B[articulation_points])


G = nx.Graph()
for i in range(11):
    G.add_node(i)

G.add_edges_from(
    [
        (0, 2),
        (0, 8),
        (0, 5),
        (0, 3),
        (3, 1),
        (5, 2),
        (8, 9),
        (9, 2),
        (2, 1),
        (2, 4),
        (2, 6),
        (1, 7),
        (1, 6),
        (7, 4),
        (7, 6),
        (4, 6),
        (6, 10),
    ]
)

D = sageall.Graph(G).to_directed()
D.is_planar(set_embedding=True)
print(D._embedding)

print(find_planar_code(G))
