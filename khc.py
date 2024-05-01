import itertools
import networkx as nx

from collections import defaultdict

import sage.all as sageall
from sage.graphs import connectivity

import warnings

warnings.filterwarnings("ignore")


def get_twin_edge(edge, tree, node):
    """
    Let T be SPQR tree, and let e be a virtual edge
    """

    u, v, _ = edge
    for neigh in tree[node]:
        neigh_skeleton = set(neigh[1].vertex_iterator())
        if u in neigh_skeleton and v in neigh_skeleton:
            # get the edge from u to v
            for u_, v_, label in get_skeleton(neigh).edge_iterator():
                if u_ == u and v_ == v:
                    return (u, v, label), neigh

    raise Exception("Edge not found")


def get_skeleton(node):
    return node[1].to_directed()


def find_code(edge, neigh, A, tree):
    return ""


def is_virtual(edge):
    return edge[-1] is not None


def code_of_S_root_node(root, A, tree):
    CV = {}
    skeleton = get_skeleton(root)
    for edge in filter(is_virtual, skeleton.edge_iterator()):
        u, v, label = edge
        twin_edge, neigh = get_twin_edge((u, v, label), tree, root)

        CV[edge] = find_code(twin_edge, neigh, A, tree)

    def get_code_for_edge(edge):
        u, v, label = edge
        result = "(S"
        result += str(skeleton.size())

        e = next(filter(lambda x: x[1] != u and x[0] == v, skeleton.edges_incident(v)))
        e_in = e
        tour_counter = 1
        while True:
            if is_virtual(e):
                result += str(tour_counter)
                result += CV[e]
            if e[1] in A:
                result += str(tour_counter)
                result += "*"
                result += str(A[e[1]])
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

        result += ")S"
        print("Result:", result)
        return result

    return min(
        (
            get_code_for_edge(edge)
            for edge in filter(is_virtual, skeleton.edge_iterator())
        )
    )


def find_biconnected_codes_from_root(root, A, tree):
    code = "(B"

    node_type = root[0]

    print("Root:", root)

    if node_type == "S" or node_type == "Q":
        code += code_of_S_root_node(root, A, tree)
        print("Code:", code)
    elif node_type == "P":
        code += code_of_P_root_node(root, A, tree)
    elif node_type == "R":
        code += code_of_R_root_node(root, A, tree)
    else:
        raise ValueError(f"Unknown node type: {node_type}")

    code += ")B"

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


def find_biconnected_code(G, A, B):
    print("Subgraph:", G.nodes, G.edges)
    G_sage = sageall.Graph(G)
    tree = connectivity.TriconnectivitySPQR(G_sage).get_spqr_tree()

    min_code, min_center = min(
        (
            (find_biconnected_codes_from_root(center, A, tree), center)
            for center in tree.center()
        ),
        key=lambda x: x[0],
    )

    return tree, min_code, min_center


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

find_planar_code(G)
