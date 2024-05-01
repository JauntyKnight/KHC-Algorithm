import itertools
import networkx as nx

from collections import defaultdict

import sage.all as sageall
from sage.graphs import connectivity

import warnings

warnings.filterwarnings("ignore")


def code_of_S_root_node(skeleton, A, tree):
    code = "S"

    for child in tree.children(skeleton):
        code += A[child]

    return code


def find_biconnected_codes_from_root(root, A, tree):
    code = "(B"

    node_type, skeleton = root

    print("Root:", root)
    skeleton = skeleton.networkx_graph()
    print("Skeleton:", skeleton.nodes, skeleton.edges)

    if node_type == "S" or node_type == "Q":
        code += code_of_S_root_node(skeleton, A, tree)
    elif node_type == "P":
        code += code_of_P_root_node(skeleton, A, tree)
    elif node_type == "R":
        code += code_of_R_root_node(skeleton, A, tree)
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
