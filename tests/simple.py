import sys

sys.path.append(".")
sys.path.append("..")

import sage.all as sageall
from khc import find_planar_code, code_to_string


G = sageall.Graph()
for i in range(19):
    G.add_vertex(i)

G.add_edges(
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
