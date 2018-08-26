"""A* atom mapper.

This script uses the A* atom mapping algorithm to generate all possible mappings between two molecules.

This script accepts two XYZ ﬁles with chemical elements and their cartesian coordinates as input.

The `plotly` package is required for visualisation of molecular graphs.
"""
from sys import argv

from a_star import a_star, list_difference
from molecular_graph import Graph
from plot import plot, gen_trace

# The algorithm takes two molecular graphs as input.
script, file1, file2 = argv

molecule1, molecule2 = Graph(), Graph()

molecule1.read_file(file1)
molecule2.read_file(file2)

trace1 = gen_trace(adjacency_list=molecule1.adjacency_list, elements=molecule1.elements,
                   x_coordinates=molecule1.x_coordinates, y_coordinates=molecule1.y_coordinates,
                   z_coordinates=molecule1.z_coordinates)

trace2 = gen_trace(adjacency_list=molecule2.adjacency_list, elements=molecule2.elements,
                   x_coordinates=molecule1.x_coordinates, y_coordinates=molecule2.y_coordinates,
                   z_coordinates=molecule2.z_coordinates)

element_labels1 = [dict(text=f'{el} ({num})', x=x, y=y, z=z, showarrow=False, yshift=15)
                   for num, (el, (x, y, z)) in enumerate(molecule1)]
element_labels2 = [dict(text=f'{el} ({num})', x=x, y=y, z=z, showarrow=False, yshift=15)
                   for num, (el, (x, y, z)) in enumerate(molecule2)]

# The graph with the largest number of edges is used as the reference graph.
if len(list(molecule2.edges())) > len(list(molecule1.edges())):
    plot(trace2, trace1, annotations1=element_labels2, annotations2=element_labels1)
else:
    plot(trace1, trace2, annotations1=element_labels1, annotations2=element_labels2)

text = """Please enter atom numbers of the first mapping pair, separated by a comma
(atom number from the left molecule, atom number from the right molecule): """

if len(molecule1) == len(molecule2) and not list_difference(molecule1.elements, molecule2.elements):
    while True:
        start_nodes = input(text).split(',')
        if len(start_nodes) == 2:
            # A* needs a starting point in the form of two mapped atoms in both graphs,
            # e.g. node_ref of molecule1 which can be mapped to the node_map of molecule2.
            node_ref = int(start_nodes[0])
            node_map = int(start_nodes[1])
            print()
            break
    mappings = a_star(molecule1, molecule2, [node_ref, node_map])
    print()
    for i, mapping in enumerate(mappings, start=1):
        print(f'{i}: {mapping}')
else:
    print('The A* algorithm can not map molecules with different numbers of atoms.')