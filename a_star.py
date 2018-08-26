import heapq
import sys
from collections import Counter, deque
from operator import itemgetter


def bfs_tree(adjacency_list, start):
    """Sorts nodes in breadth-ﬁrst search order starting from the start node"""
    visited_nodes, queue = [start], deque([start])
    while queue:
        unvisited_node = queue.popleft()
        for neighbour in adjacency_list[unvisited_node]:
            if neighbour in visited_nodes:
                continue
            visited_nodes.append(neighbour)
            queue.append(neighbour)
    return visited_nodes


def list_difference(list_1, list_2):
    counter_1 = Counter(list_1)
    counter_2 = Counter(list_2)
    difference = counter_2 - counter_1 | counter_1 - counter_2
    return sum(difference.values())


def a_star(reference_molecule, molecule_to_map, start_nodes):
    def extend_mapping(mapping):
        """Expands the mapping by combining of the next atom from the breadth-ﬁrst order of the reference graph
        with all atoms of another graph which have the same chemical element.
        """
        try:
            next_node_ref = bfs_order[len(mapping)]
        except IndexError:
            return []
        subgraph_rest_map = molecule_to_map.subgraph(mapping, diff=True)
        candidate_nodes = filter(
            lambda node: reference_molecule.atomic_numbers[next_node_ref] == molecule_to_map.atomic_numbers[node],
            subgraph_rest_map.adjacency_list.keys())
        for candidate in candidate_nodes:
            extended_mapping = mapping[:] + [candidate]
            yield extended_mapping

    def accumulated_cost(nodes_map):
        """The accumulated cost g is calculated using the mapped part of the molecule,
        i.e. the subgraph induced by mapped nodes nodes_map."""
        nodes_ref = bfs_order[:len(nodes_map)]

        last_node_ref = nodes_ref[-1]
        last_node_map = nodes_map[-1]

        subgraph_ref = reference_molecule.subgraph(nodes_ref)
        subgraph_map = molecule_to_map.subgraph(nodes_map)

        neighbs_of_last_ref = subgraph_ref.adjacency_list[last_node_ref]
        neighbs_of_last_map = subgraph_map.adjacency_list[last_node_map]

        mapping = dict(zip(nodes_ref, nodes_map))

        common_neighbs = {mapping[node] for node in neighbs_of_last_ref} & neighbs_of_last_map

        subgraph_nodes_ref = neighbs_of_last_ref | {last_node_ref, }
        subgraph_nodes_map = common_neighbs | {last_node_map, }

        subgraph_ref = reference_molecule.subgraph(subgraph_nodes_ref)
        subgraph_map = molecule_to_map.subgraph(subgraph_nodes_map)

        bond_str_rest_ref = [''.join(sorted(itemgetter(i, j)(reference_molecule.atomic_numbers))) for i, j in
                             subgraph_ref.edges()]
        bond_str_rest_map = [''.join(sorted(itemgetter(i, j)(molecule_to_map.atomic_numbers))) for i, j in
                             subgraph_map.edges()]

        g = list_difference(bond_str_rest_ref, bond_str_rest_map)

        return g

    def future_cost(nodes_map):
        """The future bond cost h_bonds and the future atom cost h_atoms are estimated
        using the unmapped part of the molecule, i.e. the subgraph induced by unmapped nodes.
        nodes_map are mapped nodes.
        """
        subgraph_rest_ref = reference_molecule.subgraph(bfs_order[:len(nodes_map)], diff=True)
        subgraph_rest_map = molecule_to_map.subgraph(nodes_map, diff=True)

        atom_str_ref = [reference_molecule.atom_string(node) for node in subgraph_rest_ref.adjacency_list.keys()]
        atom_str_map = [molecule_to_map.atom_string(node) for node in subgraph_rest_map.adjacency_list.keys()]

        h_atoms = list_difference(atom_str_ref, atom_str_map) / 2

        # bonds cost
        bond_str_rest_ref = [''.join(sorted(itemgetter(i, j)(reference_molecule.atomic_numbers))) for i, j in
                             subgraph_rest_ref.edges()]
        bond_str_rest_map = [''.join(sorted(itemgetter(i, j)(molecule_to_map.atomic_numbers))) for i, j in
                             subgraph_rest_map.edges()]

        h_bonds = list_difference(bond_str_rest_ref, bond_str_rest_map)

        return max(h_atoms, h_bonds)

    # The A* Algorithm
    queue = []
    node_ref, node_map = start_nodes
    bfs_order = bfs_tree(reference_molecule.adjacency_list, start=node_ref)

    # The ﬁrst mapping pair with the path cost f = 0 is added to the priority queue.
    heapq.heappush(queue, (0, ([node_map], 0)))

    # The upper bound ub_f for the total path cost is initialised to inﬁnity.
    ub_f = float('inf')
    length_of_complete_mapping = len(bfs_order)
    possible_mappings = []
    
    # The future cost of complete mapping (both graphs are unmapped) is calculated
    # and the upper bound for the accumulated cost ub_g is set to this value.
    ub_g = future_cost([])

    # The algorithm stops if the priority queue is empty.
    while queue:

        # The mapping withdf the smallest total cost f (and corresponding is taken from the priority queue and considered.
        f, (mapping, g) = heapq.heappop(queue)

        # If total cost f of the current mapping is bigger than the upper bound for the total cost ub_f.
        # the algorithm stops, because only worse solutions (with higher costs) are left in the priority queue.
        if f > ub_f:
            break

        # The algorithm determines if the current mapping is complete
        elif len(mapping) == length_of_complete_mapping:
            ub_f = f # The upper bound for total cost ub_f is set to the total cost f of the current mapping.
            possible_mappings.append(mapping) # The complete mapping is added to the set of optimal mappings
            sys.stdout.write(f'\r{len(possible_mappings)} possible mappings:')
            sys.stdout.flush()

        # The incomplete mapping is expanded.
        else:
            for mapping_ext in extend_mapping(mapping):
                g_exp, h_exp = accumulated_cost(mapping_ext), future_cost(mapping_ext)
                g_exp += g # The accumulated cost of the newly mapped atoms g_exp is added to g
                f_exp = g_exp + h_exp

                # If the accumulated cost g_exp of the expanded mapping does not exceed the upper bound
                # for accumulated cost ub_g the mapping and its corresponding f and g values
                # are added to the priority queue. Otherwise, the mapping is pruned from the mapping process.
                if g_exp <= ub_g:
                    heapq.heappush(queue, (f_exp, (mapping_ext, g_exp)))

    mappings = [dict(zip(bfs_order, mapping)) for mapping in possible_mappings]

    return mappings
