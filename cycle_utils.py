import time

from node_utils import (
    OntologyBranch,
    compute_children_overlap_ratio,
    compute_parent_overlap_ratio,
    disconnect_node,
    reconnect_node
)

class Cycle(object):
    __slots__ = ["elements", "parents"]
    def __init__(self, elements, parents=None):
        self.elements = elements
        self.parents = parents

    def __len__(self):
        return len(self.elements)

    def __str__(self):
        return "Cycle(%r)" % (self.elements)

    def __repr__(self):
        return str(self)


def find_ncycles(node, n, activated, origin=None):
    if origin is None:
        origin = node
    found = []
    for parent in node.parents:
        if parent.name not in activated:
            if parent == origin:
                found.append([node, parent])
            else:
                if n > 1:
                    found.extend([[node] + els for els in find_ncycles(parent, n - 1, activated, origin)])
    return found


def break_up_cycle(nodes, cycle, lookup_table, newname=None):
    print("num_cycles %d" % (len(cycle.parents),))

    joinable = []
    for parent_cycle in cycle.parents:
        if sum(len(n.children) for n in parent_cycle.elements) < 1000:
            joinable.append(parent_cycle)

    if len(joinable) > 0:
        print("len(joinable) = %d" % (len(joinable),))

        concerned = set([el for parent_cycle in joinable for el in parent_cycle.elements])

        merged = merge_many_nodes(
            concerned, lookup_table
        )

        nodes = nodes - concerned
        nodes = nodes | {merged}


    remainder = nodes

    # find most incestuous node:
    while True:
        found = disconnect_most_incestuous(
            remainder,
            lambda node: len(node.children) * (1.0 if any(n in remainder for n in node.children) else 0.0)
        )

        cyclical_node = cycles_exist(remainder, 10)

        if cyclical_node is None:
            break

        if not found:
            print("Could not perform cycle removal")
            break

    for node in remainder:
        assert(node.name in lookup_table), "Could not find %r in lookup_table" % (node.name,)

    return remainder


def disconnect_most_incestuous(family, order_fun):
    most_incestuous = sorted(family,
                             key=order_fun,
                             reverse=True)
    # lambda node: sum(parent in family for parent in node.parents)
    found = False
    for node in most_incestuous:
        newchildren = [child for child in node.children if child not in family]

        if newchildren == 0:
            continue
        else:
            found = True
            disconnect_node(node)

            node.children = newchildren

            reconnect_node(node)
            break
    return found


def merge_many_nodes(nodes, lookup_table, newname=None):
    if newname is None:
        newname = "join(" + ",".join([node.name for node in nodes]) + ")"
    newnode = OntologyBranch(newname)
    newnode.history = nodes

    nodes_set = set(nodes)

    all_children = list(set([child for node in nodes for child in node.children]) - nodes_set)
    all_parents = list(set([parent for node in nodes for parent in node.parents]) - nodes_set)
    newnode.children = all_children
    newnode.parents = all_parents

    for node in nodes:
        del lookup_table[node.name]

    lookup_table[newnode.name] = newnode

    for parent in newnode.parents:
        parent.children = [child for child in parent.children if child not in nodes_set] + [newnode]

    for child in newnode.children:
        child.parents = [parent for parent in child.parents if parent not in nodes_set] + [newnode]

    if any(node.lookup_table is not None for node in nodes):
        newnode.lookup_table = lookup_table

    return newnode

def merge_cycles(cycles):
    node2cycle = {}
    for cycle in cycles:
        for node in cycle.elements:
            if node in node2cycle:
                node2cycle[node].append(cycle)
            else:
                node2cycle[node] = [cycle]

    remaining_keys = list(node2cycle.keys())
    new_cycles = []

    while len(remaining_keys) > 0:
        key = remaining_keys[0]
        joined_set = set([key])
        prev_size = len(joined_set)

        involved_cycles = []

        while True:
            current_join_set = list(joined_set)
            for subkey in current_join_set:
                if subkey in node2cycle:
                    for cycle in node2cycle[subkey]:
                        joined_set = joined_set | cycle.elements
                        involved_cycles.append(cycle)
                    del node2cycle[subkey]

            new_size = len(joined_set)

            if prev_size == new_size:
                break
            else:
                prev_size = new_size

        new_cycle = Cycle(joined_set, parents=involved_cycles)
        new_cycles.append(new_cycle)

        remaining_keys = list(node2cycle.keys())

    ## Sanity check (always true)
    # for c_idx, cycle in enumerate(new_cycles):
    #     for other_cycle in new_cycles[c_idx + 1:]:
    #         assert(len(cycle.elements & other_cycle.elements) == 0), "Overlap detected between merged cycles"

    return new_cycles


def resolve_cycle(cycle, lookup_table, lower_cutoff, upper_cutoff, savepath):
    if len(cycle.elements) == 2:
        return {merge_many_nodes(cycle.elements, lookup_table)}

    before = max(len(node.children) for node in cycle.elements)
    pooled_children = set([child for node in cycle.elements for child in node.children]) - cycle.elements
    after = len(pooled_children)

    # control the maximal expansion to prevent explosions:
    if after - before > lower_cutoff:
        children_ratios = compute_children_overlap_ratio(pooled_children, cycle.elements)
        pooled_parents = set([parent for node in cycle.elements for parent in node.parents]) - cycle.elements
        parent_ratios = compute_parent_overlap_ratio(pooled_parents, cycle.elements)

        if max(min(children_ratios), min(parent_ratios)) == 0.0:
            # cannot merge
            broken_cycle = break_up_cycle(cycle.elements, cycle, lookup_table)
            return broken_cycle
        else:
            if after - before > upper_cutoff:
                nodes2png(cycle.elements, cycle.parents, savepath + "_before.pdf")
                broken_cycle = break_up_cycle(cycle.elements, cycle, lookup_table)
                print("excessive merge broken (%d => %d)" % (
                        max([len(n.children) for n in cycle.elements]),
                        max([len(n.children) for n in broken_cycle]))
                )
                nodes2png(broken_cycle, cycle.parents, savepath + "_after.pdf")
                return broken_cycle
            else:
                res = merge_many_nodes(cycle.elements, lookup_table)
                print("weird merge (%d => %d)" %
                      (max([len(n.children) for n in cycle.elements]),
                       len(res.children))
                )
                # allow merge... although this could be bad
                return {res}
    else:
        return {merge_many_nodes(cycle.elements, lookup_table)}


def merge_cyclical_nodes(nodes, lookup_table, activated, min_cycle, max_cycle, lower_cutoff, upper_cutoff, savepath):
    cycles = []
    node2cycle = {}

    merged = set()
    total_search_time = 0.0

    for cycle_size in range(min_cycle, max_cycle + 1):
        tsearch_duration = time.time()
        num_found = 0

        for el in nodes:
            found = find_ncycles(el, cycle_size, activated)
            if len(found) > 0:
                num_found += 1

            for cycle_elements in found:

                cycle = Cycle(set(cycle_elements))
                for key in cycle.elements:
                    if key in node2cycle:
                        node2cycle[key].elements = node2cycle[key].elements | cycle.elements
                    else:
                        node2cycle[key] = cycle
        tsearch_duration = time.time() - tsearch_duration
        total_search_time += tsearch_duration
        print("%d cycles found of size %d (%.2fs)" % (num_found, cycle_size, tsearch_duration))
        if num_found > 0:
            break

    found_cycles = list(set(node2cycle.values()))
    tmerge_cycle_time = time.time()
    found_cycles = merge_cycles(found_cycles)
    tmerge_cycle_time = time.time() - tmerge_cycle_time

    for node in nodes:
        if node not in node2cycle:
            merged.add(node)

    # consider the unique set of cycles
    tmerge_node_time = time.time()
    newnodes = set()

    for cycle in found_cycles:
        assert(len(cycle) > 1), "Cycle (%r) must have more than 1 element" % (cycle,)
        for node in cycle.elements:
            if node.name not in lookup_table:
                assert(node.name in lookup_table), "Could not find %s in lookup_table" % (node.name,)

    for k, cycle in enumerate(found_cycles):
        newnode = resolve_cycle(cycle, lookup_table, lower_cutoff, upper_cutoff, savepath + "_%d" % (k,))
        newnodes = newnodes | newnode

    merged = merged | newnodes
    tmerge_node_time = time.time() - tmerge_node_time
    return merged, newnodes, total_search_time, tmerge_cycle_time, tmerge_node_time

def remove_cycles(lookup_table, min_cycle_size=2, max_cycle_size=10, lower_cutoff=100, upper_cutoff=1000):
    assert(max_cycle_size >= min_cycle_size), "max_cycle_size must be larger than or equal to min_cycle_size"
    leaves = [val for val in lookup_table.values() if len(val.children) == 0]

    visited = set()
    to_process = set(leaves)

    prev_new_candidates = None

    activated = {}
    blocked = False

    iteration = 0

    while len(to_process) > 0:
        print(len(to_process))
        new_candidates = {}
        for candidate in to_process:
            if len(candidate.children) == 0:
                activated[candidate.name] = True
                for parent in candidate.parents:
                    if parent.name not in visited:
                        new_candidates[parent.name] = parent
                        visited.add(parent.name)
            else:
                if all(child.name in activated for child in candidate.children):
                    for parent in candidate.parents:
                        if parent.name not in visited:
                            new_candidates[parent.name] = parent
                            visited.add(parent.name)

                    activated[candidate.name] = True
                else:
                    new_candidates[candidate.name] = candidate


        to_process = set(new_candidates.values())
        if prev_new_candidates == new_candidates:
            savepath = "cycle_conflict_" + str(iteration)
            prev = len(to_process)
            new_to_process, merged_nodes, search_time, merge_cycle_time, merge_node_time = merge_cyclical_nodes(
                to_process,
                lookup_table,
                activated,
                min_cycle_size,
                max_cycle_size,
                lower_cutoff,
                upper_cutoff,
                savepath
            )
            now = len(new_to_process)
            print("Merging cycles: %d nodes converted into %d nodes (search: %.2fs, cycle merge: %.2fs, node merge: %.2fs)" % (
                prev, now, search_time, merge_cycle_time, merge_node_time
            ))
            to_process = new_to_process

            if len(merged_nodes) == 0:

                if not blocked:
                    # partial failure: -- try to recover by
                    # finding cycles in the unactivated and not-yet-queued pre-req elements
                    # of current queued elements:
                    external_disturbance = set()
                    for node in to_process:
                        for child in node.children:
                            if child.name not in activated and child not in to_process:
                                external_disturbance.add(child)

                    if len(external_disturbance) > 0:
                        print("Extending processing queue with %d remaining nodes" % (len(external_disturbance),))
                        to_process = to_process | external_disturbance
                        blocked = True
                    else:
                        break
                else:
                    # attempt to expand original set failed. No cycle or resolution could be made.
                    # stuck at this stage forever. Breaking out.
                    break
        else:
            blocked = False
            prev_new_candidates = new_candidates

        iteration += 1

    return to_process, activated


def findcycle(current, support, n, origin=None):
    if n == 1:
        return False
    if origin is None:
        origin = current

    for node in current.parents:
        if node == origin:
            return True
        else:
            if node in support:
                if findcycle(node, support, n-1, origin):
                    return True

    return False


def cycles_exist(support, max_size=10):
    for node in support:
        if findcycle(node, support, max_size):
            return node
    return None
