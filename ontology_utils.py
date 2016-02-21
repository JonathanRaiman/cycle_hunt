import time
import gzip
import os
import io

from IPython.display import clear_output

from node_utils import OntologyBranch
from cycle_utils import remove_cycles


def add_lattice_edge(parent, child, branch_map, parentless):
    if type(child) is str:
        if child not in branch_map:
            branch_map[child] = OntologyBranch(child)
        child = branch_map[child]
    if type(parent) is str:
        if parent not in branch_map:
            branch_map[parent] = OntologyBranch(parent)
            parentless.append(branch_map[parent])
        parent = branch_map[parent]
    if child != parent:
        child.add_parent(parent)
    return (parent, child)


def load_roots_from_stream(fp, roots, total_size):
    text_fin = io.TextIOWrapper(fp, newline='')

    parentless = []
    branch_map = {}
    right_arrow = "->"
    left_arrow = "<-"

    marked_branch = None
    last_edge_is_right_arrow = True

    for k, line in enumerate(text_fin):
        tokens = line.split(right_arrow, 1)
        if len(tokens) >= 2:
            for i in range(len(tokens)-1):
                marked_branch = add_lattice_edge(tokens[i], tokens[i+1].strip(), branch_map, parentless)[0]
                last_edge_is_right_arrow = True
        else:
            tokens = line.split(left_arrow, 1)
            if len(tokens) >= 2:
                for i in range(len(tokens)-1):
                    marked_branch = add_lattice_edge(tokens[i+1].strip(), tokens[i], branch_map, parentless)[1]
                    last_edge_is_right_arrow = False
            elif marked_branch is not None:
                if last_edge_is_right_arrow:
                    add_lattice_edge(marked_branch, tokens[0].strip(), branch_map, parentless)
                else:
                    add_lattice_edge(tokens[0].strip(), marked_branch, branch_map, parentless)
        if k % 2000 == 0:
            progress = fp.tell() / total_size
            print("█" * (int(20 * progress)) + " %.1f%%" % (100 * progress,))
            clear_output(wait=True)

    for k in parentless:
        if len(k.parents) == 0:
            roots.append(k)
            k.lookup_table = branch_map


def load_abstract_trees(path):
    roots = []
    total_size = os.stat(path).st_size

    try:
        fp = gzip.open(path, "rb")
        load_roots_from_stream(fp, roots, total_size)
        fp.close()
    except OSError:
        fpalt = open(path, "rb")
        load_roots_from_stream(fpalt, roots, total_size)
        fpalt.close()
    finally:
        fp.close()

    for key in roots[0].lookup_table.values():
        if key in key.parents:
            print("p")
        if key in key.children:
            print("c")
    #print(len(roots[0].lookup_table))
    #to_process, activated = remove_cycles(roots[0].lookup_table)
    #print(len(roots[0].lookup_table))
    #assert(len(to_process) == 0), "Some cycles remain in graph"
    #
    return roots
