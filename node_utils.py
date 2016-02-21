class OntologyBranch(object):
    __slots__ = ["parents", "children", "name", "lookup_table", "id", "history"]

    def __init__(self, name):
        self.name = name
        self.parents = []
        self.children = []
        self.history = None
        self.lookup_table = None

    def add_parent(self, parent):
        self.parents.append(parent)
        parent.add_child(self)

    def revert(self, lookup_table):
        if self.history is not None and len(self.history) > 0:
            return revert_node(self, self.history, lookup_table)
        else:
            raise ValueError("This branch's history is empty")

    def add_child(self, child):
        self.children.append(child)

    def __str__(self):
        return "<Branch name=\"%s\", num_children=%d, num_parents=%d>" % (self.name, len(self.children), len(self.parents))

    def __repr__(self):
        return str(self)

    def deep_copy(self):
        assert(self.lookup_table is not None)
        newlookup = deepcopy_tree(self.lookup_table)
        return newlookup[self.name]


def deepcopy_tree(lookup_table):
    copy_table = {}

    for key in lookup_table:
        alternate = OntologyBranch(key)
        copy_table[key] = alternate

    for key, value in lookup_table.items():
        alternate = copy_table[key]
        alternate.parents = [copy_table[parent.name]
                             for parent in value.parents]
        alternate.children = [copy_table[child.name]
                             for child in value.children]

        if value.lookup_table is not None:
            alternate.lookup_table = copy_table

    return copy_table



def compute_children_overlap_ratio(pool, nodes):
    ratios = []
    for node in nodes:
        num = len(set(node.children) & pool)
        ratio = num / len(node.children)
        ratios.append(ratio)
    return ratios


def compute_parent_overlap_ratio(pool, nodes):
    ratios = []
    for node in nodes:
        num = len(set(node.parents) & pool)
        ratio = num / len(node.parents)
        ratios.append(ratio)
    return ratios


def disconnect_node(node):
    for parent in node.parents:
        parent.children = list(set(parent.children) - {node})

    for child in node.children:
        child.parents = list(set(child.parents) - {node})


def reconnect_node(node):
    for parent in node.parents:
        parent.children = list(set(parent.children) | {node})

    for child in node.children:
        child.parents = list(set(child.parents) | {node})


def revert_node(node, other_nodes, lookup_table):
    disconnect_node(node)

    for other_node in other_nodes:
        reconnect_node(other_node)
        lookup_table[other_node.name] = other_node

    del lookup_table[node.name]
    return other_nodes

