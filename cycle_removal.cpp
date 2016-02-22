#include <iostream>
#include <vector>
#include <iomanip>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <string>

// Printing utilities
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::unordered_set<T>& v) {
        if (v.size() == 0) return os << "{}";
        os << "{";
        for (auto& f : v)
                os << std::fixed
                   << std::setw( 7 ) // keep 7 digits
                   << std::setprecision( 3 ) // use 3 decimals
                   << std::setfill( ' ' ) // pad values with blanks this->w(i,j)
                   << f << ", ";
        return os << "}";
}


template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
        if (v.size() == 0) return os << "[]";
        os << "[";
        for (auto& f : v)
                os << std::fixed
                   << std::setw( 7 ) // keep 7 digits
                   << std::setprecision( 3 ) // use 3 decimals
                   << std::setfill( ' ' ) // pad values with blanks this->w(i,j)
                   << f << " ";
        return os << "]";
}

using std::vector;

typedef std::unordered_set<int> int_set_t;
typedef std::unordered_map<int, int> int_map_t;

template<typename T>
bool is_in(const int& key, const std::unordered_map<int, T>& map) {
    return map.find(key) != map.end();
}

bool is_in(const int& key, const std::unordered_set<int>& map) {
    return map.find(key) != map.end();
}

template<typename T>
void remove_elements(std::vector<T>& v, const T& key) {
    v.erase( std::remove( v.begin(), v.end(), key ), v.end() );
}

class Graph {

    public:
        typedef int node_t;

        typedef std::tuple<std::vector<node_t>, std::vector<node_t>> connections_data_t;

        typedef std::unordered_map<node_t, connections_data_t> connections_t;

        connections_t connections;

        Graph() {};

        Graph(connections_t&& _connections) : connections(_connections) {};

        void add_node(node_t node) {
            if (!is_in(node, connections)) {
                connections[node] = connections_data_t();
            }
        }

        void add_edge(node_t from, node_t to) {
            add_node(from);
            add_node(to);
            std::get<0>(connections[from]).push_back(to);
            std::get<1>(connections[to]).push_back(from);
        }

        const std::vector<node_t>& neighbors(node_t node) const {
            return std::get<0>(connections.at(node));
        }

        Graph copy_edges() const {
            connections_t nonempty_connections;
            for (auto& node : connections) {
                if (std::get<0>(node.second).size() > 0  || std::get<1>(node.second).size() > 0) {
                    nonempty_connections[node.first] = node.second;
                }
            }
            return Graph(std::move(nonempty_connections));
        }

        void remove_node(node_t node) {
            if (is_in(node, connections)) {
                auto& node_data = connections[node];

                for (auto& other_node : std::get<0>(node_data)) {
                    remove_elements(
                        std::get<1>(connections[other_node]),
                        node
                    );
                }

                for (auto& other_node : std::get<1>(node_data)) {
                    remove_elements(
                        std::get<0>(connections[other_node]),
                        node
                    );
                }

            }
        }

        Graph subgraph(int_set_t& node_set) {
            return Graph();
        }
};

std::vector<int_set_t> strongly_connected_components(const Graph& G) {
    int_map_t preorder;
    int_map_t lowlink;
    std::unordered_map<int, bool> scc_found;
    std::vector<int> scc_queue;

    int i = 0; // Preorder counter

    std::vector<int_set_t> scc_components;

    for (auto& it : G.connections) {
        auto& source = it.first;
        if (!is_in(source, scc_found)) {
            std::vector<int> queue = {source};

            while (!queue.empty()) {
                auto v = queue.back();

                if (!is_in(v, preorder)) {
                    i = i + 1;
                    preorder[v] = i;
                }
                bool done = true;
                auto& v_nbrs = G.neighbors(v);

                for (auto& w : v_nbrs) {
                    if (!is_in(w, preorder)) {
                        queue.push_back(w);
                        done = false;
                        break;
                    }
                }
                if (done) {
                    lowlink[v] = preorder[v];
                    for (auto& w : v_nbrs) {
                        if (!is_in(w, scc_found)) {
                            if (preorder[w] > preorder[v]) {
                                lowlink[v] = std::min(lowlink[v], lowlink[w]);
                            } else {
                                lowlink[v] = std::min(lowlink[v], preorder[w]);
                            }
                        }
                    }
                    queue.pop_back();
                    if (lowlink[v] == preorder[v]) {
                        scc_found[v] = true;
                        int_set_t scc = {v};

                        while (!scc_queue.empty() && preorder[scc_queue.back()] > preorder[v]) {
                            int k = scc_queue.back();
                            scc_queue.pop_back();
                            scc_found[k] = true;
                            scc.insert(k);
                        }
                        scc_components.emplace_back(scc);
                    } else {
                        scc_queue.push_back(v);
                    }
                }
            }
        }
    }
    return scc_components;
}

template<typename T>
T&& set_pop(std::unordered_set<T>& set) {
    auto node = (*set.begin());
    set.erase(node);
    return std::move(node);
}

std::vector<std::vector<int> > simple_cycles(const Graph& G) {

    auto _unblock = [](const int& thisnode, int_set_t& blocked, std::unordered_map<int, int_set_t>& B) {
        int_set_t stack = {thisnode};
        while (!stack.empty()) {
            // pop from stack
            auto node = set_pop(stack);

            if (is_in(node, blocked)) {
                blocked.erase(node);
                for (auto& v : B[node]) {
                    stack.insert(v);
                }
                B[node].clear();
            }

        }
    };

    Graph subG = G.copy_edges();

    auto sccs = strongly_connected_components(subG);


    std::vector<int> nbrs;

    int iters = 0;

    std::vector<std::vector<int>> paths;
    std::cout << "sccs begin" << std::endl;
    while (!sccs.empty()) {
        std::cout << "sccs iter" << std::endl;
        auto scc = sccs.back();
        sccs.pop_back();

        auto startnode = set_pop(scc);

        std::vector<int> path = {startnode};

        int_set_t blocked;
        int_set_t closed;

        blocked.insert(startnode);

        std::unordered_map<int, int_set_t> B;

        std::vector<std::tuple<int, std::vector<int>> > stack;

        stack.emplace_back(
            startnode,
            subG.neighbors(startnode)
        );

        std::cout << "stack begin" << std::endl;
        iters = 0;
        while (!stack.empty()) {
            iters += 1;
            std::cout << "    stack iter" << std::endl;

            auto& stack_back = stack.back();
            auto& nbrs = std::get<1>(stack_back);
            const auto& thisnode = std::get<0>(stack_back);

            if (!nbrs.empty()) {
                auto nextnode = nbrs.back();
                nbrs.pop_back();

                if (nextnode == startnode) {
                    paths.emplace_back(path);

                    for (auto& v : path) {
                        closed.insert(v);
                    }
                } else if (!is_in(nextnode, blocked)) {
                    path.push_back(nextnode);
                    stack.emplace_back(
                        nextnode,
                        subG.neighbors(nextnode)
                    );
                    closed.erase(nextnode);
                    blocked.insert(nextnode);
                    continue;
                }
            }

            if (nbrs.empty()) {
                std::cout << "nbrs empty (stack len = " <<  stack.size() << ")" << std::endl;
                if (is_in(thisnode, closed)) {
                    _unblock(thisnode, blocked, B);
                } else {
                    for (auto& nbr : subG.neighbors(thisnode)) {
                        if (!is_in(thisnode, B[nbr])) {
                            B[nbr].insert(thisnode);
                        }
                    }
                }
                stack.pop_back();

                path.pop_back();
            }
            std::cout << "    stack len = " << stack.size() << std::endl;
            if (iters > 10) {
                return paths;
            }
        }
        std::cout << "stack end" << std::endl;

        subG.remove_node(startnode);

        auto H = subG.subgraph(scc);

        for (auto& h_scc : strongly_connected_components(H)) {
            sccs.push_back(h_scc);
        }


    }
    std::cout << "sccs end" << std::endl;

    return paths;
}

int main() {
    auto G = Graph();

    G.add_node(0);
    G.add_node(1);
    G.add_node(3);
    G.add_edge(0, 1);
    G.add_edge(1, 3);
    G.add_edge(3, 0);

    std::cout << strongly_connected_components(G) << std::endl;
    std::cout << simple_cycles(G) << std::endl;
}
