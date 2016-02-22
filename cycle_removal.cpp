#include <iostream>
#include <vector>
#include <iomanip>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <fstream>
#include <sys/stat.h>

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


template<typename A, typename B>
std::ostream& operator<<(std::ostream& os, const std::tuple<A, B>& v) {
        os << "(";
        os << std::fixed
           << std::setw( 7 ) // keep 7 digits
           << std::setprecision( 3 ) // use 3 decimals
           << std::setfill( ' ' ) // pad values with blanks this->w(i,j)
           << std::get<0>(v) << ", ";
        os << std::fixed
           << std::setw( 7 ) // keep 7 digits
           << std::setprecision( 3 ) // use 3 decimals
           << std::setfill( ' ' ) // pad values with blanks this->w(i,j)
           << std::get<1>(v);
        return os << ")";
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

template<typename K, typename T>
bool is_in(const K& key, const std::unordered_map<K, T>& map) {
    return map.find(key) != map.end();
}

template<typename K>
bool is_in(const K& key, const std::unordered_set<K>& map) {
    return map.find(key) != map.end();
}

template<typename T>
void remove_elements(std::vector<T>& v, const T& key) {
    v.erase( std::remove( v.begin(), v.end(), key ), v.end() );
}

namespace std {
    template<typename... TTypes>
    class hash<std::tuple<TTypes...>> {
        private:
            typedef std::tuple<TTypes...> Tuple;

            template<int N>
            size_t operator()(Tuple value) const {
                return 0;
            }

            template<int N, typename THead, typename... TTail>
            size_t operator()(Tuple value) const {
                constexpr int Index = N - sizeof...(TTail) - 1;
                return hash<THead>()(std::get<Index>(value)) ^ operator()<N, TTail...>(value);
            }

        public:
            size_t operator()(Tuple value) const {
                return operator()<sizeof...(TTypes), TTypes...>(value);
            }
    };
}

template<typename K, typename condition_t>
std::vector<K> copy_if(const std::vector<K>& from, condition_t& condition) {
    std::vector<K> out;
    for (auto& n : from) {
        if (condition(n)) {
            out.emplace_back(n);
        }
    }
    return out;
}


class Graph {

    public:
        typedef int node_t;
        typedef std::tuple<std::vector<node_t>, std::vector<node_t>> connections_data_t;
        typedef std::unordered_map<node_t, connections_data_t> connections_t;

        connections_t connections;

        Graph() {};

        Graph(connections_t&& _connections) : connections(_connections) {};

        Graph copy() const {
            auto connections_copy = connections;
            return Graph(std::move(connections_copy));
        }

        void add_node(node_t node) {
            if (!is_in(node, connections)) {
                connections[node] = connections_data_t();
            }
        }

        std::vector<node_t> get_roots() const {
            std::vector<node_t> roots;
            for (auto& kv : connections) {
                if (std::get<0>(kv.second).size() == 0) {
                    roots.emplace_back(kv.first);
                }
            }
            return roots;
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

        const std::vector<node_t>& predecessors(node_t node) const {
            return std::get<1>(connections.at(node));
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

        void remove_edge(node_t from, node_t to) {
            remove_elements(
                std::get<0>(connections[from]),
                to
            );
            remove_elements(
                std::get<1>(connections[to]),
                from
            );
        }

        Graph subgraph(int_set_t& nbunch) const {
            Graph H;

            auto condition = [&nbunch](const node_t& node) {
                return is_in(node, nbunch);
            };

            for (auto& n : nbunch) {
                if (is_in(n, connections)) {
                    H.connections[n] = connections_data_t(
                        copy_if(std::get<0>(connections.at(n)), condition),
                        copy_if(std::get<1>(connections.at(n)), condition)
                    );
                }
            }
            return H;
        }

        size_t size() const {
            return connections.size();
        }

        node_t max_key() const {
            return std::max_element(connections.begin(), connections.end(),
                [](const std::pair<node_t, connections_data_t>& a,
                   const std::pair<node_t, connections_data_t>& b) {
                return a.first > b.first;
            })->first;
        }

};

std::ostream& operator<<(std::ostream& os, const Graph& G) {
        if (G.size() == 0) return os << "Graph()";
        os << "Graph(\n";
        for (auto& n : G.connections) {
                os << std::fixed
                   << std::setw( 7 ) // keep 7 digits
                   << std::setprecision( 3 ) // use 3 decimals
                   << std::setfill( ' ' ) // pad values with blanks this->w(i,j)
                   << n.first << " -> " << std::get<0>(n.second) << "\n";
        }
        return os << ")";
}

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

struct SimpleCyclesResult {
    std::vector<std::vector<int> > paths;
    std::unordered_map<std::tuple<int, int>, int> common_edges;
};

SimpleCyclesResult simple_cycles(const Graph& G, int max_cycles, const std::unordered_set<std::tuple<int, int>>& blacklist) {

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
    int max_edge = 0;

    SimpleCyclesResult result;
    auto& paths = result.paths;
    auto& common_edges = result.common_edges;

    while (!sccs.empty()) {
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

        while (!stack.empty()) {

            auto& stack_back = stack.back();
            auto& nbrs = std::get<1>(stack_back);
            const auto& thisnode = std::get<0>(stack_back);

            if (!nbrs.empty()) {
                auto nextnode = nbrs.back();
                nbrs.pop_back();

                if (nextnode == startnode) {
                    paths.emplace_back(path);

                    for (int i = 0; i < path.size(); i++) {
                        std::tuple<int, int> edge(path[i], path[(i+1) % path.size()]);

                        if (is_in(edge, common_edges)) {
                            common_edges[edge] += 1;
                            max_edge = std::max(common_edges[edge], max_edge);
                        } else {
                            if (!is_in(edge, blacklist)) {
                                common_edges[edge] = 1;
                            }
                        }
                    }

                    if (paths.size() >= max_cycles) {
                        std::cout << std::endl;
                        return result;
                    }

                    if (paths.size() % 10 == 0) {
                        std::cout << paths.size() << " paths found. "
                                  << common_edges.size() << " common_edges. "
                                  << max_edge << " max_edge\r" << std::flush;
                    }

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
        }

        subG.remove_node(startnode);

        auto H = subG.subgraph(scc);
        for (auto& h_scc : strongly_connected_components(H)) {
            sccs.push_back(h_scc);
        }
    }
    std::cout << std::endl;

    return result;
}


Graph load_connection_matrix(const std::string& fname) {
    std::ifstream fp(fname, std::ios::in | std::ios::binary);
    int32_t from, to;

    Graph G;

    while (fp) {
        fp.read((char*)&from, sizeof(from));
        fp.read((char*)&to, sizeof(to));
        G.add_edge(from, to);
    }

    return G;
}

bool file_exists (const std::string& fname) {
    struct stat buffer;
    return (stat(fname.c_str(), &buffer) == 0);
}

std::vector<std::tuple<int, int>> cycle_edges(const std::vector<int>& cycle) {
    std::vector<std::tuple<int, int>> out;
    for (int i = 0; i < cycle.size(); i++) {
        out.emplace_back(cycle[i], cycle[(i + 1) % cycle.size()]);
    }
    return out;
}

bool edge_in_cycle(const std::vector<int>& cycle, const std::tuple<int, int>& edge) {
    for (int i = 0; i < cycle.size(); i++) {
        if (cycle[i] == std::get<0>(edge) && cycle[(i + 1) % cycle.size()] == std::get<1>(edge)) {
            return true;
        }
    }
    return false;
}

struct RemoveCyclesResult {
    Graph G;
    std::vector<std::tuple<int, int>> removals;
};


bool root_can_reach_all(const Graph& G) {
    auto root = G.get_roots()[0];
    int_set_t visited = {root};
    std::vector<int> queue = {root};

    while (!queue.empty()) {
        std::vector<int> new_queue;

        for (auto& node : queue) {
            for (auto& neighbor : G.predecessors(node)) {
                if (!is_in(neighbor, visited)) {
                    visited.insert(neighbor);
                    new_queue.emplace_back(neighbor);
                }
            }
        }
        queue = std::move(new_queue);
    }

    if (visited.size() != G.size()) {
        return false;
    } else {
        return true;
    }
}


RemoveCyclesResult greedy_remove_cycles(const Graph& G, int max_cycles, std::unordered_set<std::tuple<int, int>>& blacklist, bool update_common_edges) {
    auto simple_cycles_result = simple_cycles(G, max_cycles, blacklist);

    auto& cycles = simple_cycles_result.paths;
    auto& common_edges = simple_cycles_result.common_edges;

    RemoveCyclesResult res;
    res.G = G;

    auto& removals = res.removals;

    while (cycles.size() > 0) {
        auto G_copy = res.G.copy();
        auto edge_kv = std::max_element(
            common_edges.begin(), common_edges.end(),
            [](const std::pair<std::tuple<int, int>, int>& a,
               const std::pair<std::tuple<int, int>, int>& b) {
                return a.second < b.second;
            });

        auto edge = edge_kv->first;

        G_copy.remove_edge(std::get<0>(edge), std::get<1>(edge));

        auto is_connected = root_can_reach_all(G_copy);

        common_edges.erase(edge);

        if (is_connected) {
            removals.push_back(edge);
            res.G = std::move(G_copy);

            std::vector<std::vector<int>> out_cycles;
            int in_cycle = 0;
            for (auto& cycle : cycles) {
                if (!edge_in_cycle(cycle, edge)) {
                    out_cycles.emplace_back(cycle);
                } else {
                    in_cycle += 1;
                }
            }

            if (update_common_edges) {
                if (out_cycles.size() > 0) {
                    common_edges.clear();
                    for (auto& cycle : out_cycles) {
                        for (auto& other_edge : cycle_edges(cycle)) {
                            if (!is_in(other_edge, blacklist)) {
                                common_edges[other_edge] += 1;
                            }
                        }
                    }
                }
            }

            cycles = std::move(out_cycles);
            std::cout << cycles.size() << " cycles left (" << in_cycle << ")" << std::endl;
        } else {
            blacklist.insert(edge);
        }

    }

    return res;
}


int main(int argc, char* argv[]) {
    int max_cycles = 5000000;
    bool update_common_edges = false;
    typedef std::unordered_set<std::tuple<int, int>> edge_set_t;

    if (argc < 2) {
        auto G = Graph();

        G.add_node(0);
        G.add_node(1);
        G.add_node(3);
        G.add_edge(0, 1);
        G.add_edge(1, 3);
        G.add_edge(3, 0);

        std::cout << G << std::endl;
        std::cout << "strongly_connected_components => " << strongly_connected_components(G) << std::endl;
        auto cycles = simple_cycles(G, max_cycles, edge_set_t());
        std::cout << "simple_cycles => " << cycles.paths << std::endl;
    } else {
        std::string path(argv[1]);

        if (argc > 2) {
            max_cycles = std::atoi(argv[2]);
        }

        if (argc > 3) {
            update_common_edges = std::atoi(argv[3]);
        }


        if (file_exists(path)) {
            auto G = load_connection_matrix(path);
            std::cout << "loaded graph with size " << G.size() << std::endl;
            std::cout << "strongly_connected_components len = " << strongly_connected_components(G).size() << std::endl;

            edge_set_t blacklist;

            while (true) {
                auto res = greedy_remove_cycles(G, max_cycles, blacklist, update_common_edges);
                G = std::move(res.G);
                std::cout << "res.removals.size() = " << res.removals.size() << std::endl;
                if (res.removals.size() == 0) {
                    break;
                }
            }

        } else {
            std::cout << path << ": No such file or directory" << std::endl;
        }
    }
}
