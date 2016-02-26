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


template<typename T>
T&& set_pop(std::unordered_set<T>& set) {
    auto node = (*set.begin());
    set.erase(node);
    return std::move(node);
}

using std::vector;

typedef std::unordered_set<int> int_set_t;
typedef std::unordered_map<int, int> int_map_t;
typedef std::unordered_set<std::tuple<int, int>> edge_set_t;

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

        std::vector<node_t> get_leaves() const {
            std::vector<node_t> roots;
            for (auto& kv : connections) {
                if (std::get<1>(kv.second).size() == 0) {
                    roots.emplace_back(kv.first);
                }
            }
            return roots;
        }

        void add_edge(node_t from, node_t to) {
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

        void update_subgraph(const Graph& subgraph) {
            auto condition = [&subgraph](const node_t& node) {
                return !is_in(node, subgraph.connections);
            };

            for (const auto& kv : subgraph.connections) {
                {
                    auto new_neighbors = copy_if(
                        std::get<0>(connections.at(kv.first)),
                        condition
                    );

                    new_neighbors.insert(
                        new_neighbors.end(),
                        std::get<0>(kv.second).begin(),
                        std::get<0>(kv.second).end()
                    );

                    std::get<0>(connections[kv.first]) = new_neighbors;
                }
                {
                    auto new_predecessors = copy_if(
                        std::get<1>(connections.at(kv.first)),
                        condition
                    );

                    new_predecessors.insert(
                        new_predecessors.end(),
                        std::get<1>(kv.second).begin(),
                        std::get<1>(kv.second).end()
                    );

                    std::get<1>(connections[kv.first]) = new_predecessors;
                }
            }
        }

        size_t size() const {
            return connections.size();
        }

        size_t number_of_edges() const {
            size_t num_edges = 0;
            for (auto& kv : connections) {
                num_edges += std::get<0>(kv.second).size();
            }
            return num_edges;
        }

        node_t max_key() const {
            return std::max_element(connections.begin(), connections.end(),
                [](const std::pair<node_t, connections_data_t>& a,
                   const std::pair<node_t, connections_data_t>& b) {
                return a.first > b.first;
            })->first;
        }

        std::unordered_set<node_t> reach_blocking_nodes() const {
            std::unordered_set<node_t> to_process;

            {
                auto leaves = get_leaves();
                to_process = std::unordered_set<node_t>(leaves.begin(), leaves.end());
            }

            std::unordered_set<node_t> visited;
            std::unordered_set<node_t> activated;

            while (!to_process.empty()) {
                std::cout << to_process.size() << std::endl;

                decltype(to_process) new_candidates;
                size_t num_postponed = 0;
                size_t queue_size = to_process.size();
                for (const auto& candidate : to_process) {
                    const auto& children = predecessors(candidate);
                    bool all_children_active = false;
                    if (children.size() == 0) {
                        all_children_active = true;
                    } else {
                        all_children_active = std::all_of(
                            children.begin(),
                            children.end(),
                            [&activated](const node_t& child) {
                                return is_in(child, activated);
                            }
                        );
                    }
                    if (all_children_active) {
                        activated.insert(candidate);
                        for (const auto& parent : neighbors(candidate)) {
                            if (!is_in(parent, visited)) {
                                new_candidates.insert(parent);
                                visited.insert(parent);
                            }
                        }
                    } else {
                        num_postponed += 1;
                        new_candidates.insert(candidate);
                    }
                }

                to_process = std::move(new_candidates);
                if (num_postponed == queue_size) {
                    std::cout << "Extending processing queue" << std::endl;
                    // stuck condition
                    for (auto& kv : connections) {
                        if (!is_in(kv.first, activated)) {
                            to_process.insert(kv.first);
                        }
                    }
                    std::cout << queue_size << " => " << to_process.size() << std::endl;
                    break;
                }
            }
            return to_process;
        }

        bool root_can_reach_all() const {
            auto root = get_roots()[0];
            int_set_t visited = {root};
            std::vector<int> queue = {root};

            while (!queue.empty()) {
                std::vector<int> new_queue;

                for (auto& node : queue) {
                    for (auto& neighbor : predecessors(node)) {
                        if (!is_in(neighbor, visited)) {
                            visited.insert(neighbor);
                            new_queue.emplace_back(neighbor);
                        }
                    }
                }
                queue = std::move(new_queue);
            }

            if (visited.size() != this->size()) {
                return false;
            } else {
                return true;
            }
        }

        struct RemoveCyclesResult;

        static RemoveCyclesResult greedy_remove_cycles_step(
                const Graph& G,
                int max_cycles,
                std::unordered_set<std::tuple<int, int>>& blacklist,
                bool update_common_edges);

        static Graph greedy_remove_cycles(
                const Graph& G,
                int max_cycles,
                bool update_common_edges);

        struct SimpleCyclesResult {
            std::vector<std::vector<int> > paths;
            std::unordered_map<std::tuple<int, int>, int> common_edges;
        };

        std::vector<int_set_t> strongly_connected_components() const;

        SimpleCyclesResult simple_cycles(
                int max_cycles,
                const std::unordered_set<std::tuple<int, int>>& blacklist) const;

        void save_connection_matrix(const std::string& fname) const;
        static Graph load_connection_matrix(const std::string& fname);

};

struct Graph::RemoveCyclesResult {
    Graph G;
    std::vector<std::tuple<int, int>> removals;
};

Graph::RemoveCyclesResult Graph::greedy_remove_cycles_step(
        const Graph& G,
        int max_cycles,
        std::unordered_set<std::tuple<int, int>>& blacklist,
        bool update_common_edges) {

    auto simple_cycles_result = G.simple_cycles(max_cycles, blacklist);

    auto& cycles = simple_cycles_result.paths;
    auto& common_edges = simple_cycles_result.common_edges;

    RemoveCyclesResult res;
    res.G = G;

    auto& removals = res.removals;

    std::cout << "found cycles " << cycles.size() << std::endl;

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

        auto is_connected = G_copy.root_can_reach_all();

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
            std::cout << cycles.size() << " cycles left (" << in_cycle << " removed)" << std::endl;
        } else {
            blacklist.insert(edge);
        }

    }

    return res;
}

Graph Graph::greedy_remove_cycles(
        const Graph& original,
        int max_cycles,
        bool update_common_edges) {

    edge_set_t blacklist;
    size_t total_removals = 0;

    Graph G = original;

    while (true) {
        auto res = greedy_remove_cycles_step(
            G,
            max_cycles,
            blacklist,
            update_common_edges
        );
        G = std::move(res.G);
        std::cout << "res.removals.size() = " << res.removals.size() << std::endl;
        total_removals += res.removals.size();
        if (res.removals.size() == 0) {
            break;
        }
    }
    std::cout << "total removals size = " << total_removals << std::endl;

    return G;
}

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

std::vector<int_set_t> Graph::strongly_connected_components() const {
    int_map_t preorder;
    int_map_t lowlink;
    std::unordered_map<int, bool> scc_found;
    std::vector<int> scc_queue;

    int i = 0; // Preorder counter

    std::vector<int_set_t> scc_components;

    for (auto& it : connections) {
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
                auto& v_nbrs = neighbors(v);

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

Graph::SimpleCyclesResult Graph::simple_cycles(
        int max_cycles,
        const std::unordered_set<std::tuple<int, int>>& blacklist) const {

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

    Graph subG = copy_edges();
    auto sccs = subG.strongly_connected_components();
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
        for (auto& h_scc : H.strongly_connected_components()) {
            sccs.push_back(h_scc);
        }
    }
    std::cout << std::endl;

    return result;
}


Graph Graph::load_connection_matrix(const std::string& fname) {
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

void Graph::save_connection_matrix(const std::string& fname) const {
    std::ofstream fp(fname, std::ios::out | std::ios::binary);

    for (auto& kv : connections) {
        for (auto& parent : std::get<0>(kv.second)) {
            fp.write((char*)&kv.first, sizeof(node_t));
            fp.write((char*)&parent, sizeof(node_t));
        }
    }
}

bool file_exists (const std::string& fname) {
    struct stat buffer;
    return (stat(fname.c_str(), &buffer) == 0);
}


int main(int argc, char* argv[]) {
    int max_cycles = 100;
    bool update_common_edges = true;
    std::string save_location = "";

    if (argc < 2) {
        std::cout << "Usage: "
                  << argv[0]
                  << " <load_location> <max_cycles> <update_common_edges>"
                  " <save_location>"
                  << std::endl << std::endl
                  << "load_location <string> : where to load the graph from."
                  << std::endl
                  << "max_cycles <int> : how many simple cycles "
                     "to find during each removal loop"
                  << std::endl
                  << "update_common_edges {0,1} : on each edge removal"
                  ", update the usage count for all remaining edges"
                  << std::endl
                  << "save_location <string> : where should the "
                  "cycle-removed graph be stored"
                  << std::endl;
    } else {
        std::string path(argv[1]);

        if (argc > 2) {
            max_cycles = std::atoi(argv[2]);
        }

        if (argc > 3) {
            update_common_edges = std::atoi(argv[3]);
        }

        if (argc > 4) {
            save_location = std::string(argv[4]);
        }

        if (file_exists(path)) {
            auto G = Graph::load_connection_matrix(path);
            std::cout << "loaded graph with size " << G.size() << std::endl;

            while (true) {
                // iterate through nodes, ensuring their
                // predecessors were visited beforehand. Stop
                // when no new node can be visited due to
                // inactive predecessors. Record those blocking
                // nodes:
                auto bnodes = G.reach_blocking_nodes();
                std::cout << "bnodes.size() = " << bnodes.size() << std::endl;

                if (bnodes.size() == 0) {
                    break;
                }
                // Construct a graph made from those blocking nodes: that's
                // where the cycles are hiding
                auto blockingG = G.subgraph(bnodes);

                // Now perform a greedy removal of cycles by
                // find `max_cycles` and tallying up the
                // edges used in all found cycles, and removing
                // the most common edges until all cycles are
                // gone. Process is repeated until no cycles
                // can be found in the subgraph
                auto unblockedG = Graph::greedy_remove_cycles(
                    blockingG,
                    max_cycles,
                    update_common_edges
                );

                // The newly rewired subgraph can now be
                // reinjected into the original graph. All
                // connections that were not internal to the
                // subgraph are kept, while all connections
                // internal to the subgraph are replaced by
                // the cycle-removed subgraph.
                G.update_subgraph(unblockedG);
            }

            std::cout << "G.number_of_edges() = " << G.number_of_edges() << std::endl;

            if (!save_location.empty()) {
                // The graph can now be saved in the same binary format as before
                std::cout << "saving uncycled graph to \"" << save_location << "\"" << std::endl;
                G.save_connection_matrix(save_location);
            }
        } else {
            std::cout << path << ": No such file or directory" << std::endl;
        }
    }
}
