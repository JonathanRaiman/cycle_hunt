#include "Cycle.h"

#include <iostream>
#include <iomanip>
#include <fstream>

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

bool edge_in_cycle(const std::vector<int>& cycle, const std::tuple<int, int>& edge) {
    for (int i = 0; i < cycle.size(); i++) {
        if (cycle[i] == std::get<0>(edge) && cycle[(i + 1) % cycle.size()] == std::get<1>(edge)) {
            return true;
        }
    }
    return false;
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

template<typename T>
T&& set_pop(std::unordered_set<T>& set) {
    auto node = (*set.begin());
    set.erase(node);
    return std::move(node);
}

template<typename T>
void remove_elements(std::vector<T>& v, const T& key) {
    v.erase( std::remove( v.begin(), v.end(), key ), v.end() );
}

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


std::vector<std::tuple<int, int>> cycle_edges(const std::vector<int>& cycle) {
    std::vector<std::tuple<int, int>> out;
    for (int i = 0; i < cycle.size(); i++) {
        out.emplace_back(cycle[i], cycle[(i + 1) % cycle.size()]);
    }
    return out;
}

namespace cycle_hunt {

    Graph Graph::copy() const {
        auto connections_copy = connections;
        return Graph(std::move(connections_copy));
    }

    void Graph::add_node(node_t node) {
        if (!is_in(node, connections)) {
            connections[node] = connections_data_t();
        }
    }

    std::vector<typename Graph::node_t> Graph::get_roots() const {
        std::vector<node_t> roots;
        for (auto& kv : connections) {
            if (std::get<0>(kv.second).size() == 0) {
                roots.emplace_back(kv.first);
            }
        }
        return roots;
    }

    std::vector<typename Graph::node_t> Graph::get_leaves() const {
        std::vector<node_t> roots;
        for (auto& kv : connections) {
            if (std::get<1>(kv.second).size() == 0) {
                roots.emplace_back(kv.first);
            }
        }
        return roots;
    }

    void Graph::add_edge(node_t from, node_t to) {
        std::get<0>(connections[from]).push_back(to);
        std::get<1>(connections[to]).push_back(from);
    }

    size_t Graph::size() const {
        return connections.size();
    }

    size_t Graph::number_of_edges() const {
        size_t num_edges = 0;
        for (auto& kv : connections) {
            num_edges += std::get<0>(kv.second).size();
        }
        return num_edges;
    }

    typename Graph::node_t Graph::max_key() const {
        return std::max_element(connections.begin(), connections.end(),
            [](const std::pair<node_t, connections_data_t>& a,
               const std::pair<node_t, connections_data_t>& b) {
            return a.first > b.first;
        })->first;
    }

    const std::vector<typename Graph::node_t>& Graph::neighbors(node_t node) const {
        return std::get<0>(connections.at(node));
    }

    const std::vector<typename Graph::node_t>& Graph::predecessors(node_t node) const {
        return std::get<1>(connections.at(node));
    }

    void Graph::remove_node(node_t node) {
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

    void Graph::remove_edge(node_t from, node_t to) {
        remove_elements(
            std::get<0>(connections[from]),
            to
        );
        remove_elements(
            std::get<1>(connections[to]),
            from
        );
    }

    Graph Graph::copy_edges() const {
        connections_t nonempty_connections;
        for (auto& node : connections) {
            if (std::get<0>(node.second).size() > 0  || std::get<1>(node.second).size() > 0) {
                nonempty_connections[node.first] = node.second;
            }
        }
        return Graph(std::move(nonempty_connections));
    }

    Graph Graph::subgraph(int_set_t& nbunch) const {
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

    void Graph::update_subgraph(const Graph& subgraph) {
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

    Graph::RemoveCyclesResult Graph::greedy_remove_cycles_step(
            const Graph& G,
            int max_cycles,
            std::unordered_set<std::tuple<int, int>>& blacklist,
            bool update_common_edges,
            bool lazy) {

        auto simple_cycles_result = G.simple_cycles(max_cycles, blacklist, lazy);

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

    std::unordered_set<typename Graph::node_t> Graph::reach_blocking_nodes() const {
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

    bool Graph::root_can_reach_all() const {
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

    Graph Graph::greedy_remove_cycles(
            const Graph& original,
            int max_cycles,
            bool update_common_edges) {

        edge_set_t blacklist;
        size_t total_removals = 0;

        Graph G = original;

        bool lazy = true;

        while (true) {
            auto res = greedy_remove_cycles_step(
                G,
                max_cycles,
                blacklist,
                update_common_edges,
                lazy
            );
            G = std::move(res.G);
            std::cout << "res.removals.size() = " << res.removals.size() << std::endl;
            total_removals += res.removals.size();
            if (res.removals.size() == 0) {
                if (!lazy) {
                    break;
                } else {
                    lazy = false;
                }
            }
        }
        std::cout << "total removals size = " << total_removals << std::endl;

        return G;
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
            const std::unordered_set<std::tuple<int, int>>& blacklist,
            bool lazy) const {

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

            if (scc.size() == 1 && lazy) {
                auto startnode = set_pop(scc);
                subG.remove_node(startnode);
            } else {
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
        }
        std::cout << std::endl;

        return result;
    }


    Graph Graph::load_connection_matrix(const std::string& fname) {
        std::ifstream fp(fname, std::ios::in | std::ios::binary);
        int32_t from, to;

        Graph G;

        std::unordered_set<std::tuple<int32_t, int32_t>> visited;

        while (fp) {
            fp.read((char*)&from, sizeof(from));
            fp.read((char*)&to, sizeof(to));
            auto insertable = std::tuple<int32_t, int32_t>(from, to);
            if (visited.find(insertable) == visited.end()) {
                G.add_edge(from, to);
                visited.insert(insertable);
            }
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

}

std::ostream& operator<<(std::ostream& os, const cycle_hunt::Graph& G) {
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

