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

class Graph {

    public:
        std::unordered_map<int, std::vector<int>> connections;

        Graph() {};

        void add_node(int node) {
            if (!is_in(node, connections)) {
                connections[node] = {};
            }
        }

        void add_edge(int from, int to) {
            add_node(from);
            add_node(to);
            connections[from].push_back(to);
        }

        const std::vector<int>& neighbors(int node) const {
            return connections.at(node);
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

int main() {
    auto G = Graph();

    G.add_node(0);
    G.add_node(1);
    G.add_node(3);
    G.add_edge(0, 1);
    G.add_edge(1, 3);
    G.add_edge(3, 0);

    std::cout << strongly_connected_components(G) << std::endl;
}
