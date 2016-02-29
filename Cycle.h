#ifndef CYCLE_H
#define CYCLE_H

#include <vector>
#include <unordered_map>
#include <unordered_set>


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


class Graph {

    public:
        struct RemoveCyclesResult;

        typedef int node_t;
        typedef std::tuple<std::vector<node_t>, std::vector<node_t>> connections_data_t;
        typedef std::unordered_map<node_t, connections_data_t> connections_t;

        connections_t connections;

        Graph() {};

        Graph(connections_t&& _connections) : connections(_connections) {};

        Graph copy() const;

        void add_node(node_t node);

        std::vector<node_t> get_roots() const;

        std::vector<node_t> get_leaves() const;

        void add_edge(node_t from, node_t to);

        const std::vector<node_t>& neighbors(node_t node) const;

        const std::vector<node_t>& predecessors(node_t node) const;

        Graph copy_edges() const;

        void remove_node(node_t node);

        void remove_edge(node_t from, node_t to);

        Graph subgraph(std::unordered_set<node_t>& nbunch) const;

        void update_subgraph(const Graph& subgraph);

        size_t size() const;

        size_t number_of_edges() const;

        node_t max_key() const;

        std::unordered_set<node_t> reach_blocking_nodes() const;

        bool root_can_reach_all() const;


        static RemoveCyclesResult greedy_remove_cycles_step(
                const Graph& G,
                int max_cycles,
                std::unordered_set<std::tuple<int, int>>& blacklist,
                bool update_common_edges,
                bool lazy);

        static Graph greedy_remove_cycles(
                const Graph& G,
                int max_cycles,
                bool update_common_edges);

        struct SimpleCyclesResult {
            std::vector<std::vector<int> > paths;
            std::unordered_map<std::tuple<int, int>, int> common_edges;
        };

        std::vector<std::unordered_set<node_t>> strongly_connected_components() const;

        SimpleCyclesResult simple_cycles(
                int max_cycles,
                const std::unordered_set<std::tuple<int, int>>& blacklist,
                bool lazy) const;

        void save_connection_matrix(const std::string& fname) const;
        static Graph load_connection_matrix(const std::string& fname);

};

struct Graph::RemoveCyclesResult {
    Graph G;
    std::vector<std::tuple<int, int>> removals;
};

#endif
