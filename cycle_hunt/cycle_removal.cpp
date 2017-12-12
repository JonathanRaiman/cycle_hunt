#include "Cycle.h"

#include <iostream>
#include <sys/stat.h>

using cycle_hunt::Graph;

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
            std::cout << "loaded graph with " << G.size() << " nodes" << std::endl;
            auto nroots = G.get_roots().size();
            if (nroots != 1) {
                std::cout << "Error: graph should only have 1 root, but found "
                          << nroots << "." << std::endl;
                return 1;
            }
            // iterate through nodes, ensuring their
            // predecessors were visited beforehand. Stop
            // when no new node can be visited due to
            // inactive predecessors. Record those blocking
            // nodes:
            auto bnodes = G.reach_blocking_nodes();
            std::cout << bnodes.size()
                      << " blocking nodes (circular dependencies) found"
                      << " when building DAG." << std::endl;
            // Construct a graph made from those blocking nodes:
            // that's where the cycles are hiding
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

            std::cout << "after cycle removal, graph has "
                      << G.number_of_edges() << " edges." << std::endl;

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
