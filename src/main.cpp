// Command-line driver: maximum flexi-clique with FPA or NPA (heuristics) or
// EBA (exact).
//
//   flexi --algo <fpa|npa|eba> --file <edge list> --tau <t in (0,1)>
//         [--seed <node id>] [--init <id,id,...>] [--out <file>]
//
// The variant binaries of the experiments (CMakeLists.txt) share this driver
// and print their name on the line "Variant".
//
// --seed (EBA only): the maximum flexi-clique that contains the given node.
// --init (EBA only): a known flexi-clique used as the initial incumbent when it
//                    is valid and larger than the FPA solution.
// Node ids in --seed, --init and the output are the ids of the input file.
#include "graph.h"
#include "fpa.h"
#include "eba.h"
#include "npa.h"
#include <climits>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

// Parse a whole token as a number; false on trailing characters or out-of-range values.
static bool parseDouble(const char* tok, double& x) {
    try { size_t used = 0; x = std::stod(tok, &used); return used == std::strlen(tok); }
    catch (const std::exception&) { return false; }
}
static bool parseInt(const std::string& tok, int& x) {
    try {
        size_t used = 0;
        long long v = std::stoll(tok, &used);
        if (used != tok.size() || v < INT_MIN || v > INT_MAX) return false;
        x = static_cast<int>(v);
        return true;
    } catch (const std::exception&) { return false; }
}

static void printUsage(const char* prog) {
    std::cerr << "Usage: " << prog << " --algo <fpa|npa|eba> --file <edge list> --tau <t in (0,1)>"
              << " [--seed <node id>] [--init <id,id,...>] [--out <file>]\n";
}

int main(int argc, char** argv) {
    std::string algo, file_path, out_file, init_list;
    double tau = -1.0;
    bool has_seed = false, has_init = false;
    int seed_orig = 0;

    for (int i = 1; i < argc; ++i) {
        const bool has_val = i + 1 < argc;
        if (std::strcmp(argv[i], "--algo") == 0 && has_val)      algo = argv[++i];
        else if (std::strcmp(argv[i], "--file") == 0 && has_val) file_path = argv[++i];
        else if (std::strcmp(argv[i], "--tau") == 0 && has_val) {
            if (!parseDouble(argv[++i], tau)) { std::cerr << "--tau: invalid number\n"; return 1; }
        }
        else if (std::strcmp(argv[i], "--out") == 0 && has_val)  out_file = argv[++i];
        else if (std::strcmp(argv[i], "--seed") == 0 && has_val) {
            if (!parseInt(argv[++i], seed_orig)) { std::cerr << "--seed: invalid node id\n"; return 1; }
            has_seed = true;
        }
        else if (std::strcmp(argv[i], "--init") == 0 && has_val) { init_list = argv[++i]; has_init = true; }
        else { std::cerr << "Unknown or incomplete option: " << argv[i] << "\n"; printUsage(argv[0]); return 1; }
    }
    if ((algo != "fpa" && algo != "npa" && algo != "eba") || file_path.empty()
        || !(tau > 0.0 && tau < 1.0)) {
        printUsage(argv[0]);
        return 1;
    }
    if ((has_seed || has_init) && algo != "eba") {
        std::cerr << "--seed and --init are only supported with --algo eba\n";
        return 1;
    }
#ifdef FLAT_SYMSE
    if (has_seed) {
        std::cerr << "--seed is not supported by the flat Sym-SE variant\n";
        return 1;
    }
#endif

    Graph G;
    std::unordered_map<int,int> id_map;              // original id -> internal index
    try {
        G = Graph::loadFromFile(file_path, &id_map);
    } catch (const std::exception& e) {
        std::cerr << "Error loading graph: " << e.what() << "\n";
        return 1;
    }
    auto lookup = [&](const std::string& tok, int& idx) {
        int raw;
        if (!parseInt(tok, raw)) return false;
        auto it = id_map.find(raw);
        if (it == id_map.end()) return false;
        idx = it->second;
        return true;
    };

    AlgoResult result;
    if (algo == "fpa") {
        result = runFPA(G, tau);
    } else if (algo == "npa") {
        result = runNPA(G, tau);
    } else {
        int seed = -1;
        if (has_seed && !lookup(std::to_string(seed_orig), seed)) {
            std::cerr << "--seed " << seed_orig << " is not a node of the graph\n";
            return 1;
        }
        std::vector<int> init;
        if (has_init) {
            std::stringstream ss(init_list);
            std::string tok;
            while (std::getline(ss, tok, ',')) {
                int idx;
                if (!lookup(tok, idx)) {
                    std::cerr << "--init: '" << tok << "' is not a node of the graph\n";
                    return 1;
                }
                init.push_back(idx);
            }
            if (init.empty() || init_list.back() == ',') {
                std::cerr << "--init: empty or malformed node list\n";
                return 1;
            }
        }
        result = runEBA(G, tau, seed, init.empty() ? nullptr : &init);
    }

    const bool valid = !result.nodes.empty() && isFlexiClique(G, result.nodes, tau);
    std::vector<int> orig_of(G.n, -1);
    for (const auto& [raw, idx] : id_map) orig_of[idx] = raw;
    std::string nodes;
    for (size_t k = 0; k < result.nodes.size(); ++k)
        nodes += (k ? "," : "") + std::to_string(orig_of[result.nodes[k]]);

#ifdef FLEXI_VARIANT
    std::cout << "Variant   : " << FLEXI_VARIANT << "\n";
#endif
    std::cout << "Algorithm : " << algo << "\n"
              << "Tau       : " << tau << "\n"
              << "Size      : " << result.nodes.size() << "\n"
              << "Time(s)   : " << result.time_sec << "\n"
              << "Valid     : " << (valid ? "yes" : "no") << "\n";
    if (algo == "eba")
        std::cout << "FPA_Time(s): " << result.fpa_sec << "\n"
                  << "Branches  : " << result.branches << "\n";
    std::cout << "Nodes     : " << nodes << "\n";

    if (!out_file.empty()) {
        std::error_code ec;
        fs::path p(out_file);
        if (p.has_parent_path()) fs::create_directories(p.parent_path(), ec);
        std::ofstream fout(out_file);
#ifdef FLEXI_VARIANT
        fout << "variant=" << FLEXI_VARIANT << "\n";
#endif
        fout << "algorithm=" << algo << "\ntau=" << tau << "\nsize=" << result.nodes.size()
             << "\ntime=" << result.time_sec << "\nvalid=" << (valid ? 1 : 0)
             << "\nbranches=" << result.branches << "\nnodes=" << nodes << "\n";
        fout.close();
        if (!fout) { std::cerr << "Cannot write to: " << out_file << "\n"; return 1; }
    }
    return valid || result.nodes.empty() ? 0 : 2;
}
