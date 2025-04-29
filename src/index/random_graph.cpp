#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <queue>
#include <limits>
#include <numeric>
#include <random>

std::vector<std::vector<float>> load_vectors(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<std::vector<float>> data;
    while (getline(file, line)) {
        std::istringstream iss(line);
        std::vector<float> data_vec;
        float x;
        while (iss >> x) {
            data_vec.push_back(x);
        }
        data.push_back(data_vec);
    }
    std::cout << data.size() << " " << data[0].size() << std::endl;
    return data;
}

void save_graph(const std::vector<std::vector<int>>& G, const std::string& filename) {
    std::ofstream out(filename);
    for (const auto& neighbors : G) {
        for (size_t i = 0; i < neighbors.size(); ++i) {
            out << neighbors[i];
            if (i + 1 < neighbors.size()) {
                out << ' ';
            }
        }
        out << '\n';
    }
}

int main() {
    auto data = load_vectors("../../datasets/train.txt");
    std::vector<std::vector<int>> graph(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        std::vector<int> cand(data.size());
        std::iota(cand.begin(), cand.end(), 0);
        std::shuffle(cand.begin(), cand.end(), std::mt19937{std::random_device{}()});
        graph[i] = std::vector<int>(cand.begin(), cand.begin() + 20);
    }
    save_graph(graph, "random.txt");
    std::cout << "Random graph built and saved." << std::endl;
    return 0;
}
