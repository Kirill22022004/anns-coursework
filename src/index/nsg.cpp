#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include <limits>

float euclidean_distance(const std::vector<float>& a, const std::vector<float>& b) {
    float dist = 0;
    for (size_t i = 0; i < a.size(); i++) {
        dist += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return sqrt(dist);
}

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

std::vector<float> compute_centroid(const std::vector<std::vector<float>>& data) {
    size_t N = data.size(), D = data[0].size();
    std::vector<float> centroid(D, 0);
    for (const auto& data_vec : data) {
        for (size_t i = 0; i < D; i++) {
            centroid[i] += data_vec[i];
        }
    }
    for (auto& x : centroid) {
        x /= N;
    }
    return centroid;
}

std::vector<int> search_on_graph(const std::vector<std::vector<float>>& data, const std::vector<std::vector<int>>& G, int entry, const std::vector<float>& query, int l) {
    std::priority_queue<std::pair<float, int>> candidates;
    std::unordered_set<int> visited;
    float dist = euclidean_distance(data[entry], query);
    candidates.emplace(-dist, entry);
    visited.insert(entry);
    std::vector<std::pair<float, int>> result;
    while (!candidates.empty() && result.size() < (size_t)l) {
        auto [neg_d, cur] = candidates.top();
        candidates.pop();
        result.emplace_back(-neg_d, cur);
        for (int neighbor : G[cur]) {
            if (visited.count(neighbor)) {
                continue;
            }
            visited.insert(neighbor);
            candidates.emplace(-euclidean_distance(data[neighbor], query), neighbor);
        }
    }
    sort(result.begin(), result.end());
    std::vector<int> neighbors;
    for (auto& [_, idx] : result) {
        neighbors.push_back(idx);
    }
    return neighbors;
}

bool has_conflict(const std::vector<std::vector<float>>& data, int v, int p, const std::vector<int>& R) {
    const float angle_threshold = 0.97;
    for (int r : R) {
        float dp = 0;
        for (size_t i = 0; i < data[v].size(); i++) {
            dp += (data[p][i] - data[v][i]) * data[r][i] - data[v][i];
        }
        if (dp / (euclidean_distance(data[p], data[v]) * euclidean_distance(data[r], data[v]) + 1e-9) > angle_threshold) {
            return true;
        }
    }
    return false;
}

std::vector<std::vector<int>> build_nsg(const std::vector<std::vector<float>>& data, const std::vector<std::vector<int>>& knn_graph, int l = 40, int m = 20) {
    int N = data.size();
    std::vector<std::vector<int>> nsg(N);
    std::vector<float> centroid = compute_centroid(data);
    int random_entry = rand() % N;
    int n = search_on_graph(data, knn_graph, random_entry, centroid, l)[0];
    for (int v = 0; v < N; v++) {
        std::vector<int> E = search_on_graph(data, knn_graph, n, data[v], l);
        E.insert(E.end(), knn_graph[v].begin(), knn_graph[v].end());
        std::sort(E.begin(), E.end());
        E.erase(std::unique(E.begin(), E.end()), E.end());
        std::sort(E.begin(), E.end(), [&](int a, int b) {
            return euclidean_distance(data[v], data[a]) < euclidean_distance(data[v], data[b]);
        });
        std::vector<int> R;
        R.push_back(E[0]);
        for (size_t i = 1; i < E.size() && R.size() < (size_t)m; i++) {
            int p = E[i];
            if (!has_conflict(data, v, p, R)) {
                R.push_back(p);
            }
        }
        nsg[v] = R;
    }
    while (1) {
        std::vector<bool> visited(N, false);
        std::vector<int> stack;
        stack.push_back(n);
        visited[n] = true;
        while (!stack.empty()) {
            int v = stack.back();
            stack.pop_back();
            for (int u : nsg[v]) {
                if (!visited[u]) {
                    visited[u] = true;
                    stack.push_back(u);
                }
            }
        }
        int i = 0;
        while (i < N && visited[i]) {
            i++;
        }
        if (i == N) {
            break;
        }
        float best_dist = std::numeric_limits<float>::max();
        int best_j = -1;
        for (int j = 0; j < N; j++) {
            if (visited[j]) {
                float d = euclidean_distance(data[i], data[j]);
                if (d < best_dist) {
                    best_dist = d;
                    best_j = j;
                }
            }
        }
        nsg[i].push_back(best_j);
        nsg[best_j].push_back(i);
    }
    return nsg;
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
    std::vector<std::vector<int>> knn_graph(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        std::vector<std::pair<float, int>> dists;
        for (size_t j = 0; j < data.size(); ++j) {
            if (i == j) {
                continue;
            }
            dists.emplace_back(euclidean_distance(data[i], data[j]), j);
        }
        sort(dists.begin(), dists.begin() + 40);
        for (int k = 0; k < 40; ++k) {
            knn_graph[i].push_back(dists[k].second);
        }
    }
    std::vector<std::vector<int>> nsg = build_nsg(data, knn_graph);
    save_graph(nsg, "nsg.txt");
    std::cout << "NSG built and saved." << std::endl;
    return 0;
}
