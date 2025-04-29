#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <unordered_set>
#include <queue>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <set>

using namespace std::chrono;

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
    return data;
}

std::vector<std::vector<int>> load_graph(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<std::vector<int>> G;
    while (getline(file, line)) {
        std::istringstream iss(line);
        std::vector<int> neighbors;
        int id;
        while (iss >> id) {
            neighbors.push_back(id);
        }
        G.push_back(neighbors);
    }
    return G;
}

std::vector<std::vector<int>> load_answers(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<std::vector<int>> answers;
    while (getline(file, line)) {
        std::istringstream iss(line);
        std::vector<int> row;
        int id;
        while (iss >> id) {
            row.push_back(id);
        }
        answers.push_back(row);
    }
    return answers;
}

std::vector<int> search_bruteforce(const std::vector<std::vector<float>>& data, const std::vector<std::vector<int>>& G, const std::vector<float>& query, int K, int L) {
    std::vector<std::pair<float, int>> res(data.size());
    for (int i = 0; i < data.size(); i++) {
        res[i] = {euclidean_distance(data[i], query), i};
    }
    std::sort(res.begin(), res.end());
    std::vector<int> result;
    for (int i = 0; i < K; i++) {
        result.push_back(res[i].second);
    }
    return result;
}

std::vector<int> search_on_graph(const std::vector<std::vector<float>>& data, const std::vector<std::vector<int>>& G, const std::vector<float>& query, int K, int L) {
    std::vector<std::pair<int, float>> S = {{0, euclidean_distance(data[0], query)}};
    std::unordered_set<int> checked, used = {0};
    while (true) {
        int i = 0;
        while (i < int(S.size()) && checked.find(S[i].first) != checked.end()) {
            i++;
        }
        if (i >= int(S.size())) {
            break;
        }
        checked.insert(S[i].first);
        for (auto& u : G[S[i].first]) {
            if (used.find(u) == used.end()) {
                used.insert(u);
                S.push_back({u, euclidean_distance(data[u], query)});
            }
        }
        std::sort(S.begin(), S.end(), [&](const auto& a, const auto& b) {
            return a.second < b.second;
        });
        if (int(S.size()) > L) {
            S.resize(L);
        }
    }
    std::vector<int> result(std::min(K, (int)S.size()));
    for (size_t i = 0; i < std::min(K, (int)S.size()); i++) {
        result[i] = S[i].first;
    }
    return result;
}

std::vector<int> search_on_graph_with_set(const std::vector<std::vector<float>>& data, const std::vector<std::vector<int>>& G, const std::vector<float>& query, int K, int L) {
    std::set<std::pair<float, int>> S = {{euclidean_distance(data[0], query), 0}};
    std::vector<int> checked(data.size());
    std::vector<int> used(data.size());
    used[0] = 1;
    while (true) {
        auto it = S.begin();
        while (it != S.end() && checked[it->second]) {
            it++;
        }
        if (it == S.end()) {
            break;
        }
        int v = it->second;
        checked[v] = 1;
        for (auto& u : G[v]) {
            if (!used[u]) {
                used[u] = 1;
                S.insert({euclidean_distance(data[u], query), u});
            }
        }
        while (int(S.size()) > L) {
            S.erase(--S.end());
        }
    }
    std::vector<int> result(std::min(K, (int)S.size()));
    auto it = S.begin();
    for (size_t i = 0; i < std::min(K, (int)S.size()); i++) {
        result[i] = it->second;
        it++;
    }
    return result;
}

struct SegmentTree {
    struct Node {
        int v;
        float dv;
        bool checked;
        float max_dv;
        std::pair<float, int> to_check;

        Node() : v(-1), dv(std::numeric_limits<float>::max()), checked(true), max_dv(std::numeric_limits<float>::max()),
        to_check({-1, -1}) {}
    };

    std::vector<Node> Nodes;

    void build(int L) {
        Nodes.resize(4 * L);
    }

    void pull(int v) {
        Nodes[v].max_dv = std::max(Nodes[v * 2 + 1].max_dv, Nodes[v * 2 + 2].max_dv);
        if (Nodes[v * 2 + 1].to_check.second == -1) {
            Nodes[v].to_check = Nodes[v * 2 + 2].to_check;
        } else if (Nodes[v * 2 + 2].to_check.second == -1) {
            Nodes[v].to_check = Nodes[v * 2 + 1].to_check;
        } else {
            Nodes[v].to_check = std::min(Nodes[v * 2 + 1].to_check, Nodes[v * 2 + 2].to_check);
        }
    }

    void update(int v, int l, int r, std::pair<float, int> info) {
        if (Nodes[v].max_dv <= info.first) {
            return;
        }
        if (l + 1 == r) {
            Nodes[v].v = info.second;
            Nodes[v].dv = info.first;
            Nodes[v].checked = false;
            Nodes[v].max_dv = info.first;
            Nodes[v].to_check = info;
        } else {
            int m = (l + r) / 2;
            if (Nodes[v * 2 + 1].max_dv > Nodes[v * 2 + 2].max_dv) {
                update(v * 2 + 1, l, m, info);
            } else {
                update(v * 2 + 2, m, r, info);
            }
            pull(v);
        }
    }

    void checked(int v, int l, int r, int info) {
        if (l + 1 == r) {
            Nodes[v].checked = true;
            Nodes[v].to_check = {-1, -1};
            return;
        }
        int m = (l + r) / 2;
        if (Nodes[v * 2 + 1].to_check.second == info) {
            checked(v * 2 + 1, l, m, info);
        } else {
            checked(v * 2 + 2, m, r, info);
        }
        pull(v);
    }

    void dfs(int v, int l, int r, std::vector<std::pair<float, int>>& res) {
        if (l + 1 == r) {
            if (Nodes[v].v != -1) {
                res.push_back({Nodes[v].dv, Nodes[v].v});
            }
            return;
        }
        int m = (l + r) / 2;
        dfs(v * 2 + 1, l, m, res);
        dfs(v * 2 + 2, m, r, res);
    }

    std::vector<int> get(const int K) {
        std::vector<std::pair<float, int>> res;
        dfs(0, 0, Nodes.size() / 4, res);
        std::sort(res.begin(), res.end());
        std::vector<int> result;
        for (size_t i = 0; i < std::min(res.size(), size_t(K)); i++) {
            result.push_back(res[i].second);
        }
        return result;
    }
};

std::vector<int> search_on_graph_with_segment_tree(const std::vector<std::vector<float>>& data, const std::vector<std::vector<int>>& G, const std::vector<float>& query, int K, int L) {
    SegmentTree S;
    S.build(L);
    S.update(0, 0, L, {euclidean_distance(data[0], query), 0});
    std::vector<int> used(data.size());
    used[0] = 1;
    while (true) {
        int v = S.Nodes[0].to_check.second;
        if (v == -1) {
            break;
        }
        S.checked(0, 0, L, v);
        for (auto& u : G[v]) {
            if (!used[u]) {
                used[u] = 1;
                float du = euclidean_distance(data[u], query);
                S.update(0, 0, L, {du, u});
            }
        }
    }
    return S.get(K);
}

std::vector<int> search_on_graph_with_pool(const std::vector<std::vector<float>>& data, const std::vector<std::vector<int>>& G, const std::vector<float>& query, int K, int L) {
    std::vector<std::pair<float, int>> S = {{euclidean_distance(data[0], query), 0}};
    std::vector<int> checked(data.size());
    std::vector<int> used(data.size());
    used[0] = 1;
    while (true) {
        size_t it = 0;
        while (it != S.size() && checked[S[it].second]) {
            it++;
        }
        if (it == S.size()) {
            break;
        }
        int v = S[it].second;
        checked[v] = 1;
        for (auto& u : G[v]) {
            if (!used[u]) {
                used[u] = 1;
                float du = euclidean_distance(data[u], query);
                if (int(S.size()) == L && S.back().first <= du) {
                    continue;
                }
                int j = int(S.size());
                S.push_back({du, u});
                while (j > 0 && S[j].first < S[j - 1].first) {
                    std::swap(S[j], S[j - 1]);
                    j--;
                }
                if (int(S.size()) == L + 1) {
                    S.pop_back();
                }
            }
        }
    }
    std::vector<int> result(std::min(K, (int)S.size()));
    auto it = S.begin();
    for (size_t i = 0; i < std::min(K, (int)S.size()); i++) {
        result[i] = it->second;
        it++;
    }
    return result;
}

std::vector<int> search_on_graph_with_pool_stopping(const std::vector<std::vector<float>>& data, const std::vector<std::vector<int>>& G, const std::vector<float>& query, int K, int L, int Stopping = 0) {
    if (!Stopping) {
        Stopping = L;
    }
    std::vector<std::pair<float, int>> S = {{euclidean_distance(data[0], query), 0}};
    std::vector<int> checked(data.size());
    std::vector<int> used(data.size());
    used[0] = 1;
    while (true) {
        size_t it = 0;
        while (it != S.size() && checked[S[it].second]) {
            it++;
        }
        if (it == S.size() || it >= Stopping) {
            break;
        }
        int v = S[it].second;
        checked[v] = 1;
        for (auto& u : G[v]) {
            if (!used[u]) {
                used[u] = 1;
                float du = euclidean_distance(data[u], query);
                if (int(S.size()) == L && S.back().first <= du) {
                    continue;
                }
                int j = int(S.size());
                S.push_back({du, u});
                while (j > 0 && S[j].first < S[j - 1].first) {
                    std::swap(S[j], S[j - 1]);
                    j--;
                }
                if (int(S.size()) == L + 1) {
                    S.pop_back();
                }
            }
        }
    }
    std::vector<int> result(std::min(K, (int)S.size()));
    auto it = S.begin();
    for (size_t i = 0; i < std::min(K, (int)S.size()); i++) {
        result[i] = it->second;
        it++;
    }
    return result;
}

float compute_recall(const std::vector<int>& result, const std::vector<int>& answers, int K) {
    std::unordered_set<int> truth(answers.begin(), answers.begin() + K);
    int cnt = 0;
    for (int id : result) {
        cnt += int(truth.count(id));
    }
    return static_cast<float>(cnt) / K;
}

auto calculate(const std::string& index_filename, const std::string& method, int K, int L, int Stopping = 0) {
    auto base = load_vectors("../../datasets/train.txt");
    auto queries = load_vectors("../../datasets/queries.txt");
    auto G = load_graph(index_filename);
    auto answers = load_answers("../../datasets/answers.txt");
    double total_recall = 0.0;
    auto start_time = high_resolution_clock::now();
    for (size_t i = 0; i < queries.size(); i++) {
        std::vector<int> result;
        if (method == "bruteforce") {
            result = search_bruteforce(base, G, queries[i], K, L);
        }
        if (method == "search on graph") {
            result = search_on_graph(base, G, queries[i], K, L);
        }
        if (method == "search on graph with set") {
            result = search_on_graph_with_set(base, G, queries[i], K, L);
        }
        if (method == "search on graph with segment tree") {
            result = search_on_graph_with_segment_tree(base, G, queries[i], K, L);
        }
        if (method == "search on graph with pool") {
            result = search_on_graph_with_pool(base, G, queries[i], K, L);
        }
        if (method == "search on graph with pool and stopping") {
            result = search_on_graph_with_pool_stopping(base, G, queries[i], K, L, Stopping);
        }
        float recall = compute_recall(result, answers[i], K);
        total_recall += recall;
    }
    auto end_time = high_resolution_clock::now();
    double elapsed = duration_cast<milliseconds>(end_time - start_time).count();
    elapsed = std::max(elapsed, 1e-9);
    double qps = 1000.0 * queries.size() / elapsed;
    return std::pair<float, float>{total_recall / queries.size(), qps};
}

int main() {
    std::cout << std::fixed << std::setprecision(3);
    std::vector<std::string> paths = {"../../src/index/nssg.txt", "../../src/index/nsg.txt", "../../src/index/random.txt"};
    for (const std::string& index_filename : paths) {
        for (const int K : std::vector<int>{1, 10}) {
            std::cout << index_filename << " " << K << std::endl;
            std::string method = "search on graph";
            std::string name = "random";
            if (index_filename == paths[0]) {
                name = "nssg";
            } else if (index_filename == paths[1]) {
                name = "nsg";
            }
            std::ofstream file("../results/" + name + " + K:" + std::to_string(K) + ", article baseline.txt");
            for (int L = 1; L <= 2000; L += 100) {
                auto [recall, qps] = calculate(index_filename, method, K, L);
                file << "L:" << L << ",recall:" << recall << ",qps:" << qps << std::endl;
            }
            file.close();
        }
    }
    return 0;
//    std::vector<std::string> paths = {"../../src/index/nssg.txt", "../../src/index/nsg.txt", "../../src/index/random.txt"};
//    for (const std::string& index_filename : paths) {
//        for (const int K : std::vector<int>{1, 10}) {
//            std::cout << index_filename << " " << K << std::endl;
//            std::string method = "search on graph with set";
//            std::string name = "random";
//            if (index_filename == paths[0]) {
//                name = "nssg";
//            } else if (index_filename == paths[1]) {
//                name = "nsg";
//            }
//            std::ofstream file("../results/" + name + " + K:" + std::to_string(K) + ", set.txt");
//            for (int L = 1; L <= 2000; L += 100) {
//                auto [recall, qps] = calculate(index_filename, method, K, L);
//                file << "L:" << L << ",recall:" << recall << ",qps:" << qps << std::endl;
//            }
//            file.close();
//        }
//    }
//    return 0;
//    std::vector<std::string> paths = {"../../src/index/nssg.txt", "../../src/index/nsg.txt", "../../src/index/random.txt"};
//    for (const std::string& index_filename : paths) {
//        for (const int K : std::vector<int>{1, 10}) {
//            std::cout << index_filename << " " << K << std::endl;
//            std::string method = "search on graph with pool";
//            std::string name = "random";
//            if (index_filename == paths[0]) {
//                name = "nssg";
//            } else if (index_filename == paths[1]) {
//                name = "nsg";
//            }
//            std::ofstream file("../results/" + name + " + K:" + std::to_string(K) + ", pool.txt");
//            for (int L = 1; L <= 2000; L += 100) {
//                auto [recall, qps] = calculate(index_filename, method, K, L);
//                file << "L:" << L << ",recall:" << recall << ",qps:" << qps << std::endl;
//            }
//            file.close();
//        }
//    }
//    return 0;
//    std::vector<std::string> paths = {"../../src/index/nssg.txt", "../../src/index/nsg.txt", "../../src/index/random.txt"};
//    for (const std::string& index_filename : paths) {
//        for (const int K : std::vector<int>{1, 10}) {
//            std::cout << index_filename << " " << K << std::endl;
//            std::string method = "search on graph with pool and stopping";
//            std::string name = "random";
//            if (index_filename == paths[0]) {
//                name = "nssg";
//            } else if (index_filename == paths[1]) {
//                name = "nsg";
//            }
//            std::ofstream file("../results/" + name + " + K:" + std::to_string(K) + ", pool + stopping.txt");
//            for (int L = 1; L <= 2000; L += 100) {
//                auto [recall, qps] = calculate(index_filename, method, K, L, 5000);
//                file << "L:" << L << ",recall:" << recall << ",qps:" << qps << std::endl;
//            }
//            file.close();
//        }
//    }
//    return 0;
//    std::vector<std::string> paths = {"../../src/index/nssg.txt", "../../src/index/nsg.txt", "../../src/index/random.txt"};
//    for (const std::string& index_filename : paths) {
//        for (const int K : std::vector<int>{1, 10}) {
//            std::string method = "search on graph with segment tree";
//            std::string name = "random";
//            if (index_filename == paths[0]) {
//                name = "nssg";
//            } else if (index_filename == paths[1]) {
//                name = "nsg";
//            }
//            std::ofstream file("../results/" + name + " + K:" + std::to_string(K) + ", segtree.txt");
//            for (int L = 1; L <= 2000; L += 100) {
//                auto [recall, qps] = calculate(index_filename, method, K, L);
//                file << "L:" << L << ",recall:" << recall << ",qps:" << qps << std::endl;
//            }
//            file.close();
//        }
//    }
//    return 0;
//    {
//        std::string index_filename = "../../src/index/nssg.txt";
//        std::string method = "search on graph with segment tree";
//        int K = 10;
//        for (int L = 1; L <= 2000; L += 100) {
//            auto [recall, qps] = calculate(index_filename, method, K, L);
//            std::cout << recall << " " << qps << std::endl;
//        }
//    }
//    return 0;
//    {
//        std::string index_filename = "../../src/index/nssg.txt";
//        std::string method = "search on graph with pool and stopping";
//        int K = 10;
//        for (int L = 1; L <= 2000; L += 100) {
//            int stp = 1;
//            for (; stp <= L; stp *= 2) {
//                auto [recall, qps] = calculate(index_filename, method, K, L, stp);
//                if (recall >= 0.95) break;
//            }
//            auto [recall, qps] = calculate(index_filename, method, K, L, stp);
//            std::cout << recall << " " << qps << std::endl;
//        }
//    }
//    return 0;
//    {
//        std::string index_filename = "../../src/index/nssg.txt";
//        std::string method = "search on graph with pool";
//        int K = 10;
//        for (int L = 1; L <= 2000; L += 100) {
//            auto [recall, qps] = calculate(index_filename, method, K, L);
//            std::cout << recall << " " << qps << std::endl;
//        }
//    }
//    return 0;
//    {
//        std::string index_filename = "../../src/index/nssg.txt";
//        std::string method = "search on graph with set";
//        int K = 10;
//        for (int L = 1; L <= 2000; L += 100) {
//            auto [recall, qps] = calculate(index_filename, method, K, L);
//            std::cout << recall << " " << qps << std::endl;
//        }
//    }
//    return 0;
//    {
//        std::string index_filename = "../../src/index/nsg.txt";
//        std::string method = "search on graph";
//        int K = 10;
//        for (int L = 1; L <= 2000; L += 100) {
//            auto [recall, qps] = calculate(index_filename, method, K, L);
//            std::cout << recall << " " << qps << std::endl;
//        }
//    }
//    {
//        std::string index_filename = "../../src/index/nsg.txt";
//        auto [recall, qps] = calculate(index_filename, "bruteforce", 10, 0);
//        std::cout << recall << " " << qps << std::endl;
//    }
}