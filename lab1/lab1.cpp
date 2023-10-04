#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

bool Comparator(const std::pair<int, std::string>& a, const std::pair<int, std::string>& b) {
    if (a.first == b.first) {
        return a.second < b.second;
    }

    return a.first > b.first;
}

int get_som(const std::string& name, int weight) {
    int som = name.length();

    for (char c : name) {
        if (c >= 'A' && c <= 'Z') {
            som += c - 'A' + 1;
        }

        if (c >= 'a' && c <= 'z') {
            som += c - 'a' + 1;
        }
    }

    return som * weight;
}

int main() {
    std::string names[6] = {"COLIN", "AMANDBA", "AMANDAB", "CAROL", "PauL", "JOSEPH"};
    int weights[6] = {1, 4, 4, 5, 2, 1};

    int n = 3;

    std::vector<std::pair<int, std::string>> somVector;

    for (int i = 0; i < 6; ++i) {
        int som = get_som(names[i], weights[i]);
        somVector.push_back(std::make_pair(som, names[i]));
    }

    std::sort(somVector.begin(), somVector.end(), Comparator);

    for (const auto& entry : somVector) {
        std::cout << "SOM for " << entry.second << ": " << entry.first << "\n";
    }

    std::cout << "\n\n" + somVector.at(n - 1).second << "\n";

    return 0;
}

