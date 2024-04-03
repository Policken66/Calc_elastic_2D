#include <fstream>
#include <string>
#include <vector>

std::vector<std::vector<double>> get_nodes_coord(std::string file_name) {
    std::vector<std::vector<double>> nodes_coord;
    std::ifstream in(file_name);
    std::string line;
    if(in.is_open()) {
        while(std::getline(in, line)) {
            size_t found = 
            if()
        }

    }
    in.close();
    return nodes_coord;
}

std::vector<std::vector<double>> get_elements(std::string file_name) {
    std::vector<std::vector<double>> elements;

    return elements;
}