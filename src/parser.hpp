/*Парсирование данных с файла осуществляется по 
ключеввым словам *NODE и *ELEMENT_SHELL. Описаны две функции: 
get_nodes_coord -- которая возвращает матрицу содержащую координаты узлов
get_elements -- которая возвращает матрицу содержащую номера узлов из которых состоит элемент*/
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>

//Получить список координат узлов из файла (*NODE - ключевое слово)
std::vector<std::vector<double>> get_nodes_coord(std::string file_name) {
    std::vector<std::vector<double>> nodes_coord; //матрица координат узлов
    std::ifstream in(file_name); //открываем файл для чтения
    std::string line; //будем считывать информацию с файла построчно
    if(in.is_open()) {
        while(std::getline(in, line)) {
            //Если встретили ключевое слово *NODE
            if(line.find("*NODE") != std::string::npos) { 
                while(std::getline(in, line)) {    
                    if(line.find("*") != std::string::npos) {
                        return nodes_coord;
                    }
                    if(line.find("$") == std::string::npos) {
                        std::istringstream ss(line); //добавляем в поток ввода-вывода строку line
                        int NID;
                        double X, Y, Z;
                        ss >> NID >> X >> Y >> Z; //считываем данные из потока ss
                        nodes_coord.push_back({X, Y}); //положим координаты в матрицу
                    }
                }
            }
        }
    }
    in.close(); //закроем файл
    return nodes_coord;
}

//Получить список элементов, состоящих из номеров узлов
std::vector<std::vector<int>> get_elements(std::string file_name) {
    std::vector<std::vector<int>> elements;
    std::ifstream in(file_name);
    std::string line;
    if(in.is_open()) {
        while(std::getline(in, line)) {
            if(line.find("*ELEMENT_SHELL") != std::string::npos) {
                while(std::getline(in, line)) {
                    if(line.find("*") != std::string::npos) {
                        return elements;
                    }
                    if(line.find("$") == std::string::npos) {
                        std::istringstream ss(line);
                        int EID;
                        int PID;
                        int N1;
                        int N2;
                        int N3;
                        int N4;
                        ss >> EID >> PID >> N1 >> N2 >> N3 >> N4;
                        elements.push_back({N1, N2, N3, N4});
                    }
                }
            }
        }
    }
    in.close();
    return elements;
}