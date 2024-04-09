/*#include "m_matrix.hpp"
#include "parser.hpp"

//Размерность задачи
int dimension = 2;

//Аббревиатуры для компонентов x, y
int comp_x = 0;
int comp_y = 1;

//Координаты узлов элементов
std::vector<std::vector<double>> nodes_coord = get_nodes_coord("C:\\Projects\\Calc_elastic_2D\\src\\data\\mesh_for_calc.inc");

//Формирование элементов из узлов
std::vector<std::vector<int>> elements = get_elements("C:\\Projects\\Calc_elastic_2D\\src\\data\\mesh_for_calc.inc");

//Аббревиатуры для компонентов xi, eta
int comp_xi = 0;
int com_eta = 1;    

//Количество элементов в сеточной модели
int elements_count = 2;

//Количество узлов в сеточной модели
int nodes_count = 6;

//Количество узлов в элементе
int nodes_per_element = 4;

//Количество степеней свободы (DOF) в узле:
int DOFs_per_node = dimension;

//Количество независимых компонентов тензора деформаций
int tensor_comp_count = (pow(dimension, 2) + dimension)/2;

//Толщина пластины
double delta = 0.001;

//Интерполяционные функции для четырех узлов
double phi(int node, double xi, double eta) {
    double c = 0;
    switch (node)
    {
    case 0:
        c = 1.0/4 * (1 - xi) * (1 - eta); //-1,-1
        break;
    case 1:
        c = 1.0/4 * (1 + xi) * (1 - eta); //1, -1
        break;
    case 2:
        c = 1.0/4 * (1 + xi) * (1 + eta); //1, 1
        break;
    case 3:
        c = 1.0/4 * (1 - xi) * (1 + eta); //-1, 1    
        break;
    default:
        break;
    }
    return c;
}

//Производные интерполяционных функций для четырех узлов
std::vector<double> dPhidEta(int node, double xi, double eta) {
    double c0 = 0;
    double c1 = 0;
    switch (node)
    {
    case 0:
        c0 = 1.0/4 * (eta - 1);
        c1 = 1.0/4 * (xi - 1); 
        break;
    case 1:
        c0 = -1 * 1.0/4 * (eta - 1);
        c1 = -1 * 1.0/4 * (xi + 1); 
        break;
    case 2:
        c0 = 1.0/4 * (eta + 1); 
        c1 = 1.0/4 * (xi + 1);
        break;
    case 3:
        c0 = -1 * 1.0/4 * (eta + 1);   
        c1 = -1 * 1.0/4 * (xi - 1);
        break;
    default:
        break;
    }
    std::vector<double> c = {c0, c1};
    return c;
}

//Количество точек интегрирования в элементе
int IP_count = 16;

//Количество точек интегрирования в элементе по каждому направлению 
int IP_per_direction = pow(IP_count, 1.0/dimension);

//Координаты точек интегрирования в нормализованной системе координат элемента
std::vector<double> IP_coord = {-0.861136312, -0.339981046, 0.339981046, 0.861136312};

//Весовые коэффициенты для точек интегрирования
std::vector<double> IP_h = {0.3478548451, 0.6521451548, 0.6521451548, 0.3478548451};

//Вычисление матрицы Якоби
std::vector<std::vector<double>> calc_J(double xi, double eta, int elem) {
    std::vector<std::vector<double>> J(dimension, std::vector<double>(dimension, 0.0));

    for(int node = 0; node < nodes_per_element; ++node) {
        std::vector<double> dphi = dPhidEta(node, xi, eta);
        for(int i = 0; i < dimension; ++i) {
            for(int j = 0; j < dimension; ++j) {
                J[i][j] = J[i][j] + dphi[i] * nodes_coord[elements[elem][node]][j];
            }
        }
    }

    return J;
}

//Вычисление производных базисных функций по dx и dy 
std::vector<double> dPhidX(int node, double xi, double eta, int elem) {
    std::vector<std::vector<double>> J = calc_J(xi, eta, elem);
    std::vector<std::vector<double>> J_inv = matrixInverse(J);
    std::vector<double> dphi = dPhidEta(node, xi, eta);
    std::vector<double> H = matrixMultiplicationVector(J_inv, dphi);
    return H;
}

//Вычисление матицы частных производных базисных функций
std::vector<std::vector<double>> calc_B(double xi, double eta, int elem) {
    std::vector<std::vector<double>> B(tensor_comp_count, std::vector<double>(nodes_per_element*DOFs_per_node, 0.0));
    int i1 = 0;
    int i2 = 0;

    for(int i = 0; i < nodes_per_element; ++i) {
        std::vector<double> dphi = dPhidX(i, xi, eta, elem);
        i1 = 2*i;
        B[0][i1] = dphi[comp_x];
        B[1][i1] = 0;
        B[2][i1] = dphi[comp_y];

        i2 = i1 + 1;
        B[0][i2] = 0;
        B[1][i2] = dphi[comp_y];
        B[2][i2] = dphi[comp_x];
    }

    return B;
}

//Модуль упругости
double E = 2e+11;

//Коэффициент Пуассона
double nu = 0.3;

//Вычисление матрицы упругих свойств [D]
//Плоское напряженное состояние
double alpha = E/(1 - pow(nu, 2));

std::vector<std::vector<double>> D_1 = {{alpha, alpha * nu, 0},
                                        {alpha * nu, alpha, 0},
                                        {0, 0, alpha * (1 - nu)/2}};

//Плоское деформированное состояние
double beta = E * (1 - nu)/((1 + nu) * (1 - 2*nu));

std::vector<std::vector<double>> D_2 = {{beta, beta * nu/(1 - nu), 0},
                                        {beta * nu/(1 - nu), beta, 0},
                                        {0, 0, beta*(1 - 2 * nu)/(2*(1 - nu))}};

//Выбор постановки
std::vector<std::vector<double>> D = D_1;

//Вычисление матрицы жесткости конечного элемента [Ke]
std::vector<std::vector<double>> calc_Ke(int elem) {
    std::vector<std::vector<double>> K(nodes_per_element * DOFs_per_node,
                                       std::vector<double>(nodes_per_element * DOFs_per_node, 0.0));

    for(int IP_xi = 0; IP_xi < IP_per_direction; ++IP_xi) {
        for(int IP_eta = 0; IP_eta < IP_per_direction; ++IP_eta) {
            std::vector<std::vector<double>> B = calc_B(IP_coord[IP_xi], 
                                                        IP_coord[IP_eta],
                                                        elem);   

            std::vector<std::vector<double>> K_IP = matrixMultiplicationDouble(matrixMultiplication(matrixMultiplication(matrixTransposition(B), D), B), delta);

            for(int i = 0; i < nodes_per_element * DOFs_per_node; ++i) {
                for(int j = 0; j < nodes_per_element * DOFs_per_node; ++j) {
                    K[i][j] = K[i][j] + K_IP[i][j] * (IP_h[IP_xi] * IP_h[IP_eta]) * matrixDeterminant(calc_J(IP_coord[IP_xi], IP_coord[IP_eta], elem));
                }
            }     
        }
    }

    return K;
}

//Вычисление матрицы жесткости системы [K]:
std::vector<std::vector<double>> calc_K(int param) {
    std::vector<std::vector<double>> K(nodes_count * DOFs_per_node,
                                       std::vector<double>(nodes_count * DOFs_per_node, 0.0));

    for(int elem = 0; elem < elements_count; ++elem) {
        std::vector<std::vector<double>> Ke = calc_Ke(elem);
        for(int i_node = 0; i_node < nodes_per_element; ++i_node) {
            for(int j_node = 0; j_node < nodes_per_element; ++j_node) {
                for(int i_dof = 0; i_dof < DOFs_per_node; ++i_dof) {
                    for(int j_dof = 0; j_dof < DOFs_per_node; ++j_dof) {
                        int ig = (elements[elem][i_node]) * DOFs_per_node + i_dof;
                        int jg = (elements[elem][j_node]) * DOFs_per_node + j_dof;
                        int ie = (i_node * DOFs_per_node) + i_dof;
                        int je = (j_node * DOFs_per_node) + j_dof;

                        K[ig][jg] = K[ig][jg] + Ke[ie][je];
                    }
                }
            }
        }
    }

    return K;
}

//Вектор внешних сил {F}
std::vector<std::vector<double>> create_F(int param) {
    return std::vector<std::vector<double>>(nodes_count*DOFs_per_node,std::vector<double>(1, 0.0));
}

std::vector<std::vector<double>> F = create_F(0);

//Граничные условия(ГУ)
void Add_BC_U(double value, 
              int dof, 
              std::vector<std::vector<double>> NS,
              std::vector<std::vector<double>> &K,
              std::vector<std::vector<double>> &F) {
    for(int i = 0; i < NS.size(); ++i) {
        int ii = (NS[i][0] - 1) * DOFs_per_node + dof;
        for(int j = 0; j < nodes_count*DOFs_per_node; ++j) {
            F[j][0] = F[j][0] - K[j][ii] * value;
        }
    }

    for(int i = 0; i < NS.size(); ++i) {
        int ii = (NS[i][0] - 1) * DOFs_per_node + dof;
        for(int j = 0; j < nodes_count*DOFs_per_node; ++j) {
            K[ii][j] = 0;
            K[j][ii] = 0;
        }
        K[ii][ii] = 1;
        F[ii][0] = value;
    }
}

void Add_BC_F(double value,
              int dof, 
              std::vector<std::vector<double>> NS,
              std::vector<std::vector<double>> &F) {
    for(int i = 0; i < NS.size(); ++i) {
        int ii = (NS[i][0] - 1) * DOFs_per_node + dof;
        F[ii][0] = value;
    }
}

//Именованные выборки узлов (NS):
std::vector<std::vector<double>> NS_1 = {{1}, {4}};

std::vector<std::vector<double>> NS_2 = {{1}, {2}, {5}};

std::vector<std::vector<double>> NS_3 = {{5}, {6}};



int main() {
    
    std::cout << "Hello world" << std::endl;
    //Вычислим матрицу жесткости
    std::vector<std::vector<double>> K = calc_K(0);

    //Добавим граничые условия
    Add_BC_F(100, comp_x, NS_3, F);
    Add_BC_U(0, comp_x, NS_1, K, F);
    Add_BC_U(0, comp_y, NS_1, K, F);

    //Матрица жесткости
    std::cout << "K = " << std::endl;
    printMatrix(K);
    //Вектор внешних сил
    std::cout << "F = " << std::endl;
    printMatrix(F);
    
    //Решение СЛАУ
    std::vector<std::vector<double>> U = gauss(K, F);
    std::cout << "U = " << std::endl;
    printMatrix(U);
    
    
    //printMatrix(get_elements("data\\example_mesh.inc"));
    return 0;
}
*/

/*
#include <GLFW/glfw3.h>

int main(void)
{
    GLFWwindow* window;

    // Initialize the library 
    if (!glfwInit())
        return -1;

    // Create a windowed mode window and its OpenGL context 
    window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current 
    glfwMakeContextCurrent(window);

    // Loop until the user closes the window 
    while (!glfwWindowShouldClose(window))
    {
        // Render here 
        //glClear(GL_COLOR_BUFFER_BIT);

        // Swap front and back buffers 
        glfwSwapBuffers(window);

        // Poll for and process events 
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
*/