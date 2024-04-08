#include <iostream>
#include <vector>
#include <cmath>

//Вывод вектора в консоль
void printVector(std::vector<double> A) {
    int n = A.size();

    for(int i = 0; i < n; ++i) {
        std::cout << A[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

//Вывод матрицы в консоль
void printMatrix(std::vector<std::vector<double>> A) {
    int n = A.size();
    int m = A[0].size();

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < m; ++j) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

//Умножение матриц
std::vector<std::vector<double>> matrixMultiplication(const std::vector<std::vector<double>>& A, 
                                                      const std::vector<std::vector<double>>& B) {
    int n = A.size();
    int m = A[0].size();
    int p = B[0].size();
    

    std::vector<std::vector<double>> result(n, std::vector<double>(p, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            for (int k = 0; k < m; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

//Умножение матрицы на скаляр
std::vector<std::vector<double>> matrixMultiplicationDouble(const std::vector<std::vector<double>>& A,
                                                            double scalar) {
    int n = A.size();
    int m = A[0].size();

    std::vector<std::vector<double>> result(n, std::vector<double>(m, 0.0));

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < m; ++j) {
            result[i][j] = A[i][j]*scalar;
        }
    }

    return result;
}

//Умножение матрицы на вектор 
std::vector<double> matrixMultiplicationVector(const std::vector<std::vector<double>> A,
                                               const std::vector<double> &x) {
    int n = A.size();
    int m = A[0].size();

    std::vector<double> result(n, 0.0);   

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < m; ++j) {
            result[i] += A[i][j] * x[j];
        }
    }    

    return result;                                        
}

//Сложение матриц
std::vector<std::vector<double>> matrixAddition(const std::vector<std::vector<double>>& A,
                                           const std::vector<std::vector<double>>& B) {
    int Na = A.size();
    int Nb = B.size();
    int Ma = A[0].size();
    int Mb = B[0].size();
    
    std::vector<std::vector<double>> result(Na, std::vector<double>(Ma, 0.0));
    //проверка на одинаковые размерности слагаемых матриц
    if(Na != Nb || Ma != Mb) { 
        std::cout << "Erorr addition" << std::endl;
        return result;
    }

    for(int i = 0; i < Na; ++i) {
        for(int j = 0; j < Ma; ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result; 
}

//Вычитание матриц
std::vector<std::vector<double>> matrixSubtraction(const std::vector<std::vector<double>>& A,
                                           const std::vector<std::vector<double>>& B) {
    int Na = A.size();
    int Nb = B.size();
    int Ma = A[0].size();
    int Mb = B[0].size();
    
    std::vector<std::vector<double>> result(Na, std::vector<double>(Ma, 0.0));
    //проверка на одинаковые размерности слагаемых матриц
    if(Na != Nb || Ma != Mb) { 
        std::cout << "Erorr subtraction" << std::endl;
        return result;
    }

    for(int i = 0; i < Na; ++i) {
        for(int j = 0; j < Ma; ++j) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result; 
}

//Транспонирование матриц
std::vector<std::vector<double>> matrixTransposition(const std::vector<std::vector<double>>& A) {
    int n = A.size();
    int m = A[0].size();
    std::vector<std::vector<double>> result(m, std::vector<double>(n, 0.0));
    
    for(int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            result[j][i] = A[i][j];
        }
    }
    
    return result;

}

//Метод простых итераций
std::vector<double> simpleIterations(const std::vector<std::vector<double>>& A, 
                                     const std::vector<double>& b, 
                                     int maxIterations = 100, 
                                     double precision = 1e-4) {
    int n = b.size();
    std::vector<double> x(n, 0.0);
    std::vector<double> x_new(n);

    for (int k = 0; k < maxIterations; ++k) {
        for (int i = 0; i < n; ++i) {
            x_new[i] = b[i];
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    x_new[i] -= A[i][j] * x[j];
                }
            }
            x_new[i] /= A[i][i];
        }

        double error = 0.0;
        for (int i = 0; i < n; ++i) {
            error += std::abs(x_new[i] - x[i]);
            x[i] = x_new[i];
        }

        if (error < precision) {
            break;
        }
    }

    return x;
}

//Поиск обратной матрицы
std::vector<std::vector<double>> matrixInverse(std::vector<std::vector<double>>& A,
                                               int maxIterations = 200,
                                               double precision = 0.00001) {
    int n = A.size();
    std::vector<std::vector<double>> X(n, std::vector<double>(n, 0.0));
    std::vector<double> y(n, 0.0);

    for(int i = 0; i < n; ++i) {
        X[i][i] = 1.0;
    }

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            if (i == j) {
                y[j] = 1;
            } else {
                y[j] = 0;
            }
        }
        y = simpleIterations(A, y, maxIterations, precision);
        for(int j = 0; j < n; ++j) {
            X[i][j] = y[j];
        }
    }

    return matrixTransposition(X);
}

//Определитель матрицы 
double matrixDeterminant (const std::vector<std::vector<double>> &A) {
    int n = A.size();
    int m = A[0].size();

    double result = 0;

    if(n != m) {   
        std::cout << "Erorr: matrixDeterminant, different size" << std::endl;
        return result;
    }

    if(n == 2) {
        result = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    }

    if(n == 3) {
        result = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2];
    }

    return result;
}

//Метод Гаусса
std::vector<std::vector<double>> gauss(std::vector<std::vector<double>> A,
                                       std::vector<std::vector<double>> B) {
    int n = A.size();
    std::vector<std::vector<double>> M(n, std::vector<double>(n + 1));

    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            M[i][j] = A[i][j];
        }
        M[i][n] = B[i][0];
    }                     

    for(int i = 0; i < n; ++i) {
        int pivot = i;
        for(int j = i + 1; j < n; ++j) {
            if (abs(M[j][i]) > abs(M[pivot][i])) {
                pivot = j;
            }
        }

        std::swap(M[i], M[pivot]);

        for (int j = i + 1; j < n; ++j) {
            double factor = M[j][i] / M[i][i];
            for (int k = i; k < n + 1; ++k) {
                M[j][k] -= factor * M[i][k];
            }
        }
    }

    std::vector<std::vector<double>> x(n, std::vector<double>(1, 0.0));
    for(int i = n - 1; i >= 0; --i) {
        x[i][0] = M[i][n] / M[i][i];
        for (int j = i + 1; j < n; ++j) {
            x[i][0] -= M[i][j] * x[j][0] / M[i][i]; 
        }
    }    

    return x;           
                  
}
