// Lab1-2_Matrix_Determination.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <stdio.h>
#include <random>
#include <fstream>
#include "SqMatrixCalculator.h"
#include <time.h>
#include <chrono>

int main()
{

    int n = 3;

   
    double Matrix[3][3] = { {6.25,-1,0.5},{-1,5,2.12},{0.5,2.12,3.6} };
    double** Matrix_ = new double* [3];
    for (int i = 0; i < 3; i++) {
        Matrix_[i] = new double[3];
        for (int j = 0; j < 3; j++) {
            Matrix_[i][j] = Matrix[i][j];
        }
    }

    double addition[3] = { 7.5,-8.68,-0.24 };
    double* addition_ = new double[3];
    for (int i = 0; i < 3; i++) {
        addition_[i] = addition[i];
    }
    
    auto start_time = std::chrono::steady_clock::now();
    double* mur = SqMatrixCalculator::CholeskyMethod(Matrix_, addition_, n);
    auto end_time = std::chrono::steady_clock::now();
    std::cout << "\n\n";
    std::cout << "Cholesky method: " << "\n";
    for (int i = 0; i < n; i++) {
            std::cout << mur[i] << " ";
        }
    std::cout << "\n Used time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();

    start_time = std::chrono::steady_clock::now();
    double* meow = SqMatrixCalculator::GaussMethod(Matrix_, addition_, n);
    end_time = std::chrono::steady_clock::now();
    std::cout << "\n\n";
    std::cout << "Gauss method: " << "\n";
    for (int i = 0; i < n; i++) {
        std::cout << meow[i] << " ";
    }
    std::cout << "\n Used time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();

    std::cout << "\n\n";
    std::cout << "Iterations method: " << "\n";
    start_time = std::chrono::steady_clock::now();
    double* meowkaMurka = SqMatrixCalculator::IterationsMethod(Matrix_, addition, n, 0.00000001);
    end_time = std::chrono::steady_clock::now();
    for (int i = 0; i < n; i++) {
        std::cout << meowkaMurka[i] << " ";
    }
    std::cout << "\n Used time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();

}

