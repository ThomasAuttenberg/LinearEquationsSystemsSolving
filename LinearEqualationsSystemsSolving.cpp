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

    int n = 4;

    // double Matrix[3][3] = { {6.25,-1,0.5},{-1,5,2.12},{0.5,2.12,3.6} };
   
    double Matrix[4][4] = { {99,14,-15,23},{16,128,-22, 29},{18,20,-167,32}, {10,12,-16,100} };
   //double Matrix[4][4] = { {20,2,3,7},{1,12,-2, -5},{5,-3,13,0}, {0,0,-1,15} };
    double** Matrix_ = new double*[n];
    for (int i = 0; i < n; i++) {
        Matrix_[i] = new double[n];
        for (int j = 0; j < n; j++) {
            Matrix_[i][j] = Matrix[i][j];
        }
    }

    // double addition[3] = { 7.5,-8.68,-0.24 };
    double addition[4] = { 5,8,9,4 }; //для метод
    //double addition[4] = { 5,4,-3,7 };
    double* addition_ = new double[n];
    for (int i = 0; i < n; i++) {
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
    double* meowkaMurka = SqMatrixCalculator::IterationsMethod(Matrix_, addition, n, 0.0000001);
    end_time = std::chrono::steady_clock::now();
    for (int i = 0; i < n; i++) {
        std::cout << meowkaMurka[i] << " ";
    }
    std::cout << "\n Used time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();

}

