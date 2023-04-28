#include <iostream>
#include <math.h>

static class SqMatrixCalculator
{
public:

	typedef double** Matrix;
	typedef double* MRow;
	typedef double MElement;

	enum {
		UPPER,
		LOWER
	};



	static MRow IterationsMethod(Matrix matrix, MRow additions, int size, double accuracy) {

		MRow xi = new MElement[size]; // storage for x_i;
		for (int i = 0; i < size; i++) xi[i] = 0;

		MRow deltaxi = new MElement[size]; //storage for accuracy
		int iterations = 0;
		while (1) {
			iterations++;
			for (int i = 0; i < size; i++) {
				deltaxi[i] = xi[i];
				xi[i] = additions[i];
				for (int j = 0; j < size; j++) {
					if (j == i) continue;
					xi[i] -= matrix[i][j] * xi[j];
				}
				xi[i] /= matrix[i][i];
				deltaxi[i] = abs(deltaxi[i] - xi[i]);
			}

			MElement max_delta = DBL_MIN;
			for (int i = 0; i < size; i++) {
				if (deltaxi[i] > max_delta) max_delta = deltaxi[i];
			}
			if (max_delta < accuracy) {
				//std::cout <<std::endl<< iterations << std::endl;
				return xi;
			}
		}

	}




	static MRow GaussMethod(Matrix matrix, MRow additions, int size) {
		Matrix tMatrix;
		tMatrix = copyMatrix(matrix, size);
		for (int coloumn = 0; coloumn < size; coloumn++) {

			for (int string = coloumn + 1; string < size; string++) {
				MElement koeff = tMatrix[string][coloumn] / tMatrix[coloumn][coloumn];
				sumStrings(tMatrix[string], tMatrix[coloumn], -koeff, size);
				additions[string] -= additions[coloumn] * koeff;
				tMatrix[string][coloumn] = 0;
			}
		}
		return solveByTriangular(tMatrix, additions, size, UPPER);
	}



	static MRow CholeskyMethod(Matrix matrix, MRow addition, int size) {
		Matrix LMatrix = new MRow[size];

		for (int i = 0; i < size; i++) {
			LMatrix[i] = new MElement[size];
			for (int j = 0; j < size; j++) {
				LMatrix[i][j] = 0;
			}
		}



		LMatrix[0][0] = sqrt(matrix[0][0]);

		for (int string = 0; string < size; string++) {
			//L[i][0]			
			LMatrix[string][0] = matrix[string][0] / LMatrix[0][0];

			//L[i][j], j<i
			for (int coloumn = 1; coloumn < string; coloumn++) {

				double buffer = matrix[string][coloumn];
				for (int k = 0; k < coloumn; k++) {
					buffer -= (LMatrix[coloumn][k] * LMatrix[string][k]);
				}
				LMatrix[string][coloumn] = buffer / LMatrix[coloumn][coloumn];

			}

			// L[i][i] 

			double buffer = matrix[string][string];
			for (int k = 0; k < string; k++) {
				buffer -= pow(LMatrix[string][k], 2);
			}
			LMatrix[string][string] = sqrt(buffer);


		}
		// Solving


		MRow results = solveByTriangular(TransposeMatrix(LMatrix, size), solveByTriangular(LMatrix, addition, size, LOWER), size, UPPER);


		return results;

	}



























	static MElement Determine(Matrix matrix, int size) {
		bool sign = 1;
		Matrix triangularMatrix = Triangulation(matrix, size, &sign);
		MElement determinator = 1;

		for (int string = 0; string < size; string++) {
			determinator *= triangularMatrix[string][string];
		}
		delete[] triangularMatrix;
		return sign == 1 ? determinator : -determinator;
	}


	static Matrix* LUFactorization(Matrix matrix, int size) {
		Matrix LMatrix = new MRow[size];
		Matrix UMatrix = new MRow[size];

		for (int i = 0; i < size; i++) {
			LMatrix[i] = new MElement[size];
			UMatrix[i] = new MElement[size];
		}

		//preparing the L-Matrix by stating it's upper part as zeros and diagonal elements as 1

		for (int i = 0; i < size; i++) {
			LMatrix[i][i] = 1;
			for (int j = i + 1; j < size; j++) {
				LMatrix[i][j] = 0;
			}
		}
		// Triangulation using Gauss Method with the step-by-step elements pushing into L-matrix

		UMatrix = copyMatrix(matrix, size);

		for (int coloumn = 0; coloumn < size; coloumn++) {
			for (int string = coloumn + 1; string < size; string++) {
				MElement koeff = UMatrix[string][coloumn] / UMatrix[coloumn][coloumn];
				
				sumStrings(UMatrix[string], UMatrix[coloumn], -koeff, size);
				LMatrix[string][coloumn] = koeff;
				UMatrix[string][coloumn] = 0;
			}
		}

		return new Matrix[2]{LMatrix, UMatrix};

	}



	static MRow solveByTriangular(Matrix matrix, MRow additions, int size, short type = UPPER) {
		
		MRow results = new MElement[size];
		for (int i = 0; i < size; i++) results[i] = additions[i];

		if (type == LOWER) {

			for (int i = 0; i < size; i++) {
				for (int j = 0; j < i; j++) {
					results[i] -= matrix[i][j] * results[j];
				}
				results[i] /= matrix[i][i];
			}

		}

		if (type == UPPER) {

			for (int i = size - 1; i >= 0; i--) {
				for (int j = size - 1; j > i; j--) {
					results[i] -= matrix[i][j] * results[j];
				}
				results[i] /= matrix[i][i];
			}

		}

		return results;

	}




	static Matrix TransposeMatrix(Matrix matrix, int size) {
		
		Matrix transposedM = new MRow[size];

		for (int i = 0; i < size; i++) {
			transposedM[i] = new MElement[size];
			for (int j = 0; j < size; j++) {
				transposedM[i][j] = matrix[j][i];
			}
		}
		return transposedM;
	}




	static Matrix MultiplyMatrix(Matrix one, Matrix two, int size) {
		Matrix result = new MRow[size];

		//result[0][0] = one[0][0] * two[0][0] + one[0][1] * two[1][0];
		//result[0][1] = one[0][0] * two[0][1] + one[0][1] * two[1][1]

		for (int i = 0; i < size; i++)
		{
			result[i] = new MElement[size];
			for (int j = 0; j < size; j++)
			{
				result[i][j] = 0;
				for (int k = 0; k < size; k++)
					result[i][j] += one[i][k] * two[k][j];
			}
		}


		return result;
	}






	static Matrix Triangulation(Matrix matrix, int size, bool* sign = nullptr) {

		Matrix tMatrix;
			tMatrix = copyMatrix(matrix, size);
			for (int coloumn = 0; coloumn < size; coloumn++) {

				for (int string = coloumn + 1; string < size; string++) {
					MElement koeff = tMatrix[string][coloumn] / tMatrix[coloumn][coloumn];
					sumStrings(tMatrix[string], tMatrix[coloumn], -koeff, size);
					tMatrix[string][coloumn] = 0;
				}
			}
		return tMatrix;
	}




private:
	static Matrix copyMatrix(Matrix matrix, int size) {
		Matrix newMatrix = new MRow[size];
		for (int i = 0; i < size; i++) {
			newMatrix[i] = new MElement[size];
			for (int j = 0; j < size; j++) {
				newMatrix[i][j] = matrix[i][j];
			}
		}
		return newMatrix;
	}

	static int findRowWithMaxElem(Matrix matrix, int size, int from, int coloumn = -1) {
		int Row_i = 0;
		MElement maxEl = -std::numeric_limits<MElement>::infinity();

		if (coloumn < 0 || coloumn >= size)
			for (int string = from; string < size; string++) {
				for (int col = 0; col < size; col++) {
					if (abs(matrix[string][col]) > maxEl) {
						Row_i = string;
						maxEl = abs(matrix[string][col]);
					}
				}
			}

		if (coloumn >= 0 && coloumn < size)
			for (int string = from; string < size; string++) {
				if (abs(matrix[string][coloumn]) > maxEl) {
					Row_i = string;
					maxEl = abs(matrix[string][coloumn]);
				}
			}

		return Row_i;

	}

	static void swapStrings(MRow& a, MRow& b) {
		double* temp = a;
		a = b;
		b = temp;
	}

	static void sumStrings(MRow& sum, MRow& summand, MElement k, int size) {
		if (sum != summand) {
			for (int i = 0; i < size; i++) {
				sum[i] += (summand[i] * k);
			}
		}
	}


};
