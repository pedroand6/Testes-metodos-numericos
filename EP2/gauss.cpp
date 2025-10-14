#include <iostream>
#include <vector>
#include <cmath>
#include <string.h>
#include <experimental/simd> //Biblioteca com vetorização
using namespace std;
namespace stdx = std::experimental;
using simd_vec = stdx::native_simd<double>;

class Matrix{
public:
    vector<double> *data;
    int32_t nlin;
    int32_t ncol;

    Matrix(vector<double> *thisData, int lin, int col){
        data = thisData;
        nlin = lin;
        ncol = col;
    }

    double& at(int i, int j){
        return (*data)[i * ncol + j];
    }
};

void printMatrix(Matrix& matrix) {
    for (int i = 0; i < matrix.nlin; ++i) {
        for (int j = 0; j < matrix.ncol; ++j) {
            cout << matrix.at(i, j) << "\t";
        }
        cout << endl;
    }
    cout << "--------------------" << endl;
}

double pivoting(Matrix &matrix, int index){

    double maxPivot = index;

    // Itera em todas as linhas para ver o maior pivot da coluna index
    for(int i = index+1; i < matrix.nlin; ++i){
        if(abs(matrix.at(i, index)) > abs(matrix.at(i, maxPivot))){
            maxPivot = i;
        }
    }

    // Troca as linhas se achou um pivot diferente do inicial
    if(maxPivot != index){
        double* helper = (double*)(malloc(sizeof(double) * matrix.ncol));

        memcpy(&matrix.at(0, index), helper, sizeof(double) * matrix.ncol);
        memcpy(&matrix.at(0, maxPivot), &matrix.at(0, index), sizeof(double) * matrix.ncol);
        memcpy(&helper, &matrix.at(0, maxPivot), sizeof(double) * matrix.ncol);

        free(helper);
    }

    return matrix.at(index, index);
}

bool isBadConditioned(vector<vector<double>> &matrix){
    return false;
}

void elimination(Matrix &matrix){

    int nlin = matrix.nlin;

    // Itera em todas as linhas para zerar todos os elementos da coluna do pivot abaixo da linha atual
    for(int i = 0; i < nlin; ++i){
        double pivot = pivoting(matrix, i);

        if (pivot == 0) { // Sistema singular ou com múltiplas soluções
            cout << "Erro: Pivô zero encontrado. O sistema pode não ter solução única." << endl;
            return;
        }

        // Itera em todas as linhas abaixo da atual
        for (int j = i + 1; j < nlin; ++j){
            double mult = -matrix.at(j, i) / pivot;

            // Itera em todas as colunas da linha atual, começando pela
            // coluna da diagonal atual, pois antes dela foi tudo zerado
            for(int k = i; k < nlin + 1; ++k){
                matrix.at(j, k) += matrix.at(i, k) * mult;
            }
        }
    }

}

vector<double> substitution(Matrix &matrix){

    int nlin = matrix.nlin;
    vector<double> result(nlin);

    // Itera em todas as linhas, começando pela última
    for(int i = nlin - 1; i >= 0; --i){
        double sum = 0.0;

        // Itera em cada coluna após a diagonal atual (tirando o termo independente)
        for(int j = i + 1; j < nlin; ++j){
            sum += result[j] * matrix.at(i, j);
        }

        //Calcula o resultado final daquela linha, sendo o termo independente menos a soma divididos pelo coeficiente da diagonal
        result[i] = (matrix.at(i, nlin) - sum) / matrix.at(i, i);
    }

    return result;
}

int main(){

    vector<double> matrixVector = {
        10,     -7,     0,  7, 
        -3,     2.099,  6,  3.901, 
        5,      -1,     5,  6
    };

    Matrix matrix = Matrix(&matrixVector, 3, 4);

    cout << "Matriz Original:" << endl;
    printMatrix(matrix);

    elimination(matrix);

    cout << "Matriz Triangular Superior:" << endl;
    printMatrix(matrix);

    vector<double> result = substitution(matrix);

    cout << "Solucao (x, y, z):" << endl;
    for (size_t i = 0; i < result.size(); ++i) {
        cout << "x" << i << " = " << result[i] << endl;
    }

    return 0;
}