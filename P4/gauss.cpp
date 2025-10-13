#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

void printMatrix(const vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << val << "\t";
        }
        cout << endl;
    }
    cout << "--------------------" << endl;
}

double pivoting(vector<vector<double>> &matrix, int index){

    double maxPivot = index;

    // Itera em todas as linhas para ver o maior pivot da coluna index
    for(int i = index+1; i < matrix.size(); ++i){
        if(abs(matrix[i][index]) > abs(matrix[i][maxPivot])){
            maxPivot = i;
        }
    }

    // Troca as linhas se achou um pivot diferente do inicial
    if(maxPivot != index){
        vector<double> lastLine = matrix[index];
        matrix[index] = matrix[maxPivot];
        matrix[maxPivot] = lastLine;
    }

    return matrix[index][index];
}

void elimination(vector<vector<double>> &matrix){

    int size = matrix.size();

    // Itera em todas as linhas para zerar todos os elementos da coluna do pivot abaixo da linha atual
    for(int i = 0; i < size; ++i){
        double pivot = pivoting(matrix, i);

        if (pivot == 0) { // Sistema singular ou com múltiplas soluções
            cout << "Erro: Pivô zero encontrado. O sistema pode não ter solução única." << endl;
            return;
        }

        // Itera em todas as linhas abaixo da atual
        for (int j = i + 1; j < size; ++j){
            double mult = -matrix[j][i] / pivot;

            // Itera em todas as colunas da linha atual, começando pela
            // coluna da diagonal atual, pois antes dela foi tudo zerado
            for(int k = i; k < size + 1; ++k){
                matrix[j][k] += matrix[i][k] * mult;
            }
        }
    }

}

vector<double> substitution(vector<vector<double>> &matrix){

    int size = matrix.size();
    vector<double> result(size);

    // Itera em todas as linhas, começando pela última
    for(int i = size - 1; i >= 0; --i){
        double sum = 0.0;

        // Itera em cada coluna após a diagonal atual (tirando o termo independente)
        for(int j = i + 1; j < size; ++j){
            sum += result[j] * matrix[i][j];
        }

        //Calcula o resultado final daquela linha, sendo o termo independente menos a soma divididos pelo coeficiente da diagonal
        result[i] = (matrix[i][size] - sum) / matrix[i][i];
    }

    return result;
}

int main(){

    vector<vector<double>> matrix = {
        {10, -7, 0, 7}, 
        {-3, 2.099, 6, 3.901}, 
        {5, -1, 5, 6}
    };

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