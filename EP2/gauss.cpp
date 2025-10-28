#include <iostream>
#include <vector>
#include <ranges>
#include <cmath>
#include <string.h>
#include <experimental/simd> //Biblioteca com vetorização
using namespace std;
namespace stdx = std::experimental;
using simd_vec = stdx::native_simd<float>;

int batch_size = simd_vec::size();

// Algoritmo para resolver o inverso da raiz quadrada 4x mais 
// rápido que a operação comum (Quake III - Fast inverse square root)
float Q_rsqrt( float number )
{
	long i;
	float x2, y;
	const float threehalfs = 1.5F;

	x2 = number * 0.5F;
	y  = number;
	i  = *(long*) &y;                       // evil floating point bit level hacking
	i  = 0x5f3759df - ( i >> 1 );               // what the fuck?
	y  = *(float*) &i;
	y  = y *(threehalfs - ( x2 * y * y ) );   // 1st iteration - método de newton

	return y;
}

class LinearSystem{
public:
    vector<float> *data;
    int32_t size;

    vector<int> linePivot;
    vector<int> columnPivot;

    LinearSystem(vector<float> *thisData){
        data = thisData;
        size = sqrt(data->size());
        linePivot = vector<int>(size, -1);
        columnPivot = vector<int>(size, -1);
    }

    float& at(int i, int j){
        return (*data)[i * size + j];
    }

    void printLinearSystem() {
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                cout << this->at(i, j) << "\t";
            }
            cout << endl;
        }
        cout << "--------------------" << endl;
    }

    void normalize(){
        for(int j = 0; j < size; j++){
            simd_vec sqr_sum_vec(0.0f);

            int i = 0;
            for(; i + batch_size <= size; i += batch_size){
                float* row_j_ptr = &this->at(j, i);
                simd_vec rowVec;
                rowVec.copy_from(row_j_ptr, stdx::element_aligned);

                sqr_sum_vec += rowVec * rowVec;

                rowVec.copy_to(row_j_ptr, stdx::element_aligned);
            }

            float sqr_sum = stdx::reduce(sqr_sum_vec);

            for(; i < size; i++){
                sqr_sum += this->at(j, i) * this->at(j, i);
            }

            if(sqr_sum == 0){
                cout << "A row of this matrix has magnitude zero, it can't be normalized." << endl;
                return;
            }
            else if(sqr_sum == 1){
                // Row already normalized, continue to next row
                continue;
            }

            float invSqr = Q_rsqrt(sqr_sum);

            i = 0;
            for(; i + batch_size <= size; i += batch_size){
                float* row_j_ptr = &this->at(j, i);
                simd_vec rowVec;
                rowVec.copy_from(row_j_ptr, stdx::element_aligned);

                rowVec *= invSqr;

                rowVec.copy_to(row_j_ptr, stdx::element_aligned);
            }

            for(; i < size; i++){
                this->at(j, i) *= invSqr;
            }
        }
    }
};

float pivoting(LinearSystem &matrix, int index){

    float maxPivotLine = index;

    // Itera em todas as linhas para ver o maior valor da coluna index
    for(int i = index+1; i < matrix.size; ++i){
        if(abs(matrix.at(i, index)) > abs(matrix.at(maxPivotLine, index))){
            maxPivotLine = i;
        }
    }

    float maxPivotCol = index;

    // Itera em todas as colunas para ver o maior valor da linha index
    for(int j = index+1; j < matrix.size; ++j){
        if(abs(matrix.at(index, j)) > abs(matrix.at(index, maxPivotCol))){
            maxPivotCol = j;
        }
    }

    // Encerra o pivotamento total caso este seja o maior valor possivel
    if(maxPivotLine == index && maxPivotCol == index) return matrix.at(index, index);

    // Caso o maior valor das colunas seja maior que o maior valor das linhas analisadas,
    // troca as colunas, senão, troca as linhas
    if(abs(matrix.at(index, maxPivotCol)) > abs(matrix.at(maxPivotLine, index))){
        // Realiza a troca de colunas por iteração
        for(int i = 0; i < matrix.size; ++i){
            float helper = matrix.at(i, index);
            matrix.at(i, index) = matrix.at(i, maxPivotCol);
            matrix.at(i, maxPivotCol) = helper;
        }

        matrix.columnPivot[index] = maxPivotCol;
    }
    else{
        //Realiza a troca das linhas na memória
        float* helper = (float*)(malloc(sizeof(float) * matrix.size));

        memcpy(helper, &matrix.at(index, 0), sizeof(float) * matrix.size);
        memcpy(&matrix.at(index, 0), &matrix.at(maxPivotLine, 0), sizeof(float) * matrix.size);
        memcpy(&matrix.at(maxPivotLine, 0), helper, sizeof(float) * matrix.size);

        free(helper);

        matrix.linePivot[index] = maxPivotLine;
    }

    // Realiza o pivotamento novamente até que o maior valor total seja encontrado
    return pivoting(matrix, index);
}

void elimination(LinearSystem &matrix){
    int size = matrix.size;

    // Itera em todas as linhas para zerar todos os elementos da coluna do pivot abaixo da linha atual
    for(int i = 0; i < size-1; ++i){
        float pivot = pivoting(matrix, i);

        if (pivot == 0) { // Sistema singular ou com múltiplas soluções
            cout << "Erro: Pivô zero encontrado. O sistema pode não ter solução única." << endl;
            return;
        }

        // Itera em todas as linhas abaixo da atual
        for (int j = i + 1; j < size; ++j){
            float mult = -matrix.at(j, i) / pivot;

            simd_vec mult_vec(mult);

            // Realiza a eliminação de gauss em cada linha de forma parelizada com SIMD
            matrix.at(j, i) = mult;
            int k = i + 1;
            for(; k + batch_size <= size; k += batch_size){
                float* row_i_ptr = &matrix.at(i, k);
                float* row_j_ptr = &matrix.at(j, k);

                simd_vec rowVec_i;
                rowVec_i.copy_from(row_i_ptr, stdx::element_aligned);

                simd_vec rowVec_j;
                rowVec_j.copy_from(row_j_ptr, stdx::element_aligned);

                rowVec_j += rowVec_i * mult_vec;

                rowVec_j.copy_to(row_j_ptr, stdx::element_aligned);
            }

            for(; k < size; ++k){
                matrix.at(j, k) += matrix.at(i, k) * mult;
            }
        }
    }

    // Checa se a matriz é mal condicionada ou impossível
    // Primeiro normalizamos a matriz
    vector<float> normalizedData = vector<float>(&matrix.at(0,0), &matrix.at(size-1, matrix.size));

    LinearSystem normalizedLinearSystem(&normalizedData);
    normalizedLinearSystem.normalize();

    // Por ela ser triangular agora podemos calcular o determinante apenas multiplicando os elementos da diagonal
    // O determinante tem o sinal invertido dependendo do número de pivotamentos mas isso não importa em nada aqui pois usamos abs()
    float det = 1.0;
    for(int i = 0; i < size; ++i){
        det *= normalizedLinearSystem.at(i, i);
    }

    // Se o determinante for muito perto de zero (menor que 0.01) ela é mal-condicionada ou impossível
    if(abs(det) < 0.01f) {
        cout << "Matriz mal-condicionada ou impossível, os resultados podem não condizer com a realidade." << endl;
    }

}

vector<float> substitution(LinearSystem &matrix, vector<float> &coeff){
    int size = matrix.size;
    vector<float> result(size);

    // Arruma os coeficientes de acordo com os pivotamentos de linhas e multiplicadores
    for(int i = 0; i < size; ++i){
        if(matrix.linePivot[i] != -1){
            float helper = coeff[i];
            coeff[i] = coeff[matrix.linePivot[i]];
            result[matrix.linePivot[i]] = helper;
        }

        // Itera nas linhas abaixo de i, aplicando os multiplicadores guardados
        for(int j = i+1; j < size; ++j){
            coeff[j] += matrix.at(j, i) * coeff[i];
        }
    }

    // AQUI COMEÇA A SUBSTITUIÇÃO
    // Itera em todas as linhas, começando pela última
    for(int i = size - 1; i >= 0; --i){
        float sum = 0.0;

        // Itera em cada coluna após a diagonal atual
        for(int j = i + 1; j < size; ++j){
            sum += result[j] * matrix.at(i, j);
        }

        //Calcula o resultado final daquela linha, sendo ele o termo independente menos a soma divididos pelo coeficiente da diagonal
        result[i] = (coeff[i] - sum) / matrix.at(i, i);
    }

    // Troca as linhas do resultado de acordo com os pivotamentos de colunas
    for(int i = 0; i < size; ++i){
        if(matrix.columnPivot[i] != -1){
            float helper = result[i];
            result[i] = result[matrix.columnPivot[i]];
            result[matrix.columnPivot[i]] = helper;
        }
    }

    return result;
}

int main(){
    float input = 0;
    vector<float> matrixVector = {
        2, 1, 0, 1, 0, 0,
        1, 2, 1, 0, 1, 0,
        0, 1, 2, 0, 0, 1,
        1, 0, 0, 2, 1, 1,
        0, 1, 0, 1, 2, 1,
        1, 1, 1, 0, 1, 2
    };
    // Input do usuario da matriz principal
    // cout << "Insira os valores da sua matriz quadrada (não inclua termos independentes): " << endl;
    // cout << "Digite 'p' para parar." << endl;
    // while ((cin >> input) && input != 'p'){
    //     matrixVector.push_back(input);
    // }

    LinearSystem myLinearSystem = LinearSystem(&matrixVector);

    cout << "Matriz Original:" << endl;
    myLinearSystem.printLinearSystem();

    elimination(myLinearSystem);

    cout << "Matriz Triangular Superior:" << endl;
    myLinearSystem.printLinearSystem();

    input = 0;
    vector<float> coefficient = { 8, 13, 14, 20, 22, 23 };
    // Input do usuario do vetor de termos independentes
    // cout << "Insira os valores do seu vetor de termos independentes: " << endl;
    // cout << "Digite 'p' para parar." << endl;
    // while ((cin >> input) && input != 'p'){
    //     coefficient.push_back(input);
    // }

    vector<float> result = substitution(myLinearSystem, coefficient);

    cout << "Solucao:" << endl;
    for (size_t i = 0; i < result.size(); ++i) {
        cout << "x" << i << " = " << result[i] << endl;
    }

    return 0;
}