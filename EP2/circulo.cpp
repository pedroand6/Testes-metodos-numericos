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

enum ConditionState {
    ILL_CONDITIONED,
    INDETERMINED,
    WELL_CONDITIONED
};

ConditionState systemState;

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
// Classe da matriz principal do sistema linear para gerenciamento e print dos dados
public:
    vector<float> *data; // Ponteiro dos dados da matriz
    int32_t size; // Largura da matriz

    vector<int> linePivot; // Vector com os valores do índice referente a cada linha, ex: [3, 2, 1, 4] -> trocou linhas 1 e 3 no pivotamento
    vector<int> columnPivot; // Vector com os valores do índice referente a cada coluna

    LinearSystem(vector<float> *thisData){
        data = thisData;
        size = sqrt(data->size());
        linePivot.resize(size);
        columnPivot.resize(size);

        // Inserindo índices padrão das linhas e colunas
        for(int i = 0; i < size; ++i){
            linePivot[i] = i;
            columnPivot[i] = i;
        }
    }

    // Método para pegar valor dos dados com os índices da matriz
    float& at(int i, int j){
        return (*data)[i * size + j];
    }

    // Output dos dados da matriz em formatação de tabela
    void printLinearSystem() {
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                cout << this->at(i, j) << "\t";
            }
            cout << endl;
        }
        cout << "--------------------" << endl;
    }

    // Método para normalizar os vetores das linhas da matriz (usado na validação do sistema linear)
    void normalize(){
        for(int j = 0; j < size; j++){
            simd_vec sqr_sum_vec(0.0f);

            int i = 0;
            for(; i + batch_size <= size; i += batch_size){
                float* row_j_ptr = &this->at(j, i);
                simd_vec rowVec;
                rowVec.copy_from(row_j_ptr, stdx::element_aligned);

                sqr_sum_vec += rowVec * rowVec;
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

// Função de pivotamento total da matrix principal
float pivoting(LinearSystem &matrix, int index){
    int n = matrix.size;
    int maxPivotLine = index;
    int maxPivotCol = index;
    float maxAbs = abs(matrix.at(index, index));

    // Itera em todas as linhas e colunas para ver qual o maior valor absoluto da 
    // submatriz abaixo da linha e da coluna atual
    for(int i = index; i < n; ++i){
        // Itera em todas as colunas para ver o maior valor da linha index
        for(int j = index; j < n; ++j){
            float val = abs(matrix.at(i, j));
            if(val > maxAbs){
                maxAbs = val;
                maxPivotLine = i;
                maxPivotCol = j;
            }
        }
    }

    // Troca as colunas necessárias e as linhas necessárias para deixar o maior valor no lugar do valora atual
    if(maxPivotCol != index){
        // Realiza a troca de colunas por iteração
        for(int i = 0; i < matrix.size; ++i){
            float helper = matrix.at(i, index);
            matrix.at(i, index) = matrix.at(i, maxPivotCol);
            matrix.at(i, maxPivotCol) = helper;
        }

        // Troca os elementos referentes às trocas
        float helper = matrix.columnPivot[index];
        matrix.columnPivot[index] = matrix.columnPivot[maxPivotCol];
        matrix.columnPivot[maxPivotCol] = helper;
    }

    if(maxPivotLine != index){
        //Realiza a troca das linhas na memória
        float* helper = (float*)(malloc(sizeof(float) * matrix.size));

        memcpy(helper, &matrix.at(index, 0), sizeof(float) * matrix.size);
        memcpy(&matrix.at(index, 0), &matrix.at(maxPivotLine, 0), sizeof(float) * matrix.size);
        memcpy(&matrix.at(maxPivotLine, 0), helper, sizeof(float) * matrix.size);

        free(helper);

        // Troca os elementos referentes às trocas
        float lineHelper = matrix.linePivot[index];
        matrix.linePivot[index] = matrix.linePivot[maxPivotLine];
        matrix.linePivot[maxPivotLine] = lineHelper;
    }

    return matrix.at(index, index);
}

// Função da etapa de eliminação de Gauss da matriz principal
void elimination(LinearSystem &matrix){
    int size = matrix.size;

    // Itera em todas as linhas para zerar todos os elementos da coluna do pivot abaixo da linha atual
    for(int i = 0; i < size; ++i){
        float pivot = pivoting(matrix, i);

        if (pivot == 0) { // Sistema singular ou com múltiplas soluções
            systemState = INDETERMINED; // Erro: Pivô zero encontrado. O sistema pode não ter solução única.
            return;
        }

        // Itera em todas as linhas abaixo da atual
        for (int j = i + 1; j < size; ++j){
            float mult = matrix.at(j, i) / pivot;

            simd_vec mult_vec(mult);

            // Realiza a eliminação de gauss em cada linha de forma parelizada com SIMD
            matrix.at(j, i) = mult;
            int k = i + 1;
            for(; k + batch_size <= size; k += batch_size){
                float* row_i_ptr = &matrix.at(i, k);
                float* row_j_ptr = &matrix.at(j, k);

                simd_vec rowVec_i;
                simd_vec rowVec_j;

                rowVec_i.copy_from(row_i_ptr, stdx::element_aligned);
                rowVec_j.copy_from(row_j_ptr, stdx::element_aligned);

                rowVec_j -= rowVec_i * mult_vec;

                rowVec_j.copy_to(row_j_ptr, stdx::element_aligned);
            }

            for(; k < size; ++k){
                matrix.at(j, k) -= matrix.at(i, k) * mult;
            }
        }
    }

    // Checa se a matriz é mal condicionada ou impossível
    // Primeiro normalizamos a matriz
    vector<float> normalizedData = *(matrix.data);

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
        systemState = ILL_CONDITIONED; // Erro: Matriz mal-condicionada ou impossível, os resultados podem não condizer com a realidade
    }

}

// Função da etapa de substituição, que pode ser reutilizada de um mesmo objeto LinearSystem para diversos vetores de termos independentes
vector<float> substitution(LinearSystem &matrix, vector<float> &coeff){
    int size = matrix.size;
    vector<float> result(size);

    vector<float> realCoeff(size);
    // Arruma os coeficientes de acordo com os pivotamentos de linhas e multiplicadores
    for(int i = 0; i < size; ++i){
        realCoeff[i] = coeff[matrix.linePivot[i]];
    }

    for(int i = 0; i < size; ++i){
        // Itera nas linhas abaixo de i, aplicando os multiplicadores guardados
        for(int j = i+1; j < size; ++j){
            realCoeff[j] -= matrix.at(j, i) * realCoeff[i];
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

        float diagonal = matrix.at(i, i);

        if(diagonal == 0){
            systemState = INDETERMINED; // Encontrada diagonal com valor zero - Sistema indeterminado
            result[i] = 0.0f;
            continue;
        }

        //Calcula o resultado final daquela linha, sendo ele o termo independente menos a soma divididos pelo coeficiente da diagonal
        result[i] = (realCoeff[i] - sum) / diagonal;
    }

    vector<float> realResult(size);
    // Troca as linhas do resultado de acordo com os pivotamentos de colunas
    for(int i = 0; i < size; ++i){
        realResult[matrix.columnPivot[i]] = result[i];
    }

    return realResult;
}

// Função para calcular a matriz principal do problema do círculo e o vetor de termos independentes
vector<float> CircleMatrix(float thisPoints[3][2], vector<float> &linearTerms){

    vector<float> outMatrix;

    for(int i = 0; i < 3; i++){
        outMatrix.push_back(2 * thisPoints[i][0]); // 2 * x_i
        outMatrix.push_back(2 * thisPoints[i][1]); // 2 * y_i
        outMatrix.push_back(-1);

        linearTerms.push_back(thisPoints[i][0] * thisPoints[i][0] + thisPoints[i][1] * thisPoints[i][1]); // x_i ^ 2 + y_i ^ 2
    }

    return outMatrix;
}

int main(){

    systemState = WELL_CONDITIONED;
    string stateMessages[3] = {
        "Your system is ill conditioned, the results may not be precise.",
        "Your system is indetermined, it may contain multiple solutions.",
        "Your system is well conditioned."
    };

    // Input dos valores dos pontos
    float input_x, input_y;
    float points[3][2];
    
    cout << "Insira os valores dos pontos: " << endl;
    for(int i = 0; i < 3; i++){
        cout << "Ponto " << i + 1 << ": " << endl;
        cout << "x_" << i + 1 << ": ";
        cin >> input_x;
        points[i][0] = input_x;

        cout << "y_" << i + 1 << ": ";
        cin >> input_y;
        points[i][1] = input_y;
    }

    vector<float> linCoeff;
    vector<float> matrixVector = CircleMatrix(points, linCoeff);

    // Resolução do sistema linear
    LinearSystem myLinearSystem = LinearSystem(&matrixVector);

    elimination(myLinearSystem);

    vector<float> result = substitution(myLinearSystem, linCoeff);
    float a = result[0];
    float b = result[1];
    float r = sqrt(a*a + b*b - result[2]);

    cout << "--------------------" << endl;
    cout << stateMessages[systemState] << endl;
    cout << "--------------------" << endl;

    cout << "Solução:" << endl;
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "r = " << r << endl;

    return 0;
}