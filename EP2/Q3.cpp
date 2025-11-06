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

class LinearSystem{
// Classe da matriz principal do sistema linear para gerenciamento e print dos dados
public:
    vector<float> *data; // Ponteiro dos dados da matriz
    int32_t size; // Largura da matriz

    LinearSystem(vector<float> *thisData){
        data = thisData;
        size = sqrt(data->size());
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
};

// Função auxiliar para pegar o maior elemento de um array
float MaxElement(float list[], int size){
    float max = list[0];

    for(int i = 0; i < size; i++){
        if(list[i] > max){
            max = list[i];
        }
    }
    
    return max;
}

// Função para determinar se o critério de Sassenfeld é satisfeito por uma matriz
bool Sassenfeld(LinearSystem &matrix){
    int size = matrix.size;
    float betas[size]; // Array de valores de Beta

    for(int i = 0; i < size; i++){
        float thisBeta = 0;

        for(int j = 0; j < i; j++){
            thisBeta += betas[j] * abs(matrix.at(i, j));
        }

        for(int j = i+1; j < size; j++){
            thisBeta += abs(matrix.at(i, j));
        }

        betas[i] = thisBeta / abs(matrix.at(i, i));
    }

    if(MaxElement(betas, size) < 1) return true;

    return false;
}

// Função para calcular as soluções do sistema linear pelo método de Gauss-Seidel
int GaussSeidel(LinearSystem &matrix, vector<float> &coeff, vector<float> &result){
    int size = matrix.size;
    vector<float> nextResult = vector<float>(result); // Vetor de resultados adiante que são calculados durante a iteração

    float convergence = 1.0f;
    float convergenceFactors[size]; // Valores para a determinação da convergência

    // Definição de um máximo de iterações para evitar loops infinitos com situações de não-convergência
    const int maxIterations = 10000;
    int iterations = 0;

    // Continua o loop do método até que o critério de convergência seja satisfeito
    while(convergence >= 1e-5 && iterations < maxIterations){
        // Copia os dados dos resultados adiante para os resultados passados 
        // (não serão usados até o cálculo convergência)
        memcpy(result.data(), nextResult.data(), sizeof(float) * size); 
        
        // Itera em cada linha da matriz
        for(int i = 0; i < size; i++){
            float thisValue = coeff[i]; // Calcula a soma total do valor da linha

            // Itera em cada coluna menos na diagonal
            for(int j = 0; j < size; j++){
                if(j == i) continue;
                thisValue -= matrix.at(i, j) * nextResult[j]; // Vai usando valores sempre atualizados dos resultados
            }
            // Atualiza os resultados no mesmo momento para a próxima iteração
            nextResult[i] = thisValue / matrix.at(i, i); // Pode dar ruim quando o valor da diagonal for muito perto de zero
        }

        // Calcula o valor do critério de convergência (erro relativo) usando vetorização SIMD
        int i = 0;
        for(; i + batch_size <= size; i += batch_size){
            simd_vec nextResultVec;
            simd_vec lastResultVec;
            simd_vec convergenceVec;
            float* nextResult_ptr = &nextResult[i];
            float* lastResult_ptr = &result[i];
            
            nextResultVec.copy_from(nextResult_ptr, stdx::element_aligned);
            lastResultVec.copy_from(lastResult_ptr, stdx::element_aligned);

            convergenceVec = abs((nextResultVec - lastResultVec) / nextResultVec);

            convergenceVec.copy_to(&convergenceFactors[i], stdx::element_aligned);
        }

        // Calcula o mesmo critério acima para os elementos que não couberam no SIMD
        for(; i < size; i++){
            convergenceFactors[i] = abs((nextResult[i] - result[i]) / nextResult[i]);
        }

        convergence = MaxElement(convergenceFactors, size);
        iterations++;
    }

    return iterations;
}


int main(){
    float input = 0;
    vector<float> matrixVector;

    // vector<float> matrixVector = {
    //     10,     -2, -2,  1,
    //     -2,      5, -1, -1,
    //      1,   0.5f, -6,  1,
    //     -1,     -1,  0, 20
    // };

    // Input do usuario da matriz principal
    int matCount;
    cout << "Insira qual a largura da sua matriz quadrada principal: ";
    cin >> matCount;

    cout << "Insira os valores da sua matriz quadrada (não inclua termos independentes): " << endl;
    for(int i = 0; i < matCount*matCount; i++){
        cin >> input;
        matrixVector.push_back(input);
    }

    cout << "--------------------" << endl;

    LinearSystem myLinearSystem = LinearSystem(&matrixVector);

    cout << "Matriz Original:" << endl;
    myLinearSystem.printLinearSystem();

    // Verificação se a matriz satisfaz o critério de Sassenfeld
    bool doesConverge = Sassenfeld(myLinearSystem);

    if(doesConverge){
        cout << "O critério de Sassenfeld foi satisfeito." << endl;
    }
    else{
        cout << "O critério de Sassenfeld não foi satisfeito, não pode se esperar que os resultados convirjam." << endl;
    }

    cout << "--------------------" << endl;

    vector<float> coefficient;
    //vector<float> coefficient = { 3, 5, -9, 17 };

    // Input do usuario do vetor de termos independentes
    cout << "Insira os valores do seu vetor de termos independentes: " << endl;
    for(int i = 0; i < myLinearSystem.size; i++){
        cin >> input;
        coefficient.push_back(input);
    }

    cout << "--------------------" << endl;

    vector<float> results;
    //vector<float> results = { 100, 100, 100, 100 };

    // Input do usuario do vetor de termos independentes
    cout << "Insira o chute inicial para as respostas: " << endl;
    for(int i = 0; i < myLinearSystem.size; i++){
        cin >> input;
        results.push_back(input);
    }

    cout << "--------------------" << endl;
    
    cout << "Insira o chute inicial para as respostas: " << endl;
    for(int i = 0; i < myLinearSystem.size; i++){
        cout << results[i] << endl;
    }

    cout << "--------------------" << endl;

    // Soluções
    int iterations = GaussSeidel(myLinearSystem, coefficient, results);

    cout << "A Solução convergiu em " << iterations << " iterações: " << endl;
    for (size_t i = 0; i < results.size(); ++i) {
        cout << "x" << i << " = " << results[i] << endl;
    }

    return 0;
}