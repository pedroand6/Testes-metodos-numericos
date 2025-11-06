// Importação de bibliotecas
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include <functional>

using namespace std;

struct InputData {
    /* Estrutura para armazenar os dados de entrada */

    function<double(double)> mathFunction;
    double start;
    double end;
    double tolerance;
    int maxIter;

    // Construtores:

    // Construtor padrão usado pra dentro da classe ZeroFunction
    InputData() : mathFunction(nullptr), start(0.0), end(0.0), tolerance(0.0), maxIter(0) {}

    InputData(function<double(double)> myFunc, double s, double e, double tol, int maxIt){
        mathFunction = myFunc;
        start = s;
        end = e;
        tolerance = tol;
        maxIter = maxIt;
    }
};

struct Result {
    /* Estrutura para armazenar os resultados */

    double root;
    double error;
    int iter;

    // Construtor para criar um resultado
    Result(double thisError, int thisIter, double thisRoot){
        root = thisRoot;
        error = thisError;
        iter = thisIter;
    }
};

class ZeroFunction{
    /* Classe que implementa os dois métodos para encontrar zeros de funções */

    private:
        InputData input; // Dados de entrada
    public:
        ZeroFunction(InputData thisInput){
            input = thisInput;
        }

    double SetError(double startInterval, double endInterval){
        /* Função que calcula o erro relativo ou absoluto dependendo se
         o valor do denominador do erro relativo for muito pequeno ou não */

        // O nextafter é o menor número representável maior que o valor dado (aqui em direção ao infinito)
        if( endInterval < 10 * nextafter(endInterval, INFINITY) ) {
            return (abs(endInterval - startInterval));
        }
        else{
            return (abs(endInterval - startInterval) / abs(endInterval));
        }
    }

    void UpdateInterval(double* ptrStart, double* ptrEnd, double intersection){
        /* Atualiza o intervalo do método de acordo com o Bolzano, respeitando
            a intersecção calculada */

        if(input.mathFunction(*ptrStart) * input.mathFunction(intersection) < 0.0){
            *ptrEnd = intersection;
        }
        else{
            *ptrStart = intersection;
        }
    }

    double MidIntersection(double start, double end){
        /* Calcula a intersecção do método
            da bissecção (meio do intervalo) */

        return (start + end) / 2.0;
    }

    Result Bissection(){
        /* Função do método da bisseção, retorna um resultado.
            Funciona dividindo o intervalo sempre na metade,
            verificando a raiz pelo método de Bolzano */

        // Inicialização de variáveis locais
        double startInterval = input.start;
        double endInterval = input.end;
        double intersection;
        double thisError;

        // Verificação do Teorema de Bolzano inicial
        if(input.mathFunction(startInterval) * input.mathFunction(endInterval) >= 0.0){
            cout << "Couldn't find any root in the interval" << endl;
            return Result(NAN, 0, NAN);
        }

        for(int i = 0; i < input.maxIter; ++i){
            intersection = MidIntersection(startInterval, endInterval);
            thisError = SetError(startInterval, endInterval);

            // Verificação de convergência
            if( thisError < 2*input.tolerance ){
                return Result(thisError, i, intersection);
            }

            UpdateInterval(&startInterval, &endInterval, intersection);
        }

        cout << "The function didn't converge" << endl;
        return Result(thisError, input.maxIter, intersection);
    }

    double FalsePointIntersection(double start, double end){
        /* Calcula a intersecção do método da posição
            falsa, seguindo a geometria da função */

        return start - input.mathFunction(start) * (end - start) / (input.mathFunction(end) - input.mathFunction(start));
    }

    Result FalsePosition(){
        /* Função do método da posição falsa, retorna um resultado.
            Funciona dividindo o intervalo pela interseção da reta
            que liga os pontos do intervalo */

        // Inicialização de variáveis locais
        double startInterval = input.start;
        double endInterval = input.end;
        double intersection;
        double lastIntersection = NAN;
        double thisError;

        // Verificação do Teorema de Bolzano inicial
        if(input.mathFunction(startInterval) * input.mathFunction(endInterval) >= 0.0){
            cout << "Couldn't find any root in the interval" << endl;
            return Result(NAN, 0, NAN);
        }

        for(int i = 0; i < input.maxIter; ++i){
            intersection = FalsePointIntersection(startInterval, endInterval);

            // Verificação de convergência (apenas se já houver uma interseção anterior)
            if( !isnan(lastIntersection) ) {

                thisError = SetError(intersection, lastIntersection);

                if( thisError < input.tolerance ){
                    return Result(thisError, i, intersection);
                }
            }

            UpdateInterval(&startInterval, &endInterval, intersection);
            lastIntersection = intersection; // Salve a última interseção
        }

        cout << "The function didn't converge" << endl;
        return Result(thisError, input.maxIter, intersection);

    }
};

// Primeira função de exemplo para o programa
// Esta função tem uma raiz em cerca de 0.7390851332
double myFunc(double x){
    return cos(x) - x;
}

// Segunda função de exemplo para o programa
// Esta função tem mais de uma raiz, uma em cerca de 0.0 e outra em cerca de 0.25
double myFunc2(double x){
    return cos(x) + sin(x) - 49*pow(x, 2) + log10(256*x);
}

InputData AskInputData(){
    /* Função que pergunta os dados de entrada para o usuário */

    // Inicialização de variáveis locais
    double (*selectedFunc)(double) = nullptr;
    double start, end, tolerance;
    int maxIter;

    // Pergunta qual das funções o usuário quer usar
    int funcChoice = 0;
    while(funcChoice != 1 && funcChoice != 2){
        cout << "Choose your function: \n \
        1 -> cos(x) - x \n \
        2 -> cos(x) + sin(x) - 49*x^2 + log(256*x) \n";
        cin >> funcChoice;
    }

    switch (funcChoice)
    {
        case 1:
            selectedFunc = &myFunc;
            break;
        case 2:
            selectedFunc = &myFunc2;
            break;
    }

    cout << "Enter the start of the interval: ";
    cin >> start;

    cout << "Enter the end of the interval: ";
    cin >> end;

    cout << "Enter the tolerance: ";
    cin >> tolerance;

    cout << "Enter the maximum number of iterations: ";
    cin >> maxIter;

    InputData input(selectedFunc, start, end, tolerance, maxIter);
    return input;
}

void ShowResults(Result res){
    /* Função que mostra os resultados requisitados no output */

    cout << "Root: " << res.root << endl;
    cout << "Error: " << res.error << endl;
    cout << "Iterations: " << res.iter << endl;
}

int main(){
    // Configuração para mostrar o máximo de casas decimais possíveis no output
    cout << fixed \
              << setprecision(numeric_limits<double>::max_digits10);

    while (true)
    {
        InputData input = AskInputData();
        ZeroFunction functionResolution(input); // Cria o objeto da função com os dados de entrada

        // Pergunta qual método o usuário quer usar
        int method = 0;
        while(method != 1 && method != 2){
            cout << "Choose the method (1 - Bissection, 2 - False Position): ";
            cin >> method;
        }

        // Chama o método escolhido e mostra os resultados
        if(method == 1){
            ShowResults(functionResolution.Bissection());
        }
        else{
            ShowResults(functionResolution.FalsePosition());
        }

        int tryAgain = -1;
        while(tryAgain != 0 && tryAgain != 1){
            cout << "Do you want to try again? (1 - Yes, 0 - No): ";
            cin >> tryAgain;
        }

        if(tryAgain == 0) break;
        cout << "------------------------------------" << endl;
    }

    return 0;
}
