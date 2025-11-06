// Importação de bibliotecas
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include <functional>

using namespace std;

// Definição de constantes físicas
const double h = 6.62607015e-34; // Constante de Planck (J·s)
const double c = 2.99792458e8; // Velocidade da luz no vácuo (m/s)
const double k = 1.380649e-23; // Constante de Boltzmann (J/K)
const double SUN_RADIUS = 6.957e8; // Raio do Sol (m)
const double PI = 3.14159265358979323846; // Valor de Pi

double BlackBody(double wavelength, double temperature){
    /* Função que calcula a radiação de corpo
    negro para uma determinada temperatura em um
    determinado comprimento de onda */

    // Fórmula de Planck para radiação de corpo negro
    double numerator = 2.0 * h * c * c;
    double denominator = pow(wavelength, 5.0) * (exp((h * c) / (wavelength * k * temperature)) - 1.0);
    return numerator / denominator;
}

class DustLayerStar{
    /* Classe de estrela com camada externa de poeira,
    contendo todos os atributos necessários para calcular
    sua densidade numérica de grãos da camada no Sistema Internacional */

    private:
        double starRadius;
        double starTemperature;
        double internRadius;
        double externRadius;
        double dustRadius;
        double dustTemperature;
        double colorIndex;
    public:
        // Construtor padrão
        DustLayerStar(double thisStarRadius, double thisStarTemperature, double thisInternRadius,
                            double thisExternRadius, double thisDustRadius, double thisDustTemperature, double thisColorIndex){

            starRadius = thisStarRadius;
            starTemperature = thisStarTemperature;
            internRadius = thisInternRadius;
            externRadius = thisExternRadius;
            dustRadius = thisDustRadius;
            dustTemperature = thisDustTemperature;
            colorIndex = thisColorIndex;
        }

    double DustLayerComponent(double wavelength, double n){
        /* Função que calcula a contribuição da camada de poeira
        para a radiação total da estrela em um determinado
        comprimento de onda e densidade de grãos n */

        double dustLayerFactor = (4.0 * PI * pow(dustRadius, 2.0) / 3.0) * (pow(externRadius, 3.0) - pow(internRadius, 3.0));

        return pow(starRadius, 2.0) * BlackBody(wavelength, starTemperature) + dustLayerFactor * BlackBody(wavelength, dustTemperature) * n;
    }

    double DustLayerFunction(double n){
        /* Função que relaciona o índice de cor observado com as
        características físicas da estrela e da poeira */

        return colorIndex + 2.5 * log10(DustLayerComponent(1.6E-6, n) / DustLayerComponent(2.2E-6, n));
    }
};

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

        cout << "The function didn't converged" << endl;
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

        cout << "The function didn't converged" << endl;
        return Result(thisError, input.maxIter, intersection);

    }
};

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

    // Criação do objeto estrela com camada de poeira com os dados do enunciado
    DustLayerStar star(20.0*SUN_RADIUS, 3000.0, 300.0*SUN_RADIUS, 1500.0*SUN_RADIUS, 0.2E-6, 1000.0, 0.5);

    // Criação da função "de envelope" para passar a função da classe como parâmetro
    auto starFunction = [&](double n) {
        return star.DustLayerFunction(n);
    };

    InputData input(starFunction, 0.0, 1000.0, 1E-5, 1000); // Dados de entrada
    ZeroFunction functionResolution(input); // Cria o objeto da função com os dados de entrada

    cout << "Método da Bissecção:" << endl;
    ShowResults(functionResolution.Bissection());
    cout << "------------------------------------" << endl;
    cout << "Método da Posição Falsa:" << endl;
    ShowResults(functionResolution.FalsePosition());

    return 0;
}
