// Importação de bibliotecas
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

class SplineSystem {
    /* 
    Classe para guardar as informações de uma spline, como
    elementos do sistema linear, dados de input e resultados das
    derivadas segundas do sistema e usar elas para a avaliação
    do valor de um ponto da função ajustada pela spline
     */
private:
    vector<double> x; // Dados de entrada de x
    vector<double> y; // Dados de entrada de f(x)
    vector<double> h; // x[i] - x[i-1], valores auxiliares
    vector<double> e; // 6 * (y[i] - y[i-1]) / h[i], valores auxiliares
    vector<double> phi; // Derivadas segundas para cada ponto x[i] de f(x)

    int size; // Quantidade de pontos usados
    int numDer; // Quantidade de derivadas segundas desconhecidas (size - 2)

    void build(){
        /* Constrói a spline, iniciando os valores de auxílio do sistema 
        ('h' e 'e') e calculando as derivadas segundas desconhecidas */

        // Vetores de auxílio do sistema, com tamanho igual a quantidade de pontos menos um
        h.assign(size - 1, 0.0);
        e.assign(size - 1, 0.0);
        for(int i = 0; i < size-1; ++i){
            h[i] = (x[i+1] - x[i]);
            e[i] = 6.0 * (y[i+1] - y[i]) / h[i];
        }

        // Número de derivadas segundas desconhecidas
        numDer = size - 2;

        // Vetor de derivadas segundas (já conhecemos a primeira e a última do vetor, iguais a zero)
        phi.assign(size, 0.0);
        SolveSpline();
    }

    void SolveSpline(){
        /* Resolve o sistema linear da Spline por decomposição LU,
        calculando os valores das derivadas segundas desconhecidas */

        // Arrays da decomposição LU ('h' + 'u' para matriz U e 'l' para matriz L), array 'z' solução de L
        // Aqui temos U composto pela diagonal principal 'u' e acima por 'h' e L composto pela diagonal abaixo por 'l'
        vector<double> u(numDer), l(numDer-1), z(numDer);

        // Cálculo da decomposição LU para os vetores 'l' e 'u' a partir dos valores de 'h'
        u[0] = 2*(h[0] + h[1]);
        int j; // Contador auxiliar para evitar confusão de índices (j = i+1)
        for(int i = 0; i < numDer-1; ++i){
            j = i+1;

            l[i] = h[j] / u[i];
            u[j] = 2*(h[j] + h[j+1]) - h[j]*l[i];
        }

        // Cálculo do vetor 'z' da solução do sistema do tipo L*z=b por substituição para frente
        z[0] = e[1] - e[0];
        for(int i = 1; i < numDer; ++i){
            z[i] = e[i+1] - e[i] - l[i - 1] * z[i - 1];
        }

        // Cálculo do vetor 'phi' de derivadas segundas (aqui calculamos só as desconhecidas e deixamos as
        // conhecidas com valor 0.0 padrão), por substituição para trás, do sistema de tipo U*phi=z
        int n = numDer - 1; // Último elemento dos resultados do sistema
        phi[n+1] = z[n] / u[n];
        for(int i = n - 1; i >= 0; --i){
            phi[i+1] = (z[i] - h[i+1] * phi[i+2]) / u[i]; // Preserva os valores de phi[0] e phi[size-1]
        }

    }

    int findInterval(double value){
        /* 
        double -> (int, int)

        Esse método recebe um número decimal e retorna, 
        pelo método da bissecção, os índices do intervalo
        que contém este número em nossa tabela
        */
        
        int low = 0; // Índice do menor elemento do intervalo da busca
        int high = size - 1; // Índice do maior elemento

        // Realiza a bissecção dos índices até que o intervalo da bissecção tenha tamanho 1
        while(high - low > 1){
            int middle = (low + high) / 2; // Cálculo do valor do índice do meio (bissecção) do intervalo
            
            if(x[middle] > value) high = middle;
            else low = middle;
        }

        return high; // Retorna o maior índice do intervalo encontrado
    }

    double cubicSplineInterp(int cur, double value){
        /* Fórmula padrão de interpolação por spline cúbica, usando
        os valores já conhecidos dos vetores de auxílio 'h', os pontos
        fornecidos de 'x' e 'y' e o vetor de derivadas segundas nestes pontos,
        para calcular um valor f(x) dado um x qualquer no intervalo da spline (value)
        e o índice do ponto da spline de maior deste intervalo (cur) */

        int prev = cur - 1; // índice anterior ao atual (cur)
        double deltaXi = x[cur] - value;
        double deltaX = value - x[prev];

        double S = (phi[prev] * pow(deltaXi, 3) + phi[cur] * pow(deltaX, 3)) / (6*h[prev])
                + ( y[prev] / h[prev] - h[prev] * phi[prev] / 6 ) * deltaXi
                + ( y[cur] / h[prev] - h[prev] * phi[cur] / 6 ) * (deltaX);

        return S;
    }

public:
    SplineSystem(const vector<double>& x_data, const vector<double>& y_data){
        // Construtor da classe copia os dados de input para o objeto e já constrói
        // a spline, deixando o usuário livre pra avaliar qualquer ponto dela depois
        x = vector<double>(x_data);
        y = vector<double>(y_data);
        size = x.size();
        build();
    }

    double evaluate(double value){
        /* Avalia f(x) para um dado valor 'x' para esta spline */

        int idx = findInterval(value); // índice do maior 'x' do intervalo em que este valor se encontra

        // Lida com os valores nos extremos da spline de forma explícita, evitando erros
        if (value == x[idx-1]) return y[idx-1];
        if (value == x[idx]) return y[idx];

        return cubicSplineInterp(idx, value); // Retorna o valor calculado pela interpolação de spline cúbica padrão
    }
};


int main(){
    // Valores de entrada da tabela
    vector<double> x_data, y_data;

    for(int i = 0; i < 21; ++i){
        double x = -1 + i * 0.1;
        x_data.push_back(x);
        y_data.push_back(1/(1 + 25*x*x)); // Valores da função real
    }

    SplineSystem spline = SplineSystem(x_data, y_data);

    // Tabela com os dados de execução do código, no total com 100 pontos entre -1 e 1
    // calculando os valores da função real e avaliados pela spline
    ofstream tableFile("q3.txt");
    tableFile << "x,f(x),spline" << endl;
    for (int i = 0; i < 101; ++i){
        double xi = -1 + i * 0.02;
        double fx = 1/(1 + 25*xi*xi);
        double sx = spline.evaluate(xi);
        tableFile << xi << "," << fx << "," << sx << endl;
    }
    tableFile.close();

    return 0;
}