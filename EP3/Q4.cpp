// Importação de bibliotecas
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

struct Tuple{
    /*
    Esta estrutura serve como uma tupla (int, int) para 
    retorno fácil de índices de pontos de uma tabela ou intervalo
    */
public:
    int x;
    int y;

    Tuple(int x, int y){
        this->x = x;
        this->y = y;
    }
};

class DataTable{
    /* Essa classe foi feita para armazenar
    os valores do vetor de entrada e o tamanho do vetor
    de forma a deixar intuitivo o uso dela para busca
    do índice de valores e visar escalabilidade do código
    (adição dos valores de "y" ordenados com "x") */
public:
    
    vector<double> x_data; // Dados guardados no objeto
    vector<double> y_data;
    int size; // Quantidade de elementos nos dados

    DataTable(const vector<double>& x, const vector<double>& y){
        this->x_data = vector<double>(x);
        this->y_data = vector<double>(y);
        size = x.size();
    }

    Tuple findInterval(double value){
        /* 
        double -> (int, int)

        Esse método recebe um número decimal e retorna, 
        pelo método da bissecção, os índices do intervalo
        que contém este número em nossa tabela
        */
        
        int low = 0; // Índice do menor elemento do intervalo da busca
        int high = size - 1; // Índice do maior elemento

        // Caso o valor de entrada seja maior que o valor do maior elemento 
        // ou menor que o valor do menor elemento da tabela, retornamos os 
        // índices do intervalo nos extremos dos pontos correspondentes
        if(value > x_data[high]){
            cout << "Tentativa de extrapolação pela direita encontrada." << endl;
            return Tuple(high - 1, high);
        }
        else if(value < x_data[low]){
            cout << "Tentativa de extrapolação pela esquerda encontrada." << endl;
            return Tuple(low, low + 1);
        }

        // Realiza a bissecção dos índices até que o intervalo da bissecção tenha tamanho 1
        while(high - low > 1){
            int middle = (low + high) / 2; // Cálculo do valor do índice do meio (bissecção) do intervalo
            
            if(x_data[middle] > value) high = middle;
            else low = middle;
        }

        return Tuple(low, high); // Retorna os índices do intervalo encontrado
    }

    vector<int> getNeighbors(double value, Tuple interval, int quantity){
        /* 
        (double, (int, int), int) -> vector<int>

        Este método pega os "m" índices dos elementos vizinhos de um
        valor, considerando o intervalo que o cobre na tabela, centralizando
        os vizinhos pegos no valor de entrada.
        */

        vector<int> neighbors = { interval.x, interval.y }; // array de saída (começa com o intervalo mínimo, evitando extrapolação)
        int leftIdx = interval.x - 1; // índice do elemento à esquerda do intervalo atual considerado
        int rightIdx = interval.y + 1; // índice do elemento à direita do intervalo atual considerado

        for(int i = 0; i < quantity - 2; ++i){
            // Verifica se os índices à esquerda e direita estão dentro dos limites da tabela e 
            bool canGoLeft = (leftIdx >= 0);
            bool canGoRight = (rightIdx < size);
            if(!canGoLeft && !canGoRight) break; // caso os dois saiam dos limites, não há para onde ir mais

            if(canGoLeft && canGoRight){
                // checa se o valor mais à esquerda do intervalo atual é mais próximo do nosso valor de entrada
                // que o valor mais à direita
                if(abs(value - x_data[leftIdx]) < abs(value - x_data[rightIdx])){
                    // Insere o índice do valor à esquerda nos vizinhos
                    neighbors.insert(neighbors.begin(), leftIdx--); // Aumenta o intervalo considerado em 1 para a esquerda
                }
                else{
                    // Caso a condição acima não seja satisfeita faz o mesmo que a condição acima só que para a direita do intervalo
                    neighbors.push_back(rightIdx++);
                }

            } else if(canGoLeft){
                // Só consegue ir para a esquerda então aumenta o intervalo até onde der
                neighbors.insert(neighbors.begin(), leftIdx--);
            } else if(canGoRight){
                // Só consegue ir para a direita então aumenta o intervalo até onde der
                neighbors.push_back(rightIdx++);
            }

        }

        return neighbors;
    }

};


double Interpolate(DataTable data, double value, int interpOrder){
    /* 
    (DataTable, double, int) -> double

    Esta função recebe uma tabela de dados pré-ordenada, um valor de "x"
    para interpolar (ou extrapolar) e a ordem da interpolação e, então
    realiza a interpolação do polinômio pelo método de Lagrange,
    devolvendo o valor de f(x).
    */

    // Encontra o intervalo que contém x e avisa se houver extrapolação
    Tuple interval = data.findInterval(value); 
    // Identifica os "n+1" pontos para interpolação
    vector<int> neighbors = data.getNeighbors(value, interval, interpOrder + 1);

    double result = 0; // valor de saída f(x)
    int idx_i; // índice i obtido dos vizinhos
    int idx_j; // índice j obtido dos vizinhos

    for(int i = 0; i < interpOrder + 1; ++i){

        double thisL = 1; // Fator de Lagrange L_i
        idx_i = neighbors[i];
        for(int j = 0; j < interpOrder + 1; ++j){
            if(j == i) continue; // Pula esta iteração se j = i

            idx_j = neighbors[j];

            // Cálculo do fator de Lagrange L_i pelos valores de x, x_i e x_j (produtório)
            thisL *= (value - data.x_data[idx_j]) / (data.x_data[idx_i] - data.x_data[idx_j]);
        }

        result += data.y_data[idx_i] * thisL; // Soma dos componentes do polinômio de Lagrange (somatório)
    }

    return result;
}


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
        y_data.push_back(1/(1 + 25*x*x));
    }

    SplineSystem spline = SplineSystem(x_data, y_data);

    DataTable table = DataTable(x_data, y_data);

    // Tabela com os dados de execução do código, no total com 100 pontos entre -1 e 1
    // calculando os valores da função real, avaliados pela spline e por interpolação polinomial por Lagrange
    ofstream tableFile("q4.txt");
    tableFile << "x,f(x),p20,spline" << endl;
    for (int i = 0; i < 101; ++i){
        double xi = -1 + i * 0.02;
        double fx = 1/(1 + 25*xi*xi);
        double p20x = Interpolate(table, xi, 20);
        double sx = spline.evaluate(xi);
        tableFile << xi << "," << fx << "," << p20x << "," << sx << endl;
    }
    tableFile.close();

    return 0;
}