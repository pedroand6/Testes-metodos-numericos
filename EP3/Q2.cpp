// Importação de bibliotecas
#include <iostream>
#include <vector>
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

    DataTable(vector<double> x_data, vector<double> y_data, int size){
        this->x_data = x_data;
        this->y_data = y_data;
        this->size = size;
    }

    Tuple findInterval(double value){
        /* 
        double -> (int, int)

        Esse método recebe um número decimal e retorna, 
        pelo método da bissecção, os índices do intervalo
        que contém este número em nossa tabela
        */
        
        int lowerIdx = 0; // Índice do menor elemento do intervalo da busca
        int higherIdx = this->size - 1; // Índice do maior elemento
        int middleIdx = higherIdx / 2; // Índice do elemento do "meio"

        // Caso o valor de entrada seja maior que o valor do maior elemento 
        // ou menor que o valor do menor elemento da tabela, retornamos os 
        // índices do intervalo nos extremos dos pontos correspondentes
        if(value > x_data[higherIdx]){
            cout << "Tentativa de extrapolação pela direita encontrada." << endl;
            return Tuple(higherIdx - 1, higherIdx);
        }
        else if(value < x_data[lowerIdx]){
            cout << "Tentativa de extrapolação pela esquerda encontrada." << endl;
            return Tuple(lowerIdx, lowerIdx + 1);
        }

        // Realiza a bissecção dos índices até que o intervalo da bissecção tenha tamanho 1
        while(higherIdx - lowerIdx > 1){
            
            if(this->x_data[middleIdx] > value){
                higherIdx = middleIdx;
            }
            else{
                lowerIdx = middleIdx;
            }

            // Cálculo do valor do índice do meio (bissecção) do intervalo
            middleIdx = (higherIdx + lowerIdx) / 2;
        }

        return Tuple(lowerIdx, higherIdx); // Retorna os índices do intervalo encontrado
    }

    vector<int> getNeighbors(double value, Tuple interval, int quantity){
        /* 
        (double, (int, int), int) -> vector<int>

        Este método pega os "m" índices dos elementos vizinhos de um
        valor, considerando o intervalo que o cobre na tabela, centralizando
        os vizinhos pegos no valor de entrada.
        */

        vector<int> neighbors = {}; // array de saída
        int leftIdx = interval.x; // índice do elemento à esquerda do intervalo atual considerado
        int rightIdx = interval.y; // índice do elemento à direita do intervalo atual considerado

        for(int i = 0; i < quantity; ++i){
            // Verifica se o índice à esquerda está dentro dos limites da tabela e checa se
            // o valor mais à esquerda do intervalo atual é mais próximo do nosso valor de entrad
            // que o valor mais à direita
            if(leftIdx >= 0 && abs(value - this->x_data[leftIdx]) < abs(value - this->x_data[rightIdx])){
                neighbors.insert(neighbors.begin(), leftIdx); // Insere o índice do valor à esquerda nos vizinhos
                leftIdx--; // Aumenta o intervalo considerado em 1 para a esquerda

                continue; // Segue para o próximo elemento, ignorando a parte seguinte reservada para o elemento à direita
            }

            // Caso a condição acima não seja satisfeita E o índice à direita do intervalo esteja
            // dentro dos limites da tabela, faz o mesmo que a condição acima só que para a direita do intervalo
            if(rightIdx < this->size){
                neighbors.push_back(rightIdx);
                rightIdx++;
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


int main(){
    // Valores de entrada da tabela
    vector<double> x_data = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
    vector<double> y_data = {-3.0f, -0.5f, -1.0f, 0.0f, 0.5f, 1.0f};

    DataTable table = DataTable(x_data, y_data, x_data.size());

    cout << Interpolate(table, 3.2f, 3) << endl;

    return 0;
}