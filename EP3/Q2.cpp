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


int main(){
    // Valores de entrada da tabela
    vector<double> x_data = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
    vector<double> y_data = {-3.0f, -0.5f, -1.0f, 0.0f, 0.5f, 1.0f};

    DataTable table = DataTable(x_data, y_data);

    cout << "Valor interpolado para a tabela fornecida, para x = 3.2 e ordem 3: f(x) = " << Interpolate(table, 3.2f, 3) << endl;

    return 0;
}