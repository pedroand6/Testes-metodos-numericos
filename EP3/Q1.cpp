// Importação de bibliotecas
#include <iostream>
#include <vector>
using namespace std;

class DataTable{
    /* Essa classe foi feita para armazenar
    os valores do vetor de entrada e o tamanho do vetor
    de forma a deixar intuitivo o uso dela para busca
    do índice de valores e visar escalabilidade do código */
public:
    
    vector<double> data; // Dados guardados no objeto
    int size; // Quantidade de elementos nos dados

    DataTable(vector<double> data, int size){
        this->data = data;
        this->size = size;
    }

    int findIndex(double value){
        /* 
        (double) -> int

        Esse método recebe um número decimal e retorna, 
        pelo método da bissecção, o índice com o valor
        mais próximo à esse número em nossa tabela 
        (arredondando o valor para baixo)  */
        
        int low = 0; // Índice do menor elemento do intervalo da busca
        int high = size - 1; // Índice do maior elemento

        // Caso o valor de entrada seja maior ou igual que o valor do maior elemento 
        // da tabela, retornamos o índice deste elemento
        if(value >= data[high]) return high;

        // Realiza a bissecção dos índices até que o intervalo da bissecção tenha tamanho 1
        while(high - low > 1){
            int middle = (low + high) / 2; // Cálculo do valor do índice do meio (bissecção) do intervalo
            
            if(data[middle] > value) high = middle;
            else low = middle;
        }

        return low; // Retorna o menor índice encontrado
    }

};

int main(){
    // Valores de entrada (retirados do Ex1 do EP1)
    vector<double> data = {0.0027f, 0.0059f, 0.0492f, 0.054f, 0.0886f, 0.1421f, 
        0.2362f, 0.2777f, 0.3426f, 0.3926f, 0.5386f, 0.5736f, 0.6649f, 0.6915f, 
        0.7763f, 0.7793f, 0.8335f, 0.869f, 0.9172f, 0.9383f};

    DataTable myTable = DataTable(data, data.size());

    cout << "Esta tabela tem 20 elementos aleatórios entre 0,0027 e 0,9383." << endl;
    
    cout << "Índice esperado para valor abaixo do menor elemento da tabela: 0" << endl;
    cout << "Índice encontrado: " << myTable.findIndex(-1.0f) << endl;

    cout << "----------------------------------------" << endl;

    cout << "Índice esperado para valor igual a 0,5386: 10" << endl;
    cout << "Índice encontrado: " << myTable.findIndex(0.5386f) << endl;

    cout << "----------------------------------------" << endl;

    cout << "Índice esperado para valor acima do maior elemento da tabela: 19" << endl;
    cout << "Índice encontrado: " << myTable.findIndex(1.0f) << endl;

     return 0;
}