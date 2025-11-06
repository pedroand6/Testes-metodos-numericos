// Importação de bibliotecas
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>

using namespace std;

class RandomData {
    /* Classe que gera um vetor de números aleatórios e os ordena
        pelo método de inserção */

    private:
        vector<float> vec;
    public:
        RandomData(int size) {
            vec = vector<float>(size);

            for (int i = 0; i < size; ++i) {
                vec[i] = (rand() % 10000) / 10000.0; // Gera numero aleatorio entre 0 e 1 com precisão de 4 casas decimais
            }
        }

    int SortElements(){
        /* Função que ordena os elementos do vetor pelo método de inserção
            e retorna o número de passos feitos */

        int npassos = 0; // Passos de memória executados pelo método

        //Ordenação por insercão
        for (int i = 1; i < vec.size(); ++i) {
            double tempVal = vec[i];
            npassos++; // Adiciona um passo de leitura do valor de vec[i]

            int j = i - 1;

            while (j >= 0 && vec[j] > tempVal) {
                npassos++; // Adiciona um passo de leitura do valor de vec[j]

                // Trova os elementos de lugar
                vec[j + 1] = vec[j];
                npassos++; // Adiciona um passo de escrita do valor de vec[j+1]

                j = j - 1; // Move para o próximo elemento à esquerda.
            }

            // Coloca o valor antigo na posição certa
            vec[j + 1] = tempVal;
            npassos++; // Adiciona um passo de escrita do valor de vec[j+1]
        }

        return npassos;
    }

    void PrintElements(){
        /* Função que imprime os elementos do vetor */
        for (int i = 0; i < vec.size(); ++i) {
            cout << "Element " << i << " = " << vec[i] << endl;
        }
    }
};

int main() {
    RandomData data(20);
    data.PrintElements();
    int npassos = data.SortElements();

    cout << "------------------------------------" << endl;

    data.PrintElements();
    cout << "Nº de passos: " << npassos << endl;

    //Tabela com os dados de execução do código
    ofstream tableFile("q1.txt");
    tableFile << "Npassos,Nelementos";
    for (int i = 10; i <= 100; i = i+5){
        RandomData newData(i);
        tableFile << "\n" << newData.SortElements() << ",";
        tableFile << i;
    }
    tableFile.close();

    return 0;
}
