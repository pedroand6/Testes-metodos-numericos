#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

class RandomData {
    private:
        std::vector<float> vec;
    public:
        RandomData(uint size) {
            vec = std::vector<float>(size);

            for (int i = 0; i < size; ++i) {
                vec[i] = (std::rand() % 10000) / 10000.0; // Gera numero aleatorio entre 0 e 1
            }
        }

    uint SortElements(){
        uint npassos = 0;

        //Ordenacao por insercao
        for (int i = 0; i < vec.size(); ++i) {
            for (int j = 0; j < i; ++j) {
                if (vec[i] < vec[j]){
                    float tempValue = vec[i];
                    vec[i] = vec[j];
                    vec[j] = tempValue;

                    npassos = npassos + 6;
                }
            }
        }

        return npassos;
    }

    void PrintElements(){
        for (int i = 0; i < vec.size(); ++i) {
            std::cout << "Element " << i << " = " << vec[i] << std::endl;
        }
    }
};

int main() {
    RandomData data(20);
    data.PrintElements();
    uint npassos = data.SortElements();

    std::cout << "------------------------------------" << std::endl;

    data.PrintElements();
    std::cout << "Nº de passos: " << npassos << std::endl;

    //Tabela com os dados de execucao do codigo
    for (int i = 10; i <= 100; i = i+5){
        RandomData newData(i);
        std::cout << newData.SortElements() << " -- " << i << std::endl;
    }

    return 0;
}