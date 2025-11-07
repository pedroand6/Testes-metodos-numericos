// Importação de bibliotecas
#include <iostream>
#include <math.h>
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

struct TridiagonalMat {
    vector<double> upperDiag;
    vector<double> middleDiag;
    vector<double> lowerDiag;

    TridiagonalMat() {}

    TridiagonalMat(int size){
        upperDiag = vector<double>(size-1, 0);
        middleDiag = vector<double>(size, 0);
        lowerDiag = vector<double>(size-1, 0);
    }
};

class SplineSystem {
private:
    vector<double> x_data;
    vector<double> y_data;

    vector<double> h_data;
    vector<double> e_data;
    vector<double> phi;

    int size;
    TridiagonalMat LU_Matrix;


    void DecomposeSystem(TridiagonalMat &matrix){
        int n = this->size - 1;
        matrix.middleDiag[0] = 2 * (this->h_data[0] + this->h_data[1]);
        
        for(int i = 0; i < n - 2; ++i){
            matrix.lowerDiag[i] = this->h_data[i+1] / matrix.middleDiag[i];
            matrix.middleDiag[i+1] = 2*(this->h_data[i+1] + this->h_data[i+2]) - (this->h_data[i+1] * this->h_data[i+1]) / matrix.middleDiag[i];
            //matrix.upperDiag[i - 1] = this->h_data[i]; // desnecessario para os calculos
        }
    }

    void SolveSpline(){
        int n = this->size - 1;

        vector<double> z = vector<double>(n - 1, 0);
        z[0] = this->e_data[1] - this->e_data[0];

        for(int i = 1; i < n - 1; ++i){
            z[i] = (this->e_data[i+1] - this->e_data[i]) - this->LU_Matrix.lowerDiag[i-1] * z[i-1];
        }

        this->phi[n - 1] = z[n - 2] / this->LU_Matrix.middleDiag[n - 2];

        for(int i = n - 2; i >= 1; --i){
            this->phi[i] = (z[i] - this->h_data[i+1] * this->phi[i+1]) / this->LU_Matrix.middleDiag[i];
        }

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
        if(value > this->x_data[higherIdx]){
            cout << "Tentativa de extrapolação pela direita encontrada." << endl;
            return Tuple(higherIdx - 1, higherIdx);
        }
        else if(value < this->x_data[lowerIdx]){
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

    double polynomial(int idx, double x){

        return (this->phi[idx - 1] * pow(x_data[idx] - x, 3) + this->phi[idx] * pow(x - this->x_data[idx-1], 3)) / (6 * this->h_data[idx])
                + ( this->y_data[idx-1] / this->h_data[idx] - this->h_data[idx] * this->phi[idx-1] / 6 ) * (this->x_data[idx] - x) 
                + ( this->y_data[idx] / this->h_data[idx] - this->h_data[idx] * this->phi[idx] / 6 ) * (x - this->x_data[idx-1]);

    }

public:
    SplineSystem(vector<double> x_data, vector<double> y_data){
        this->x_data = x_data;
        this->y_data = y_data;
        this->size = x_data.size();

        for(int i = 1; i < this->size; ++i){
            double h = x_data[i] - x_data[i-1];
            this->h_data.push_back(h);
            this->e_data.push_back(6 * (this->y_data[i] - this->y_data[i-1]) / h);
        }

        this->LU_Matrix = TridiagonalMat(this->size - 2);
        this->DecomposeSystem(this->LU_Matrix);

        this->phi = vector<double>(this->size, 0);
        this->SolveSpline();
    }

    double Evaluate(double value){
        Tuple interval = findInterval(value);
        return polynomial(interval.y, value);
    }

};


int main(){
    // Valores de entrada da tabela
    vector<double> x_data = {};
    vector<double> y_data = {};

    for(int i = 0; i < 21; ++i){
        double x = -1 + i * 0.1f;
        x_data.push_back(x);
        y_data.push_back(1/(1 + 25*x*x));
    }

    SplineSystem spline = SplineSystem(x_data, y_data);

    cout << spline.Evaluate(0.5) << endl; // Valor esperado: 0.137931

    return 0;
}