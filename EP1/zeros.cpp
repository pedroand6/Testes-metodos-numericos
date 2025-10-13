#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>

struct Result {
    double root;
    double error;
    int iter;
};

class ZeroFunction{
    private:
        double (*mathFunction)(double);
    public:
        ZeroFunction(double (*thisFunction)(double)){
            mathFunction = thisFunction;
        }

    double MidIntersection(double start, double end){
        return (start + end) / 2.0;
    }
    
    Result Bissection(double start, double end, double tolerance, int maxIter){
        Result result;

        if(mathFunction(start) * mathFunction(end) >= 0.0){
            std::cout << "Coudn't find any root in the interval" << std::endl;
            return result;
        }

        double intersection = 0.0;
        double thisError;

        for(int i = 0; i < maxIter; ++i){
            intersection = MidIntersection(start, end);

            thisError = (std::abs(end - start) / std::abs(end));
            if( thisError < tolerance ){
                result.iter = i;
                result.error = thisError;
                result.root = intersection;

                return result;
            }

            if(mathFunction(start) * mathFunction(intersection) < 0.0){
                end = intersection;
            }
            else{
                start = intersection;
            }
        }

        result.iter = maxIter;
        result.error = thisError;
        result.root = intersection;

        std::cout << "The function didn't converged" << std::endl;
        return result;

    }

    double FalsePointIntersection(double start, double end){
        return start - mathFunction(start) * (end - start) / (mathFunction(end) - mathFunction(start));
    }

    Result FalsePosition(double start, double end, double tolerance, int maxIter){
        Result result;

        if(mathFunction(start) * mathFunction(end) >= 0.0){
            std::cout << "Coudn't find any root in the interval" << std::endl;
            return result;
        }

        double intersection = 0.0;
        double lastIntersection = NAN;
        double thisError;

        for(int i = 0; i < maxIter; ++i){
            intersection = FalsePointIntersection(start, end);

            if( lastIntersection != NAN ) {

                thisError = std::abs(lastIntersection - intersection) / std::abs(lastIntersection);
                if( thisError < tolerance ){
                    
                    result.iter = i;
                    result.error = thisError;
                    result.root = intersection;

                    return result;
                }

            }

            if(mathFunction(start) * mathFunction(intersection) < 0.0){
                end = intersection;
            }
            else{
                start = intersection;
            }

            lastIntersection = intersection;
        }

        result.iter = maxIter;
        result.error = thisError;
        result.root = intersection;

        std::cout << "The function didn't converged" << std::endl;
        return result;

    }
};

double myFunc(double x){
    return cos(x) - x;
}

int main(){
    std::cout << std::fixed \
              << std::setprecision(std::numeric_limits<double>::max_digits10);

    ZeroFunction blob(&myFunc);
    Result bissecResResult = blob.Bissection(0.0, 1.0, 1E-15, 10000);
    Result falsePosResult = blob.FalsePosition(0.0, 1.0, 1E-15, 10000);

    std::cout << bissecResResult.root << std::endl;

    std::cout << falsePosResult.root << std::endl;

    return 0;
}