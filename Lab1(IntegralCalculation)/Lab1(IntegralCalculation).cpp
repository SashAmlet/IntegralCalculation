#include <iostream>
#include <cmath>
#include <vector>
using namespace std;


double Function(double x) {
    if (x < 0.0 /*|| x == 1.0*/) {
        // Функция не определена при x <= 0 или x == 1
        cerr << "Ошибка: Функция не определена для x <= 0 или x == 1." << endl;
        return NAN; // Возврат NaN (не число) в случае ошибки
    }

    double result = 1.0 / (1.0 - log(x));
    return result;
}

// Функция для вычисления таблицы Ромберга
vector<vector<double>> RombergIntegration(int numOfSteps, const vector<double>* Xs) {
    vector<vector<double>> R(numOfSteps, vector<double>(numOfSteps, 0.0));
    double a = Xs->front(), b = Xs->back(), h = abs(b - a);
    int numOfTerms = 1;

    for (int n = 0; n < numOfSteps; n++) {
        for (int i = 0; i < numOfTerms; i++) {
            double tempA = a + i * h;
            R[n][0] += h*Function((tempA + tempA + h) / 2);
        }

        for (int m = 1; m <= n; m++) {
            double factor = pow(4, m);
            R[n][m] = (factor * R[n][m - 1] - R[n - 1][m - 1]) / (factor - 1);
        }

        h /= 2.0;
        numOfTerms *= 2;
    }

    return R;
}



int main()
{
    // User's variables
    const int numOfSteps = 10;
    const double start = 0, finish = 1;
    //

    const double numOfVars = pow(2, numOfSteps - 1) + 1, intervalLength = abs(finish - start);
    const double delta = intervalLength / (numOfVars - 1); // except 0
    vector<double> Xs;
    vector<vector<double>> ResultMatrix(numOfSteps, vector<double>(numOfSteps));


    for (int i = 0; i < numOfVars; i++)
    {
        Xs.push_back(start + delta * i);
    }
    ResultMatrix = RombergIntegration(numOfSteps, &Xs);
    cout << ResultMatrix.back().back();
}