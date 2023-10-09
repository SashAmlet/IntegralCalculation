#include <iostream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <iomanip>
using namespace std;


#pragma region 1)

unordered_map<int, vector<pair<double, double>>>/*n, x, c*/ LEGENDRE_POLYNOMIALS_TABLE = {
        {1, {
                {0.0, 2.000000000000000}
            }},
        {2, {
                {0.577350269189626, 1.000000000000000},
                {-0.577350269189626, 1.000000000000000}
            }},
        {3, {
                {0.000000000000000, 0.888888888888889},
                {0.774596669241483, 0.555555555555556},
                {-0.774596669241483, 0.555555555555556}
            }},
        {4, {
                {0.339981043584856, 0.652145154862546},
                {0.861136311594053, 0.347854845137454},
                {-0.339981043584856, 0.652145154862546},
                {-0.861136311594053, 0.347854845137454}
            }},
        {5, {
                {0.000000000000000, 0.568888888888889},
                {0.538469310105683, 0.478628670499366},
                {0.906179845938664, 0.236926885056189},
                {-0.538469310105683, 0.478628670499366},
                {-0.906179845938664, 0.236926885056189}
            }},
        {6, {
                {0.238619186083197, 0.467913934572691},
                {0.661209386466265, 0.360761573048139},
                {0.932469514203152, 0.171324492379170},
                {-0.238619186083197, 0.467913934572691},
                {-0.661209386466265, 0.360761573048139},
                {-0.932469514203152, 0.171324492379170}
            }},
        {7, {
                {},
                {},
                {},
                {},
                {},
                {},
                {},
            }}
};
#pragma region old version
double Function__1_TO_1(double x) {
    if (x < -1.0 /*|| x == 1.0*/) {
        return NAN;
    }

    double result = 1.0 / (2.0 - 2.0 * log((x + 1.0) / 2.0));
    return result;
}

double GaussMethod(int NumOfTerms) {
    double result = 0.0;
    for (const auto& term : LEGENDRE_POLYNOMIALS_TABLE[NumOfTerms]) {
        result += term.second * Function__1_TO_1(term.first);
    }
    return result;
}
#pragma endregion


double Function__1_TO_1(double x, double dx, double xi) {
    if (x < -1.0 /*|| x == 1.0*/) {
        return NAN;
    }
    double c = dx / 2.0;
    double result = (c * (1.0 / (1.0 - log(xi + (x + 1.0) * c))));
    return result;
}

double GaussMethod(int NumOfTerms, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        for (const auto& term : LEGENDRE_POLYNOMIALS_TABLE[NumOfTerms]) {
            result += term.second * Function__1_TO_1(term.first, 1.0 / n, static_cast<double>(i) / n);
        }
    }
    return result;
}

void Output(int numOdDig, int polynomialOrder, int numOfFunc) {
    std::cout << std::fixed << std::setprecision(numOdDig); // Устанавливаем формат вывода для double

    std::cout << std::setw(15) << "N/N";
    for (int j = 1; j < numOfFunc; j = j * 10) {
        std::cout << std::setw(15) << j << ' ';
    }
    cout << endl;
    for (int i = 0; i < polynomialOrder; i++)
    {
        cout << "----------------";
    }
    cout << endl;
    for (int i = 1; i < polynomialOrder; i++) {
        std::cout << std::setw(15) << i << "| ";
        for (int j = 10; j < numOfFunc; j = j * 10) {
            std::cout << std::setw(15) << GaussMethod(i, j) << " ";
        }
        std::cout << std::endl;
    }
}



#pragma endregion

#pragma region 2)

double Function_0_TO_1(double x) {
    if (x < 0.0 /*|| x == 1.0*/) {
        return NAN;
    }

    double result = 1.0 / (1.0 - log(x));
    return result;
}

vector<vector<double>> RombergIntegration(int numOfSteps, const vector<double>* Xs) {
    vector<vector<double>> R(numOfSteps, vector<double>(numOfSteps, 0.0));
    double a = Xs->front(), b = Xs->back(), h = abs(b - a);
    int numOfTerms = 1;

    for (int n = 0; n < numOfSteps; n++) {
        for (int i = 0; i < numOfTerms; i++) {
            double tempA = a + i * h;
            R[n][0] += h * Function_0_TO_1((tempA + tempA + h) / 2);
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


#pragma endregion







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
    cout << ResultMatrix.back().back() << endl;


    Output(7, 7, pow(10, 7));

}