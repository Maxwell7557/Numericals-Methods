#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

class ZeidelAndJakoby{
private:
    int systemSize;
    double** coeff;
    double* freeTerm;
    double eps;
    string filename;
public:
    ZeidelAndJakoby();
    ~ZeidelAndJakoby();
    bool converge(double*,double*);
    void initializeFromFile();
    void change();
    void showSystem();
    void applyZeidelMethod();
    void applyJakobyMethod();
};

ZeidelAndJakoby::ZeidelAndJakoby()
{
    eps = 0.001;
    filename = "mat.txt";
}

ZeidelAndJakoby::~ZeidelAndJakoby()
{
    for (int i = 0; i < systemSize; i++)
        delete [] coeff[i];
    delete [] coeff;
    delete [] freeTerm;
}

bool ZeidelAndJakoby::converge(double *x, double *p)
{
    double* diff = new double[systemSize];
    double maxVal = 0;
    for (int i = 0; i < systemSize; i++)
    {
        diff[i] = fabs(x[i] - p[i]);
        if (diff[i] > maxVal)
            maxVal = diff[i];
    }
    return (maxVal<eps);
}

void ZeidelAndJakoby::initializeFromFile()
{
    ifstream fout(filename);

    fout >> systemSize;

    coeff = new double*[systemSize];
    freeTerm = new double[systemSize];

    for (int i = 0; i < systemSize; i++)
    {
        coeff[i] = new double[systemSize+1];
        for (int j = 0; j < systemSize; j++)
        {
            fout >> coeff[i][j];
        }

        fout >> freeTerm[i];
    }

    fout.close();
}

void ZeidelAndJakoby::showSystem()
{
    cout << "System: \n";
    for (int i = 0; i < systemSize; i++)
    {
        for (int j = 0; j < systemSize; j++)
        {
            if (coeff[i][j] < 0)
                cout << "(" << coeff[i][j] << ")*x" << j;
            else
                cout << coeff[i][j] << "*x" << j;
            if (j < systemSize - 1)
                cout << " + ";
        }
        cout << " = " << freeTerm[i] << endl;
    }
    return;
}

void ZeidelAndJakoby::applyZeidelMethod()
{
    cout << "--------ZEIDEL--------" << endl;

    double* p = new double[systemSize];
    double* x = new double[systemSize];

    for (int j = 0; j<systemSize; j++)
        x[j] = 1;

    int iterations = 0;
    do
    {
        for (int i = 0; i < systemSize; i++)
            p[i] = x[i];
        for (int i = 0; i < systemSize; i++)
        {
            double var = 0;
            for (int j = 0; j < i; j++)
                var += (coeff[i][j] * x[j]);
            for (int j = i + 1; j < systemSize; j++)
                var += (coeff[i][j] * p[j]);
            x[i] = (freeTerm[i] - var) / coeff[i][i];
        }
        iterations++;
    }
    while (!converge(x, p));

    cout << "Iterations:    " << iterations << endl;

    for (int i = 0; i<systemSize; i++)
        cout << "x[" << i << "] = " << x[i] << endl;

    delete [] p;
    delete [] x;
}

void ZeidelAndJakoby::applyJakobyMethod()
{
    cout << "--------JAKOBY--------" << endl;

    double* TempX = new double[systemSize];
    double norm;
    double* x = new double[systemSize];

    for (int j = 0; j<systemSize; j++)
        x[j] = 1;

    int iterations = 0;
    do {
        for (int i = 0; i < systemSize; i++) {
            TempX[i] = freeTerm[i];
            for (int g = 0; g < systemSize; g++) {
                if (i != g)
                    TempX[i] -= coeff[i][g] * x[g];
            }
            TempX[i] /= coeff[i][i];
        }

        norm = fabs(x[0] - TempX[0]);
        for (int h = 0; h < systemSize; h++) {
            if (fabs(x[h] - TempX[h]) > norm)
                norm = fabs(x[h] - TempX[h]);
            x[h] = TempX[h];
        }
        iterations++;
    } while (norm > eps);

    cout << "Iterations:    " << iterations << endl;

    for (int i = 0; i<systemSize; i++)
        cout << "x[" << i << "] = " << x[i] << endl;

    delete[] TempX;
}

int main()
{
    ZeidelAndJakoby solve;

    solve.initializeFromFile();
    solve.showSystem();
//    solve.change();
    solve.applyJakobyMethod();
    solve.applyZeidelMethod();

    return 0;
}
