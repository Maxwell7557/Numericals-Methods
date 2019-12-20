#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <utility>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

class GaussMethod{
private:
    int systemSize;
    double** coeff;
    double** tmpMatrix;
    double* freeTerm;
    double* solutions;
    string filename;
    bool isSolutionFounded;
public:
    GaussMethod();
    ~GaussMethod();
    long int getRank(double**, int, int);
    void swap(double**, int, int, int);
    void initializeFromFile();
    void initializeFromUserInput();
    void showSystem();
    void applyGausMethod();
    void showSolutions();
    void deleteTmp();
    double** getCopy(bool);
};

GaussMethod::GaussMethod()
{
    filename = "matrix.txt";
    isSolutionFounded = false;
}

GaussMethod::~GaussMethod()
{
    delete [] solutions;
    for (int i = 0; i < systemSize; i++)
        delete [] coeff[i];
    delete [] coeff;
    delete [] freeTerm;
}

long int GaussMethod::getRank(double** a, int nSize, int mSize)
{
    MatrixXd matrix(nSize, mSize);
    for(int i = 0; i < nSize; i++)
    {
        for(int j = 0; j < mSize; j++)
        {
            matrix(i,j) = a[i][j];
        }
    }

    FullPivLU<MatrixXd> lu_decomp(matrix);
    long int rank = lu_decomp.rank();

    return rank;
}

void GaussMethod::swap(double** mat, int row1, int row2, int col)
{
    for (int i = 0; i < col; i++)
        {
            double temp = mat[row1][i];
            mat[row1][i] = mat[row2][i];
            mat[row2][i] = temp;
        }
}

void GaussMethod::initializeFromFile()
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

void GaussMethod::initializeFromUserInput()
{
    cout << "Р’РІРµРґРёС‚Рµ РєРѕР»РёС‡РµСЃС‚РІРѕ СѓСЂР°РІРЅРµРЅРёР№: ";
    cin >> systemSize;
    coeff = new double*[systemSize];
    freeTerm = new double[systemSize];

    cout << "Р’РІРµРґРёС‚Рµ РєРѕСЌС„С„РёС†РёРµРЅС‚С‹: \n";
    for (int i = 0; i < systemSize; i++)
    {
        cout << i+1 << " СѓСЂР°РІРЅРµРЅРёРµ: \n";
        coeff[i] = new double[systemSize];
        for (int j = 0; j < systemSize; j++)
        {
            cout << "a[" << i << "][" << j << "]= ";
            cin >> coeff[i][j];
        }

        cout << "Р’РІРµРґРёС‚Рµ СЃРІРѕР±РѕРґРЅС‹Р№ С‡Р»РµРЅ: \n";
        cout << "y[" << i << "]= ";
        cin >> freeTerm[i];
    }

}

void GaussMethod::showSystem()
{
    cout << "Р’РІРµРґРµРЅРЅР°СЏ СЃРёСЃС‚РµРјР°: \n";
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

void GaussMethod::applyGausMethod()
{
    solutions = new double[systemSize];

    long int rankOfUsualMatrix = getRank(getCopy(0),systemSize, systemSize);
    deleteTmp();
//    cout << rankOfUsualMatrix << endl;

    double x[systemSize];

    for (int i=0; i<systemSize; i++)
        coeff[i][systemSize] = freeTerm[i];

    long int rankOfExtendedMatrix = getRank(getCopy(1),systemSize,systemSize+1);
    deleteTmp();
//    cout << rankOfExtendedMatrix << endl;

    if (rankOfUsualMatrix < rankOfExtendedMatrix)
    {
        cout << "РќРµС‚ СЂРµС€РµРЅРёСЏ!" << endl;
        return;
    }
    else if (rankOfUsualMatrix == rankOfExtendedMatrix && rankOfUsualMatrix < systemSize)
    {
        cout << "РЎРёСЃС‚РµРјР° РёРјРµРµС‚ Р±РµСЃРєРѕРЅРµС‡РЅРѕРµ РјРЅРѕР¶РµСЃС‚РІРѕ СЂРµС€РµРЅРёР№" << endl;
        return;
    }
    else if (rankOfUsualMatrix == rankOfExtendedMatrix && rankOfUsualMatrix == systemSize)
        cout << "РЎРёСЃС‚РµРјР° РёРјРµРµС‚ РµРґРёРЅСЃС‚РІРµРЅРЅРѕРµ СЂРµС€РµРЅРёРµ!" << endl;

    for(int i=systemSize-1; i>0; i--)
    {
        if(coeff[i-1][0]<coeff[i][0])
        {
            for(int j=0; j<=systemSize; j++)
            {
                double c = coeff[i][j];
                coeff[i][j] = coeff[i-1][j];
                coeff[i-1][j] = c;
            }
        }
    }

//    for (int i = 0; i< systemSize; i++)
//        for (int j = 0; j< systemSize; j++)
//            cout << coeff[i][j];

    for(int k=0; k<systemSize-1; k++)
        for(int i=k; i<systemSize-1; i++)
        {
            double c = (coeff[i+1][k]/coeff[k][k]) ;

            for(int j=0; j<=systemSize; j++)
                coeff[i+1][j]-=c*coeff[k][j];
        }

//    cout << endl;
//    for (int i = 0; i< systemSize; i++)
//        for (int j = 0; j< systemSize; j++)
//            cout << coeff[i][j];

    for(int i=systemSize-1; i>=0; i--)
    {
        double c = 0;
        for(int j=i; j<=systemSize-1; j++)
            c = c+coeff[i][j]*x[j];

        x[i]=(coeff[i][systemSize]-c)/coeff[i][i];
//        cout << x[i] << endl;
    }

    isSolutionFounded = true;

    for (int i=0; i<systemSize; i++)
        solutions[i] = x[i];

//    delete [] x;
}

void GaussMethod::showSolutions()
{
    if (isSolutionFounded)
    {
        cout << "Р РµС€РµРЅРёРµ РІРІРµРґРµРЅРЅРѕР№ СЃРёСЃС‚РµРјС‹: \n";
        for (int i = 0; i < systemSize; i++)
            cout << "x[" << i << "]=" << solutions[i] << endl;
    }
}

void GaussMethod::deleteTmp()
{
    for (int i = 0; i < systemSize; i++)
        delete [] tmpMatrix[i];
    delete [] tmpMatrix;
}

double **GaussMethod::getCopy(bool mode)
{
    tmpMatrix = new double*[systemSize];
    int nSize = 0;
    int mSize = systemSize;
    mode ? nSize = systemSize + 1 : nSize = systemSize;

    for (int i = 0; i < mSize; i++)
    {
        tmpMatrix[i] = new double[nSize];
        for (int j = 0; j < nSize; j++)
            tmpMatrix[i][j] = coeff[i][j];
    }

    return tmpMatrix;
}

int main(int arg, char** argv)
{
    GaussMethod gaussMethod;

    if(strcmp(argv[1],"user") == 0)
        gaussMethod.initializeFromUserInput();
    else if (strcmp(argv[1],"file") == 0)
        gaussMethod.initializeFromFile();

    gaussMethod.showSystem();
    gaussMethod.applyGausMethod();
    gaussMethod.showSolutions();

    return 0;
}
