#include <iostream>
#include <cmath>

using namespace std;

double func(double point)
{
    return ( 1/( exp(point) + exp(-point) ));
}

double simpson(int n,double a,double b)
{
    double h = (b - a)/n;
    double f0 = func(a);
    double fn = func(b);
    double firstSum = 0;
    double secondSum = 0;
    int m = int(n/2);

    double* xTmp;
    xTmp = new double[n+1];

    for (int i = 0; i < n+1; i++)
        xTmp[i] = (a + h*i);

    for (int i = 1; i < m; i++)
        firstSum += func(xTmp[2*i]);

    for (int i = 1; i < m+1; i++)
        secondSum += func(xTmp[2*i-1]);

    delete [] xTmp;
    return (h/3)*(f0 + 2*firstSum + 4*secondSum + fn);
}

double trapeze(int n,double a,double b)
{
    double h = (b - a)/n;
    double f0 = func(a);
    double fn = func(b);
    double sum = 0;

    double* xTmp;
    xTmp = new double[n+1];

    for (int i = 0; i < n+1; i++)
        xTmp[i] = (a + h*i);

    for (int i = 1; i < n; i++)
        sum += func(xTmp[i]);

    delete [] xTmp;
    return (h/2)*(f0 + 2*sum + fn);
}

double leftRect(int n,double a,double b)
{
    double h = (b - a)/n;
    double sum = 0;

    double* xTmp;
    xTmp = new double[n+1];

    for (int i = 0; i < n+1; i++)
        xTmp[i] = (a + h*i);

    for (int i = 0; i < n; i++)
        sum += func(xTmp[i]);

    delete [] xTmp;
    return h * sum;
}

double rightRect(int n,double a,double b)
{
    double h = (b - a)/n;
    double sum = 0;

    double* xTmp;
    xTmp = new double[n+1];

    for (int i = 0; i < n+1; i++)
        xTmp[i] = (a + h*i);

    for (int i = 1; i < n+1; i++)
        sum += func(xTmp[i]);

    delete [] xTmp;
    return h * sum;
}

double centralRect(int n,double a,double b)
{
    double f0 = func(a)/2;
    double fn = func(b)/2;
    double h = (b - a)/n;
    double sum = 0;

    double* xTmp;
    xTmp = new double[n+1];

    for (int i = 0; i < n+1; i++)
        xTmp[i] = (a + h*i);

    for (int i = 1; i < n; i++)
        sum += func(xTmp[i]);

    delete [] xTmp;
    return h * (f0 + sum + fn);
}

double checkForCondition(string type,double a,double b, double eps)
{
    int n = 2;
    double res2 = 0;
    double res1 = 0;
    bool tmp = true;

    while (tmp)
    {
        if (type == "simpson")
            res1 = simpson(n,a,b);
        else if (type == "trapeze")
            res1 = trapeze(n,a,b);
        else if (type == "leftRect")
            res1 = leftRect(n,a,b);
        else if (type == "rightRect")
            res1 = rightRect(n,a,b);
        else if (type == "centralRect")
            res1 = centralRect(n,a,b);

        double diff = abs(res1 - res2);

        res2 = res1;

        type == "simpson" ? n += 2 : n += 1;
//        if (type == "simpson")
//            n += 2;
//        else
//            n += 1;


        if (eps*15 >= diff && type == "simpson")
        {
            cout << "n = " << n-2 << endl;
            tmp = false;
        } else if (eps*3 >= diff && type != "simpson"){
            cout << "n = " << n-1 << endl;
            tmp = false;
        }
    }
    return res1;
}

void calculateIntegral(double a,double b, double eps)
{
    cout << "Simpson method: " << checkForCondition("simpson",a,b,eps) << endl;
    cout << "Trapeze method: " << checkForCondition("trapeze",a,b,eps) << endl;
    cout << "Left rectangular method: " << checkForCondition("leftRect",a,b,eps) << endl;
    cout << "Right rectangular method: " << checkForCondition("rightRect",a,b,eps) << endl;
    cout << "Central rectangular method: " << checkForCondition("centralRect",a,b,eps) << endl;
    cout << endl;
}

int main()
{
    double a,b,eps1,eps2;
    a = 0.0;
    b = 1.0;
    eps1 = 0.0001;
    eps2 = 0.00001;

    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "eps = " << eps1 << endl;
    calculateIntegral(a,b,eps1);
    cout << endl;
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "eps = " << eps2 << endl;
    calculateIntegral(a,b,eps2);
}
