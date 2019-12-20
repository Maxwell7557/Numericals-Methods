import numpy as np
import matplotlib.pyplot as plt
from sys import argv

class SplineCoefs :
    a = 0.0
    b = 0.0
    c = 0.0
    d = 0.0
    x = 0.0

def calculateCoefs(*args) :
    x, y, n, splineCoefs = args

    for i in range(0,n) :
        splineCoefs[i].x = x[i]
        splineCoefs[i].a = y[i]

    splineCoefs[0].c = 0

    alpha = [0.0 for i in range(0,n)]
    beta = [0.0 for i in range(0,n)]
    A = 0.0
    B = 0.0
    C = 0.0
    F = 0.0
    h_i = 0.0
    h_i1 = 0.0
    z = 0.0
    alpha[0] = 0.0
    beta[0] = 0.0
    for i in range(1,n-1) :
        h_i = x[i] - x[i - 1]
        h_i1 = x[i + 1] - x[i]
        A = h_i
        C = 2.0 * (h_i + h_i1)
        B = h_i1
        F = 6.0 * ((y[i + 1] - y[i]) / h_i1 - (y[i] - y[i - 1]) / h_i)
        z = (A * alpha[i - 1] + C)
        alpha[i] = -B / z
        beta[i] = (F - A * beta[i - 1]) / z

    splineCoefs[n-1].c = (F - A * beta[n - 2]) / (C + A * alpha[n - 2])

    for i in range(n-2,0,-1) :
        splineCoefs[i].c = alpha[i] * splineCoefs[i + 1].c + beta[i]

    for i in range(n-1,0,-1) :
        h_i = x[i] - x[i-1]
        splineCoefs[i].d = (splineCoefs[i].c - splineCoefs[i - 1].c) / h_i
        splineCoefs[i].b = h_i * (2.0 * splineCoefs[i].c + splineCoefs[i - 1].c) / 6.0 + (y[i] - y[i - 1]) / h_i

def calculateSpline (*args) :
    curX, splineCoefs = args
    n = len(splineCoefs)
    # print(len(splineCoefs))

    splineCoefsTmp = SplineCoefs()
    if curX <= splineCoefs[0].x :
        splineCoefsTmp = splineCoefs[1]
    elif curX >= splineCoefs[n-1].x :
        splineCoefsTmp = splineCoefs[n - 1]
    else :
        for j in range(0,n) :
            if curX <= splineCoefs[j].x :
                splineCoefsTmp = splineCoefs[j]
                break

    dx = (curX - splineCoefsTmp.x)
    return splineCoefsTmp.a + (splineCoefsTmp.b + (splineCoefsTmp.c / 2.0 + splineCoefsTmp.d * dx / 6.0) * dx) * dx

def show_chart(*coords) :
    x,y, xnew, ynew, lb, splineCoefs = coords

    figure, axis = plt.subplots(figsize=(10, 10))

    plt.ion()
    axis.plot(x,y,'b',label='pit')
    axis.set_ylabel('$Y axis (cm)$')
    axis.set_xlabel('$X axis (cm)$')
    axis.set_xlim(left=0,right=20)
    axis.set_ylim(bottom=0,top=25)
    axis.grid(color='orange')
    figure.tight_layout()

    axis.plot(xnew,ynew,'r',label=lb)

    axis.legend()

    # plt.show()
    plt.draw()
    plt.pause(0.1)

    while 1 :
        a = float(input('Enter left border in cm: '))
        b = float(input('Enter right border in cm: '))
        eps = float(input('Enter accuracy: '))
        print("")

        first = xnew.index(a)
        second = xnew.index(b)
        amountOfPoints = int((b - a)/0.01)

        calculateIntegral(xnew,amountOfPoints,a,b,first,splineCoefs,eps)

        axis.plot(a,ynew[first],'o')
        axis.plot(b,ynew[second],'o')
    plt.draw()

def simpson(*args) :
    xnew,n,a,b,first,splineCoefs = args
    h = (b - a)/n
    f0 = calculateSpline(a, splineCoefs)
    fn = calculateSpline(b, splineCoefs)
    firstSum = 0
    secondSum = 0
    m = int(n/2)

    xTmp = []
    for i in range(0, n+1) :
        xTmp.append(a + h*i)

    for i in range(1,m) :
        firstSum += calculateSpline(xTmp[2*i], splineCoefs)

    for i in range(1,m+1) :
        secondSum += calculateSpline(xTmp[2*i-1], splineCoefs)

    return (h/3)*(f0 + 2*firstSum + 4*secondSum + fn)

def trapeze(*args) :
    xnew,n,a,b,first,splineCoefs = args
    h = (b - a)/n
    f0 = calculateSpline(a, splineCoefs)
    fn = calculateSpline(b, splineCoefs)
    sum = 0

    xTmp = []
    for i in range(0, n+1) :
        xTmp.append(a + h*i)

    for i in range(1,n) :
        sum += calculateSpline(xTmp[i], splineCoefs)

    return (h/2)*(f0 + 2*sum + fn)


def leftRect(*args) :
    xnew,n,a,b,first,splineCoefs = args
    h = (b - a)/n
    sum = 0

    xTmp = []
    for i in range(0, n+1) :
        xTmp.append(a + h*i)

    for i in range(0,n) :
        sum += calculateSpline(xTmp[i], splineCoefs)

    return h * sum

def rightRect(*args) :
    xnew,n,a,b,first,splineCoefs = args
    h = (b - a)/n
    sum = 0

    xTmp = []
    for i in range(0, n+1) :
        xTmp.append(a + h*i)

    for i in range(1,n+1) :
        sum += calculateSpline(xTmp[i], splineCoefs)

    return h * sum

def centralRect(*args) :
    xnew,n,a,b,first,splineCoefs = args
    f0 = calculateSpline(a, splineCoefs)/2
    fn = calculateSpline(b, splineCoefs)/2
    h = (b - a)/n
    sum = 0

    xTmp = []
    for i in range(0, n+1) :
        xTmp.append(a + h*i)

    for i in range(1,n) :
        sum += calculateSpline(xTmp[i], splineCoefs)

    return h * (f0 + sum + fn)

def checkForCondition(*args) :
    type,xnew,step,a,b,first,splineCoefs,eps = args

    n = 10
    res2 = 0
    res1 = 0
    while True :
        if type == "simpson" :
            res1 = simpson(xnew,n,a,b,first,splineCoefs)
        elif type == "trapeze" :
            res1 = trapeze(xnew,n,a,b,first,splineCoefs)
        elif type == "leftRect" :
            res1 = leftRect(xnew,n,a,b,first,splineCoefs)
        elif type == "rightRect" :
            res1 = rightRect(xnew,n,a,b,first,splineCoefs)
        elif type == "centralRect" :
            res1 = centralRect(xnew,n,a,b,first,splineCoefs)

        diff = abs(res1 - res2)

        res2 = res1
        n = n + 2 if type == "simpson" else n + 1

        if eps*15 >= diff and type == "simpson" :
            print(f"n = {n-2}",end=" ")
            break
        elif eps*3 >= diff and type != "simpson" :
            print(f"n = {n-1}",end=" ")
            break

    return res1

def calculateIntegral(*args) :
    xnew,step,a,b,first,splineCoefs,eps = args

    print(f"Simpson method: {checkForCondition('simpson',xnew,step,a,b,first,splineCoefs,eps)}")
    print(f"Trapeze method: {checkForCondition('trapeze',xnew,step,a,b,first,splineCoefs,eps)}")
    print(f"Left rectangular method: {checkForCondition('leftRect',xnew,step,a,b,first,splineCoefs,eps)}")
    print(f"Right rectangular method: {checkForCondition('rightRect',xnew,step,a,b,first,splineCoefs,eps)}")
    print(f"Central rectangular method: {checkForCondition('centralRect',xnew,step,a,b,first,splineCoefs,eps)}")
    print("")

if __name__ == "__main__" :
    x = np.array([],dtype=float)
    y = np.array([],dtype=float)

    file = open("coord.csv")

    for line in file :
        sep = line.find(';')

        if line[0].isalpha() :
            continue

        x = np.append(x,float(line[0:sep]))
        y = np.append(y,float(line[sep+1:]))

    file.close()

    n1 = 95
    splineCoefs = []
    for i in range(0,95) :
        splineCoefs.append(SplineCoefs())

    calculateCoefs(x,y,n1,splineCoefs)
    xnew = [round(i*0.001,2) for i in range(300,19101)]
    ynew = [calculateSpline(t, splineCoefs) for t in xnew]
    label='Spline on all points'

    show_chart(x,y, xnew, ynew, label, splineCoefs)
