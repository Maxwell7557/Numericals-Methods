import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from scipy import interpolate

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
    curX, n, splineCoefs = args

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
    x,y, xnew, ynew, lb = coords

    figure, axis = plt.subplots(figsize=(10, 10))

    plt.ion()
    axis.plot(x,y, label='Pit')

    axis.set_ylabel('$Y axis (cm)$')
    axis.set_xlabel('$X axis (cm)$')
    axis.set_xlim(left=0,right=20)
    axis.set_ylim(bottom=0,top=25)
    axis.grid(color='orange')
    figure.tight_layout()

    axis.plot(xnew,ynew,'r',label=lb)

    axis.legend()

    plt.draw()
    plt.pause(0.1)

    while 1 :
        tmp = input("Enter X:   ")
        curX = float(tmp)
        print(f"Pit Y   :  {calculateY(x,y,curX)}")
        print(f"Spline Y    :  {calculateY(xnew,ynew,curX)}")

        axis.plot(curX,calculateY(x,y,curX),'o')
        axis.plot(curX,calculateY(xnew,ynew,curX),'o')
    plt.draw()

def calculateY(x, y, curX) :
    x1 = 0; x2 = 0
    y1 = 0; y2 = 0

    for i in range(len(x)) :
        if x[i] == curX :
            return y[i];
        if x[i] < curX and x[i+1] > curX :
            x1 = x[i]; x2 = x[i+1]
            y1 = y[i]; y2 = y[i+1]
            break;

    top1 = curX*y2 - curX*y1
    top2 = x2*y1-x1*y2
    bottom = x2 - x1
    y = (top1 + top2) / bottom
    return y

if __name__ == "__main__" :
    scriptName, amountOfPoints = argv

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

    xnew = []
    ynew = []
    n1 = 95
    n2 = int(95/3)
    n3 = int(95/5)
    splineCoefs = []
    for i in range(0,95) :
        splineCoefs.append(SplineCoefs())
    label=""

    if amountOfPoints == 'all' :
        calculateCoefs(x,y,n1,splineCoefs)
        xnew = [i*0.05 for i in range(6,383)]
        ynew = [calculateSpline(t, n1, splineCoefs) for t in xnew]
        # ynew = [calc(t,x,y) for t in xnew]
        label='Spline on all points'
    elif amountOfPoints == '3' :
        xTmp = x[::3]; xTmp = np.append(xTmp,x[-1])
        yTmp = y[::3]; yTmp = np.append(yTmp,y[-1])
        n2 = n2 + 2
        # print(xTmp)
        calculateCoefs(xTmp,yTmp,n2,splineCoefs)
        # xnew = [i*0.05 for i in range(6,383)]
        xnew = [i*0.01 for i in range(30,1911)]
        ynew = [calculateSpline(t, n2, splineCoefs) for t in xnew]
        label='Spline on each 3 point'
    elif amountOfPoints == '5' :
        xTmp = x[::5]; xTmp = np.append(xTmp,x[-1])
        yTmp = y[::5]; yTmp = np.append(yTmp,y[-1])
        n3 = n3 + 1
        calculateCoefs(xTmp,yTmp,n3,splineCoefs)
        xnew = [i*0.05 for i in range(6,383)]
        ynew = [calculateSpline(t, n3, splineCoefs) for t in xnew]
        label='Spline on each 5 point'

    show_chart(x,y, xnew, ynew, label)
