import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import interpolate

def show_chart(*coords) :
    x,y, xnew, ynew = coords

    # figure = plt.figure(figsize=(10, 10))
    # axis = plt.subplots(2,1,1)
    # ax = plt.subplots(2,1,2)
    figure, axis = plt.subplots(figsize=(10, 10))

    # plt.ion()
    axis.plot(x,y, label='Pit')

    axis.set_ylabel('$Y axis (cm)$')
    axis.set_xlim(left=0,right=20)
    axis.set_xlabel('$X axis (cm)$')
    axis.set_ylim(bottom=0,top=25)
    axis.grid(color='orange')
    axis.legend()

    figure.tight_layout()



    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_ylabel('$Y axis (cm)$')
    ax.set_xlabel('$X axis (cm)$')
    # ax.set_xlim(left=0,right=20)
    # ax.set_ylim(bottom=-10,top=10)
    ax.grid(color='orange')
    fig.tight_layout()

    ax.plot(xnew,ynew,'r',label="Differentiation")
    # ax.plot(xxnew,yynew,'b',label="Differentiation x")
    ax.legend()


    plt.show()
    # plt.pause(0.1)

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

def calculateDiff(x, y) :
    derivative = []
    first = (-3*y[0] + 4*y[1] - y[2])/(2 * 0.01)
    derivative.append(first);

    for i in range(1,len(x)-1) :
        value = (y[i+1] - y[i-1])/(2 * 0.01)
        derivative.append(value)

    last = (y[len(y)-1-2] - 4*y[len(y)-1-1] + 3*y[len(y)-1])/(2 * 0.01)
    derivative.append(last);

    return derivative

def spline(x,y,curX) :
    tck = interpolate.splrep(x,y)
    return interpolate.splev(curX,tck)

if __name__ == "__main__" :
    x = np.array([],dtype=float)
    y = np.array([],dtype=float)

    file = open("spline_coords.csv")

    for line in file :
        sep = line.find(';')

        if line[0].isalpha() :
            continue

        x = np.append(x,float(line[0:sep]))
        y = np.append(y,float(line[sep+1:]))

    file.close()

    # n = 95
    xSpl = [i*0.01 for i in range(30,1911)]
    # xSpl = [i*0.1 for i in range(3,191)]
    ySpl = [spline(x, y, t) for t in xSpl]

    xTmp, yTmp = xSpl, calculateDiff(xSpl,ySpl)
    # xTmp = np.delete(xTmp,[0,len(xTmp)-1])

    # xTmp2, yTmp2 = xSpl, calculateSecondDiff(xSpl,ySpl)
    # xTmp2, yTmp2 = xSpl, calculateDiffX(xSpl,ySpl)
    # xTmp2, yTmp2 = xSpl, calculateDiff(xTmp,yTmp)
    # xTmp2 = np.delete(xTmp2,[0,len(xTmp2)-1])

    show_chart(xSpl,ySpl, xTmp, yTmp)
