import numpy as np
import matplotlib.pyplot as plt
from sys import argv

def show_chart(*coords) :
    x, y, lx, ly, xnew, ynew, lxnew, lynew = coords

    figure, axis = plt.subplots(figsize=(10, 10))

    # plt.ion()
    axis.plot(x,y, label='Яма')

    axis.set_ylabel('$Y axis (cm)$')
    axis.set_xlabel('$X axis (cm)$')
    axis.set_xlim(left=0,right=20)
    axis.set_ylim(bottom=0,top=25)
    axis.minorticks_on()
    axis.grid(which = "major",color='orange',linewidth = 1)
    axis.grid(which = "minor",color='orange',linestyle = ":")
    figure.tight_layout()

    axis.plot(xnew,ynew,'r',label='Лагранж')
    axis.plot(lxnew,lynew,'g--',linewidth=2,label='Ньютон')

    axis.legend()

    plt.draw()
    plt.pause(0.1)

    while 1 :
        # tmp = input()
        curX = float(input("Enter X:    "))
        print(f"Pit Y:  {calculateY(x,y,curX)}")
        print(f"Lagrange Y:   {calculateY(xnew,ynew,curX)}")
        print(f"Newton Y:   {calculateY(lxnew,lynew,curX)}")
        print("")

        plt.plot(curX,calculateY(x,y,curX),'o')
        plt.plot(curX,calculateY(xnew,ynew,curX),'o')
        plt.plot(curX,calculateY(lxnew,lynew,curX),'o')
        plt.draw()

def L_polynomial(*polinomialArgs) :
    x, y, curX = polinomialArgs

    pn=0
    for k in range(len(y)):
        top=1; bottom=1
        for i in range(len(x)):
            if i == k:
                top = top*1; bottom = bottom*1
            else:
                top = top*(curX-x[i])
                bottom = bottom*(x[k]-x[i])
        pn = pn + y[k]*top/bottom
    return pn

def newtonCoefficient(x,y) :
    m = len(x)
    x = np.copy(x)
    a = np.copy(y)

    for k in range(1,m):
        a[k:m] = (a[k:m] - a[k-1])/(x[k:m] - x[k-1])

    return a

def calculateNewton(x, y, curX) :
    a = newtonCoefficient(x, y)
    n = len(x) - 1
    p = a[n]

    for k in range(1,n+1) :
        p = a[n-k] + (curX -x[n-k])*p

    return p

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
    amountOfPoints = int(amountOfPoints)

    x = np.array([],dtype=float)
    y = np.array([],dtype=float)

    lx = []
    ly = []

    if (amountOfPoints == 3) :
        lx.append(0.3); ly.append(22.7)
        lx.append(10.1); ly.append(3.3)
        lx.append(19.1); ly.append(22.8)
    elif (amountOfPoints == 5) :
        lx.append(2.1);  ly.append(21.6)
        lx.append(6.1);  ly.append(13.4)
        lx.append(10.1); ly.append(3.3)
        lx.append(14.1); ly.append(2.7)
        lx.append(18.1); ly.append(20.7)
    elif (amountOfPoints == 10) :
        lx.append(0.3); ly.append(22.7)
        lx.append(2.3);  ly.append(21.5)
        lx.append(4.3);  ly.append(18.9)
        lx.append(6.3);  ly.append(13.1)
        lx.append(8.3);  ly.append(10.7)
        lx.append(10.3); ly.append(2.8)
        lx.append(12.3); ly.append(1.0)
        lx.append(14.3); ly.append(3.0)
        lx.append(16.3); ly.append(9.4)
        lx.append(19.1); ly.append(22.8)

    file = open("coord.csv")

    for line in file :
        sep = line.find(';')

        if line[0].isalpha() :
            continue

        x = np.append(x,float(line[0:sep]))
        y = np.append(y,float(line[sep+1:]))


    xnew = [i*0.1 for i in range(3,192)]
    ynew = [L_polynomial(lx,ly,t) for t in xnew]

    lxnew = xnew
    lynew = [calculateNewton(lx,ly,t) for t in lxnew]

    file.close()

    show_chart(x,y, lx, ly, xnew,ynew, lxnew,lynew)
