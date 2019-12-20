import numpy as np
import matplotlib.pyplot as plt

def show_chart(*coords) :
    x,y, xnew, ynew = coords

    figure, axis = plt.subplots(figsize=(10, 10))

    plt.ion()
    axis.plot(x,y, label='Pit')

    axis.set_ylabel('$Y axis (cm)$')
    axis.set_xlabel('$X axis (cm)$')
    axis.set_xlim(left=0,right=20)
    axis.set_ylim(bottom=0,top=25)
    axis.grid(color='orange')
    figure.tight_layout()

    axis.plot(xnew,ynew,'r',label="Approximation")

    axis.legend()

    plt.draw()
    plt.pause(0.1)

    while 1 :
        tmp = input("Enter X:   ")
        curX = float(tmp)
        print(f"Pit Y   :  {calculateY(x,y,curX)}")
        print(f"Approximation Y    :  {calculateY(xnew,ynew,curX)}")

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

def buildSystem(x,y,m,n) :
    a, b = [], []

    for i in range(0,m+1) :
        row = []
        for j in range(0,m+1) :
            if (i == 0 and j == 0) :
                row.append(n)
                continue

            xsum = 0
            for k in range(0,n) :
                xsum += x[k]**(j+i)
            row.append(xsum)

        a.append(row)

        ysum = 0
        for c in range(0,n) :
            ysum += y[c] * x[c] ** i
        b.append(ysum)

    aTmp = np.array(a)
    bTmp = np.array(b)
    return aTmp, bTmp

def calculatePolynomial(a, curX) :
    value = 0
    for i in range(0,m+1):
        value += a[i] * curX**i

    return value

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

    n = 95; m = 20
    xTmp, yTmp = buildSystem(x,y,m,n)
    a = np.linalg.solve(xTmp,yTmp)

    xnew = [i*0.01 for i in range(30,1911)]
    ynew = [calculatePolynomial(a, t) for t in xnew]

    show_chart(x,y, xnew, ynew)
