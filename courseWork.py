import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import interpolate

def show_chart(*coords) :
    x,y, = coords

    figure, axis = plt.subplots(2,1,figsize=(10, 10))

    plt.ion()
    axis[0].plot(x,y, label='Pit')

    axis[0].set_ylabel('$Y axis (cm)$')
    axis[0].set_xlim(left=0,right=20)
    axis[0].set_xlabel('$X axis (cm)$')
    axis[0].set_ylim(bottom=0,top=25)
    axis[0].grid(color='orange')
    axis[0].legend()

    axis[1].set_ylabel('$V$')
    axis[1].set_xlabel('$X$')
    axis[1].grid(color='orange')

    figure.tight_layout()

    plt.draw()
    plt.pause(0.1)

    x0 = float(input("Enter initial X:    "))
    v0 = float(input("Enter initial V:    "))
    curX = x0
    curV = v0
    time = 0
    v = []
    coords = []
    a,b = 0.1,0.25
    c = np.linspace(0,2*math.pi,100)
    while 1 :
        # circle = plt.Circle((curX, spline(x,y,curX)+0.2), 0.2, color='red')
        circle = axis[0].plot(curX + a*np.cos(c),spline(x,y,curX) + b*np.sin(c) + b,color="red")
        # circle = axis[0].plot(curX,spline(x,y,curX),'bo')
        coords.append(curX); v.append(curV)
        print(f"Current x: {curX}", end="   ")
        print(f"Current v: {curV}")
        line = axis[1].plot(coords,v,'b')

        # axis[0].add_artist(circle)
        time += 0.5
        tmp = currentPosition(x,y,curX,curV,time)
        curX, curV = tmp

        plt.pause(0.001)
        # circle.remove()
        axis[0].lines[1].remove()
        axis[1].lines[0].remove()

        if curX < x[0] or curX > x[-1] :
            print("ERROR: The ball passed the boundary")
            break

    esc = input("Enter something to exit: ")
    plt.draw()

def calculateDiffInPoint(x,y,curX) :
    value = ( spline(x,y,curX+0.01) - spline(x,y,curX-0.01) )/(2*0.01)
    return value

def currentPosition(xArr,yArr,curX,curZ,t) :
    dt = 0.1

    k1 = dt* f(t,curX,curZ)
    m1 = dt* g(t,curX,curZ,xArr,yArr)

    k2 = dt* f(t+dt,curX+k1/2,curZ+m1/2)
    m2 = dt* g(t+dt,curX+k1/2,curZ+m1/2,xArr,yArr)

    k3 = dt* f(t+dt,curX+k2/2,curZ+m2/2)
    m3 = dt* g(t+dt,curX+k2/2,curZ+m2/2,xArr,yArr)

    k4 = dt* f(t+dt,curX+k3,curZ+m3)
    m4 = dt* g(t+dt,curX+k3,curZ+m3,xArr,yArr)

    deltaX = (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    deltaZ = (1/6)*(m1 + 2*m2 + 2*m3 + m4)

    newX = curX + deltaX
    newZ = curZ + deltaZ

    return newX, newZ

def g(t,x,z,xArr,yArr) :
    g = 9.8
    coeffOfFriction = 0.2
    m = 1

    currentValueOfSpline = calculateDiffInPoint(xArr,yArr,x)
    alpha = -g
    # beta = -g*coeffOfFriction
    beta = -(coeffOfFriction/m)

    # value = alpha*currentValueOfSpline + (beta*currentValueOfSpline)/((1+currentValueOfSpline**2)**(1/2))
    value = alpha* currentValueOfSpline + beta*z
    return value

def f(t,x,z) :
    value = z
    return value

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

    xSpl = [i*0.01 for i in range(30,1911)]
    ySpl = [spline(x, y, t) for t in xSpl]

    show_chart(xSpl,ySpl)
