import numpy as np
import math
from sys import argv

def func(point) :
    return ( 1/( math.exp(point) + math.exp(-point) ))

def simpson(*args) :
    n,a,b = args
    h = (b - a)/n
    f0 = func(a)
    fn = func(b)
    firstSum = 0
    secondSum = 0
    m = int(n/2)

    xTmp = []
    for i in range(0, n+1) :
        xTmp.append(a + h*i)

    for i in range(1,m) :
        firstSum += func(xTmp[2*i])

    for i in range(1,m+1) :
        secondSum += func(xTmp[2*i-1])

    return (h/3)*(f0 + 2*firstSum + 4*secondSum + fn)

def trapeze(*args) :
    n,a,b = args
    h = (b - a)/n
    f0 = func(a)
    fn = func(b)
    sum = 0

    xTmp = []
    for i in range(0, n+1) :
        xTmp.append(a + h*i)

    for i in range(1,n) :
        sum += func(xTmp[i])

    return (h/2)*(f0 + 2*sum + fn)


def leftRect(*args) :
    n,a,b = args
    h = (b - a)/n
    sum = 0

    xTmp = []
    for i in range(0, n+1) :
        xTmp.append(a + h*i)

    for i in range(0,n) :
        sum += func(xTmp[i])

    return h * sum

def rightRect(*args) :
    n,a,b = args
    h = (b - a)/n
    sum = 0

    xTmp = []
    for i in range(0, n+1) :
        xTmp.append(a + h*i)

    for i in range(1,n+1) :
        sum += func(xTmp[i])

    return h * sum

def centralRect(*args) :
    n,a,b = args
    f0 = func(a)/2
    fn = func(b)/2
    h = (b - a)/n
    sum = 0

    xTmp = []
    for i in range(0, n+1) :
        xTmp.append(a + h*i)

    for i in range(1,n) :
        sum += func(xTmp[i])

    return h * (f0 + sum + fn)

def checkForCondition(*args) :
    type,a,b,eps = args

    n = 2
    res2 = 0
    res1 = 0
    tmp = 1
    while tmp == 1 :
        if type == "simpson" :
            res1 = simpson(n,a,b)
        elif type == "trapeze" :
            res1 = trapeze(n,a,b)
        elif type == "leftRect" :
            res1 = leftRect(n,a,b)
        elif type == "rightRect" :
            res1 = rightRect(n,a,b)
        elif type == "centralRect" :
            res1 = centralRect(n,a,b)

        diff = abs(res1 - res2)

        res2 = res1
        n = n + 2 if type == "simpson" else n + 1

        if eps*15 >= diff and type == "simpson" :
            print(f"n = {n-2}",end=" ")
            tmp = 0
        elif eps*3 >= diff and type != "simpson" :
            print(f"n = {n-1}",end=" ")
            tmp = 0

    return res1

def calculateIntegral(*args) :
    a,b,eps = args

    print(f"Simpson method: {checkForCondition('simpson',a,b,eps)}")
    print(f"Trapeze method: {checkForCondition('trapeze',a,b,eps)}")
    print(f"Left rectangular method: {checkForCondition('leftRect',a,b,eps)}")
    print(f"Right rectangular method: {checkForCondition('rightRect',a,b,eps)}")
    print(f"Central rectangular method: {checkForCondition('centralRect',a,b,eps)}")
    print("")

if __name__ == "__main__" :
    a = 0.0
    b = 1.0
    eps1 = 0.0001
    eps2 = 0.00001

    print(f"a = {a}",f"b = {b}",f"eps = {eps1}",sep="\n")
    calculateIntegral(a,b,eps1)
    print("\n")
    print(f"a = {a}",f"b = {b}",f"eps = {eps2}",sep="\n")
    calculateIntegral(a,b,eps2)
