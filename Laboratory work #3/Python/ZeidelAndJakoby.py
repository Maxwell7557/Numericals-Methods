import numpy as np
import math

class ZeidelAndJakoby :
    def __init__ (self) :
        self.systemSize = 0
        self.coeff = []
        self.freeTerm = []
        self.filename = "mat.txt"
        self.eps = 0.001

    def initializationFromFile(self) :
        with open(self.filename) as file :
            self.coeff = [list(map(float, row.split())) for row in file.readlines()]
            self.systemSize = len(self.coeff[0])-1

        for i in range(0,self.systemSize) :
            self.freeTerm.append([self.coeff[i][-1]])
            self.coeff[i].pop()

    def showSystem(self) :
        print("Р’РІРµРґРµРЅРЅР°СЏ СЃРёСЃС‚РµРјР°: ")
        for row in range(len(self.freeTerm)):
            for col in range(len(self.coeff[row])):
                print(f"{self.coeff[row][col]}*x{col}", end='')
                if col < self.systemSize - 1 :
                    print(" + ", end='')
            print(f" = {self.freeTerm[row][0]}")


    def applyJakobyMethod(self) :
        print("--------JAKOBY--------")
        npCoeff = np.asarray(self.coeff)
        npFreeTerm = np.asarray(self.freeTerm)

        tmp = []
        for i in range(self.systemSize) :
            tmp.append([])
            for j in range(self.systemSize) :
                tmp[i].append(0)
            tmp[i][i] = self.coeff[i][i]

        diagonalMatrix = np.asarray(tmp)
        inversedDiagMatrix = np.linalg.inv(diagonalMatrix)

        # print(npCoeff)
        # print(npFreeTerm)
        # print(inversedDiagMatrix)

        b = np.eye(self.systemSize) - np.dot(inversedDiagMatrix,npCoeff)
        d = np.dot(inversedDiagMatrix,npFreeTerm)

        currentX = np.array([[0]]*self.systemSize,dtype = float)
        previousX = np.array([[1]]*self.systemSize,dtype = float)

        iterations = 0
        while True :
            currentX = np.dot(b,previousX) + d

            diff = currentX - previousX

            maxValue = math.fabs(diff.max())
            minValue = math.fabs(diff.min())

            iterations += 1
            if maxValue > minValue and maxValue < self.eps :
                print(f"Number of iterations: {iterations}")
                for i in range(self.systemSize) :
                    print(f"x[{i}] = {currentX[i]}")
                break;
            elif maxValue < minValue and minValue < self.eps :
                print(f"Number of iterations: {iterations}")
                for i in range(self.systemSize) :
                    print(f"x[{i}] = {currentX[i]}")
                break;

            previousX = currentX

    def applyZeidelMethod(self) :
        print("--------ZEIDEL--------")
        npCoeff = np.asarray(self.coeff)
        npFreeTerm = np.asarray(self.freeTerm)

        triangMatrix = np.tril(npCoeff)
        inversedMatrix = np.linalg.inv(triangMatrix)

        b = np.eye(self.systemSize) - np.dot(inversedMatrix,npCoeff)
        d = np.dot(inversedMatrix,npFreeTerm)

        currentX = np.array([[0]]*self.systemSize,dtype = float)
        previousX = np.array([[1]]*self.systemSize,dtype = float)

        iterations = 0
        while True :
            currentX = np.dot(b,previousX) + d

            diff = currentX - previousX

            maxValue = math.fabs(diff.max())
            minValue = math.fabs(diff.min())

            iterations += 1
            if maxValue > minValue and maxValue < self.eps :
                print(f"Number of iterations: {iterations}")
                for i in range(self.systemSize) :
                    print(f"x[{i}] = {currentX[i]}")
                break;
            elif maxValue < minValue and minValue < self.eps :
                print(f"Number of iterations: {iterations}")
                for i in range(self.systemSize) :
                    print(f"x[{i}] = {currentX[i]}")
                break;

            previousX = currentX

if __name__ == '__main__':
    solve = ZeidelAndJakoby()

    solve.initializationFromFile()
    solve.showSystem()
    solve.applyJakobyMethod()
    solve.applyZeidelMethod(
