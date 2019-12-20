import numpy as np
from sys import argv

class GaussMethod :

    def __init__ (self) :
        self.systemSize = 0
        self.coeff = []
        self.freeTerm = []
        self.solutions = []
        self.filename = "matrix.txt"
        self.isSolutionFounded = False

    def initializationFromUserInput(self):
        self.systemSize = int(input("Р’РІРµРґРёС‚Рµ РєРѕР»РёС‡РµСЃС‚РІРѕ СѓСЂР°РІРЅРµРЅРёР№: "))

        print("Р’РІРµРґРёС‚Рµ РєРѕСЌС„С„РёС†РёРµРЅС‚С‹: ")
        for i in range(self.systemSize) :
            print(f"{i+1} СѓСЂР°РІРЅРµРЅРёРµ: ")
            self.coeff.append([])

            for j in range(self.systemSize) :
                self.coeff[i].append(float(input(f"a[{i}][{j}]= ")))

            print("Р’РІРµРґРёС‚Рµ СЃРІРѕР±РѕРґРЅС‹Р№ С‡Р»РµРЅ: ")
            self.freeTerm.append(float(input(f"y[{i}]= ")))

        # print("Р’РІРµРґРёС‚Рµ СЃРІРѕР±РѕРґРЅС‹Р№ С‡Р»РµРЅ: ")
        # for i in range(self.systemSize) :
        #     self.freeTerm.append(int(input(f"y[{i}]= ")))

        self.solutions = [None for i in range(0,self.systemSize)]

    def initializationFromFile(self) :
        with open('matrix.txt') as file :
            self.coeff = [list(map(float, row.split())) for row in file.readlines()]
            self.systemSize = len(self.coeff[0])-1

        for i in range(0,self.systemSize) :
            self.freeTerm.append(self.coeff[i][-1])
            self.coeff[i].pop()

        self.solutions = [None for i in range(0,self.systemSize)]

    def showSystem(self):
        print("Р’РІРµРґРµРЅРЅР°СЏ СЃРёСЃС‚РµРјР°: ")
        for row in range(len(self.freeTerm)):
            for col in range(len(self.coeff[row])):
                print(f"{self.coeff[row][col]}*x{col}", end='')
                if col < self.systemSize - 1 :
                    print(" + ", end='')
            print(f" = {self.freeTerm[row]}")

    def applyGaussMethod (self):
        x = [0 for i in range(0,self.systemSize)]

        tmpNpArr = np.array(self.coeff)
        rankOfUsualMatrix = np.linalg.matrix_rank(tmpNpArr)
        # determinator = np.linalg.det(self.coeff)

        for i in range(0,self.systemSize) :
            self.coeff[i].append(self.freeTerm[i])

        tmpNpArr = np.array(self.coeff)
        rankOfExtendedMatrix = np.linalg.matrix_rank(tmpNpArr)

        if rankOfUsualMatrix < rankOfExtendedMatrix :
            print("РќРµС‚ СЂРµС€РµРЅРёСЏ!")
            return
        elif rankOfUsualMatrix == rankOfExtendedMatrix and rankOfUsualMatrix < self.systemSize :
            print("РЎРёСЃС‚РµРјР° РёРјРµРµС‚ Р±РµСЃРєРѕРЅРµС‡РЅРѕРµ РјРЅРѕР¶РµСЃС‚РІРѕ СЂРµС€РµРЅРёР№")
            return
        elif rankOfUsualMatrix == rankOfExtendedMatrix and rankOfUsualMatrix == self.systemSize :
            print("РЎРёСЃС‚РµРјР° РёРјРµРµС‚ РµРґРёРЅСЃС‚РІРµРЅРЅРѕРµ СЂРµС€РµРЅРёРµ!")

        for i in range(self.systemSize-1,0,-1) :
            if self.coeff[i-1][0] < self.coeff[i][0] :
                for j in range(0,self.systemSize+1) :
                    self.coeff[i][j], self.coeff[i-1][j] = self.coeff[i-1][j], self.coeff[i][j]

        for k in range(0,self.systemSize-1) :
            for i in range(k,self.systemSize-1) :
                tmp = (self.coeff[i+1][k] / self.coeff[k][k])

                for j in range(0,self.systemSize+1) :
                    self.coeff[i+1][j] -= tmp * self.coeff[k][j]

        for i in range(self.systemSize-1,-1,-1) :
            tmp = 0.0
            for j in range(i,self.systemSize) :
                tmp = tmp + self.coeff[i][j]*x[j]
            x[i] = (self.coeff[i][self.systemSize]-tmp) / self.coeff[i][i]

        self.isSolutionFounded = True

        for i in range(0,self.systemSize) :
            self.solutions[i] = x[i];
        # self.solutions = x

    def showSolutions(self) :
        if self.isSolutionFounded :
            print("Р РµС€РµРЅРёСЏ РІРІРµРґРµРЅРЅРѕР№ СЃРёСЃС‚РµРјС‹: ");
            for i in range(self.systemSize) :
                print(f"x[{i}]= {self.solutions[i]}")

if __name__ == "__main__" :
    scriptName, mode = argv

    gaussMethod = GaussMethod()
    if mode == 'user' :
        gaussMethod.initializationFromUserInput()
    elif mode == 'file' :
        gaussMethod.initializationFromFile()
    gaussMethod.showSystem()
    gaussMethod.applyGaussMethod()
    gaussMethod.showSolutions()
