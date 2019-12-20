import numpy as np
import matplotlib.pyplot as plt
from sys import argv

def show_chart(*coords) :
    x, y = coords

    figure, axis = plt.subplots(figsize=(10, 10))

    axis.plot(x,y, label='Яма')

    axis.set_ylabel('$Y axis (cm)$')
    axis.set_xlabel('$X axis (cm)$')
    axis.set_xlim(left=0,right=20)
    axis.set_ylim(bottom=0,top=25)
    axis.minorticks_on()
    axis.grid(which = "major",color='orange',linewidth = 1)
    axis.grid(which = "minor",color='orange',linestyle = ":")
    figure.tight_layout()

    axis.legend()

    # plt.draw()
    plt.show()

if __name__ == "__main__" :
    scriptName = argv

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

    show_chart(x,y)
