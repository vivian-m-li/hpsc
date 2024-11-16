#!/usr/bin/python3

import os

import sys
import getopt
import glob

styles = {
    'runA': 'w l lw 3.5 lc "blue"',
    'runB2': 'w l lw 1.2 lc "yellow"' 
}

def writePCfile(num_PEs, name, time):
    g = open('pc_' + name, 'w')
    print ('splot \\' , file = g)
    for folder in ['runA', 'runB2']:
        for PE in range(num_PEs):
            fileName = f'{folder}/phi_{PE}_{time}.plt'
            if os.path.isfile(fileName):
                print(f'"{fileName}" {styles[folder]}, \\' , file = g)
    g.close()


def writePlotCmd(argv):
    num_PEs = 16

    # compare runA at t=0.025 to the first plots in runB2
    writePCfile(num_PEs, 'restart', 24)

    # compare runA at t=0.05 to runB2 at t=0.05
    writePCfile(num_PEs, 'end', 49)


if __name__ == "__main__":
    writePlotCmd(sys.argv[1:])
