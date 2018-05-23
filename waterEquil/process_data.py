import sys
from subprocess import call
import numpy as np

def main():
    if(sys.argv[1] == None):
        print "Needs exactly one file name as input."
        sys.exit(1)
    with open(sys.argv[1], 'r') as program:
        data = program.readlines()

    col_num = len(data[0].split())


    if(col_num != 4):
        with open(sys.argv[1], 'w') as program:
            for (number, line) in enumerate(data):
                program.write('%d  %s' % (number + 1, line))
    

    call(["./get_TCF", sys.argv[1]])

    correlated_file = sys.argv[1] + ".CORRELATED"

    with open(correlated_file, 'r') as corr:
        correlated = corr.readlines()

    print np.fft.fft(correlated)


if __name__ == "__main__":
    main()
