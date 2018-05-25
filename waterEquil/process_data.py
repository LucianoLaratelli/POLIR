import sys
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) != 2:
        print "Needs exactly one file name as input."
        sys.exit(1)
    with open(sys.argv[1], 'r') as program:
        data = program.readlines()

    col_num = len(data[0].split())

    #POLIR code makes it a pain to output step numbers, so we
    #prepend them if necessary
    #the magic number `4` is the result of having the dipoles for each
    #of the x, y, and z directions, which make for three columns
    #the step number is the fourth
    if col_num != 4:
        with open(sys.argv[1], 'w') as program:
            for (number, line) in enumerate(data):
                program.write('%d  %s' % (number + 1, line))


    #heavy lifting handled by C here
    call(["./get_TCF", sys.argv[1]])

    correlated_file = sys.argv[1] + ".CORRELATED"

    with open(correlated_file, 'r') as corr:
        correlated = corr.readlines()

    fourier_of_correlated = np.fft.fft(correlated)
    frequency_of_fourier = np.fft.fftfreq(len(fourier_of_correlated))

    plt.plot(frequency_of_fourier, fourier_of_correlated)

    plt.show()


if __name__ == "__main__":
    main()
