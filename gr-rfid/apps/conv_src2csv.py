import sys
import csv
import numpy as np

def conv_src2csv(source_file, dest_file):

    src = open(source_file, 'r')
    dest = open(dest_file, 'w')
    csv_wr = csv.writer(dest, delimiter=',')

    signals = []
    while True:
        samples = src.readline()
        if not samples: break

        sample = samples.split(' ')
        sample.pop()
        signals.append(sample)

    for signal in signals:
        csv_wr.writerow(signal)

    src.close()
    dest.close()

if __name__ == "__main__":
    
    SOURCE_FILE_NAME = sys.argv[1]
    DEST_FILE_NAME = sys.argv[1] + ".csv"

    conv_src2csv(SOURCE_FILE_NAME, DEST_FILE_NAME)
