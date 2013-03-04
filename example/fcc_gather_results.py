import sys
import csv

if __name__ == "__main__":
    outfile = sys.argv[1]
    with open(outfile, 'wb') as csvout:
        writer = csv.writer(csvout)
        for n in range(256):
            filename = './fcc_6_{}.csv'.format(n)
            with open(filename, 'rb') as csvin:
                for row in csv.reader(csvin):
                    writer.writerow(row)
