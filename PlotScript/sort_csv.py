import sys
import csv
import operator

""" argv = sys.argv
reader = csv.reader(open(argv[1]), delimiter=";")
sortedList = sorted(reader, key=operator.itemgetter(1))
writer = csv.writer("test_sorted.csv") """


with open(sys.argv[1], 'r') as in_file:
    reader = csv.reader(in_file, delimiter=";")
    header = next(reader)
    data = sorted(reader, key=lambda row: float(row[1]))

with open('sorted.csv', 'w', newline='') as out_file:
    out_writer = csv.writer(out_file, delimiter=";")
    out_writer.writerow(header)
    out_writer.writerows(data)
