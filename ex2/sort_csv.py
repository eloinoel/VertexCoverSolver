import sys
import csv
import operator

argv = sys-argv
reader = csv.reader(open(argv[1]), delimenter=";")
sortedList = sorted(reader, key=argv[2])