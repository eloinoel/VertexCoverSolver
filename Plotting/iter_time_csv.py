import sys
import numpy as np

def main():
    # Using readlines()
    file = open(sys.argv[1], 'r')
    out_file = open("ittime_"+sys.argv[1], 'w')

    line = file.readline()
    split = line.split(';')
    out_file.write(split[0]+";"+split[1]+"\n")

    for file_line in file:
        split = file_line.split(';')
        exec_time = split[1]
        iterations = split[5]
        # timeout
        if(exec_time == ""):
            #out_file.write(file_line)
            out_file.write(split[0]+";0")
            out_file.write("\n")
            continue
        # write test file name in column 0
        #out_file.write(split[0]+";")
        try:
            #print("Number "+exec_time)
            # valid time
            iter_time = np.double(exec_time) / np.double(iterations)
            out_file.write(split[0]+";")
            out_file.write(str(iter_time) + "\n")
        except:
            #other error
            out_file.write(split[0]+";0")
            out_file.write("\n")
            continue

    file.close()
    out_file.close()

if __name__ == "__main__":
    main()