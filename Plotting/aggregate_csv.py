import sys
import numpy as np

def main():
    # Using readlines()
    file = open(sys.argv[1], 'r')
    out_file = open("agg_"+sys.argv[1], 'w')

    agg_time = np.double(0.0)
    line = file.readline()
    split = line.split(';')
    out_file.write(split[0]+";"+split[1]+"\n")

    for file_line in file:
        split = file_line.split(';')
        exec_time = split[1]
        # timeout
        if(exec_time == ""):
            #out_file.write("\n")
            continue
        # write test file name in column 0
        #out_file.write(split[0]+";")
        try:
            #print("Number "+exec_time)
            # valid time
            agg_time += np.double(exec_time)
            out_file.write(split[0]+";")
            out_file.write(str(agg_time) + "\n")
        except:
            #other error
            #out_file.write(split[0]+";")
            #out_file.write("\n")
            1 == 1

    file.close()
    out_file.close()

if __name__ == "__main__":
    main()