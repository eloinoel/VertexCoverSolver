import sys
import numpy as np

def main():
    # Using readlines()
    file1 = open(sys.argv[1], 'r')
    file2 = open(sys.argv[2], 'r')
    file3 = open(sys.argv[3], 'r')
    file4 = open(sys.argv[4], 'r')
    out_file = open("intersection.csv", 'w')

    line1 = file1.readline()
    line2 = file2.readline()
    line3 = file3.readline()
    line4 = file4.readline()
    split1 = line1.split(';')
    out_file.write(split1[0]+";"+"time1;time2;time3;time4"+"\n")

    finished = False

    while (line1):
        line1 = file1.readline()
        line2 = file2.readline()
        line3 = file3.readline()
        line4 = file4.readline()

        split1 = line1.split(';')
        split2 = line2.split(';')
        split3 = line3.split(';')
        split4 = line4.split(';')

        try:

            while (int(line1[0]) < int(line2[0]) or int(line1[0]) < int(line3[0]) or int(line1[0]) < int(line4[0])):
                line1 = file1.readline()
                split1 = line1.split(';')
                if(not line1):
                    finished = True
                    break

            while (int(line2[0]) < int(line1[0])):
                line2 = file2.readline()
                split2 = line2.split(';')
                if(not line2):
                    finished = True
                    break

            while (int(line3[0]) < int(line1[0])):
                line3 = file3.readline()
                split3 = line3.split(';')
                if(not line3):
                    finished = True
                    break

            while (int(line4[0]) < int(line1[0])):
                line4 = file4.readline()
                split4 = line4.split(';')
                if(not line4):
                    finished = True
                    break
            if(finished):
                file1.close()
                file2.close()
                file3.close()
                file4.close()
                out_file.close()
                return

            exec_time1 = split1[1]
            exec_time2 = split2[1]
            exec_time3 = split3[1]
            exec_time4 = split4[1]
            # timeout or error
            if(exec_time1 == "" or exec_time2 == "" or exec_time3 == "" or exec_time4 == ""):
                #out_file.write("\n")
                continue
            np.double(exec_time1)
            np.double(exec_time2)
            np.double(exec_time3)
            np.double(exec_time4)
            if(np.double(exec_time1) == 0 or np.double(exec_time2) == 0 or np.double(exec_time3) == 0 or np.double(exec_time4) == 0):
                #out_file.write("\n")   // TODO: this doesnt filter 0'es out, idk why
                continue
            # write test file name in column 0
            #out_file.write(split[0]+";")
            #print("Number "+exec_time)
            # valid time
            out_file.write(split1[0]+";")
            out_file.write(split1[1][:-1]+";")
            out_file.write(split2[1][:-1]+";")
            out_file.write(split3[1][:-1]+";")
            out_file.write(split4[1][:-1]+"\n")
        except:
            continue

    file1.close()
    file2.close()
    file3.close()
    file4.close()
    out_file.close()

if __name__ == "__main__":
    main()