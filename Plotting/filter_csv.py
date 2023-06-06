import sys

def main():
    # Using readlines()
    file = open(sys.argv[1], 'r')
    out_file = open("f_"+sys.argv[1], 'w')

    line = file.readline()
    split = line.split(';')
    out_file.write(line)

    for file_line in file:
        split = file_line.split(';')
        exec_time = split[1]
        # timeout
        if(exec_time == ""):
            #out_file.write("\n")
            continue
        # write test file name in column 0
        #out_file.write(line+"\n")
        try:
            # try to trigger exception
            i = float(exec_time)

            #print("Number "+exec_time)
            # valid time
            out_file.write(file_line)
        except:
            #other error
            #out_file.write(split[0]+";")
            #out_file.write("\n")
            continue

    file.close()
    out_file.close()

if __name__ == "__main__":
    main()