import datetime
import time

import utils_
def main():
    start = time.time()
    cmd="python3 SRA_Analysis.py --PDAT "
    Month=[31,28,31,30,31,30,31,31,30,31,30,31]
    for x in range(0,12):
        for d in range(1,Month[x]+1):
            ds=time.time()
            c = datetime.datetime(2020,x+1,d)
            tmp=c.strftime("%Y/%m/%d")
            cmd2=cmd+tmp
            print(cmd2)
            utils_.run_cmd3(cmd2)
            with open("./Automate_check.log","a+") as f:
                f.write("{}:{}:{}\n".format(tmp, time.time()-ds,time.time()-start))

    print('Done,total cost', time.time() - start, 'secs')
    return 0


if __name__ == '__main__':
    main()