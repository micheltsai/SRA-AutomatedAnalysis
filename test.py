import os
import subprocess
import time

import xml.etree.cElementTree as ET


def main():
    sd_Y=2008
    ed_Y=2021
    sd_M=1
    ed_M=4
    sd_D=1
    ed_D=15
    Month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    for yy in range(sd_Y,ed_Y+1):
        Month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        if (yy % 4) == 0:
            if (yy % 100) == 0:
                if (yy % 400) == 0:
                    Month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            else:
                Month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        if sd_Y == yy:
            sM = sd_M
            eM = 12
            sD = sd_D
        elif yy != ed_Y:
            sM = 1
            eM = 12
            sD = 1
        else:
            sM = 1
            eM = ed_M
            sD = 1
        for mon in range(sM,eM+1):
            if yy!=ed_Y:
                eD=Month[mon-1]
            elif mon != eM:
                eD=Month[mon-1]
            else:
                eD=ed_D
            for d in range(sD,eD+1):
                print("{}/{}/{}".format(yy,mon,d))
    return 0
if __name__ == '__main__':
    main()