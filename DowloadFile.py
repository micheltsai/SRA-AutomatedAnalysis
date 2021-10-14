from __future__ import print_function
# -*- coding: utf-8 -*-
import csv
import datetime
import shlex
import time
import os
import re
import sys
import shutil
import argparse
import subprocess
import traceback
from pathlib import Path
import pandas as pd
import utils_

### v3 is modify outdir : args.outdir/PDAT/Assembled

def main():
    # 計算時間
    start = time.time()
    utils_.progress_bar("read command")
    parser = argparse.ArgumentParser("Download_avaliable_runs_from_NCBI_sra_database_wit_PDAT.")
    parser.add_argument("--pattern",
                        default="salmonella enterica[ORGN] AND illumina[PLAT] AND wgs[STRA] AND genomic[SRC] AND paired[LAY]",
                        help="Searching condition.")
    # PDAT格式：YYYY/MM/DD
    parser.add_argument("--PDAT", required=True, help="Publication Date[PDAT] of Runs.")
    #parser.add_argument("--output", required=True, help="Folder to save Contigs.fa after assembly of eah runs.")
    parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    parser.add_argument("--n", default=3, help="count of download sra file eahc time", type=int)
    args = parser.parse_args()

    pattern = args.pattern
    # print(pattern)
    date = args.PDAT
    #output = args.output
    threads = args.threads
    gsize = args.gsize
    n = args.n
    utils_.progress_bar("read arguments")


    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")
    ## read SRAsetting.txt
    utils_.progress_bar("read SRAsetting.txt")
    setting_path = os.path.join(current_path, "SRAsettings.txt")
    with open(setting_path, "r") as f:
        setList = f.readlines()

    print(setList)
    i = 0
    for line in setList:
        line = line.strip("\n")
        line_ = line.split("=")
        if line != "" and len(line_) == 2:
            print(line_)
            print("line{}. {}:{}\n".format(i, line_[0], line_[1]))
        i += 1

    outdir = setList[10].strip("\n").split("=")[1]

    utils_.mkdir_join(outdir)
    pdat = date.replace("/", "")
    outdir = os.path.join(outdir, pdat)
    utils_.mkdir_join(outdir)

    print("output: {}\n".format(outdir))

    check_log = os.path.join(outdir, "checkDownload.log")
    # commit
    #run_cmd2("touch {}".format("check.log"))
    myfile = Path(check_log)
    myfile.touch(exist_ok=True)
    f = open(check_log, 'a+')
    t = str(datetime.datetime.now()).split(".")[0]
    f.write(t)
    f.write("\n")
    f.close()

    #mkdir_join(log_dir)
    #check_log = os.path.join(log_dir, "check.log")


    # print(date)
    pattern, count = utils_.count_egquery(pattern, date, date)
    print ("pattern: {}\ncount: {}\n".format(pattern,count))

    i_e_=time.time()
    idlist = utils_.IdList_esearch(pattern, 'sra', count)

    print(idlist)

    runinfo = utils_.Get_RunInfo(idlist)
    #progress_bar("get SRAfile name List stored in run_list")
    run_list = list(runinfo['Run']) #get SRAfile nameList stored in run_list
    print("runinfo: {}\n run_list: {}\n".format(runinfo, run_list))


    sra_dir = os.path.join(outdir, "sra")  # .sra file
    utils_.mkdir_join(sra_dir)

    read_log_=time.time()
    f = open(check_log, 'r+')
    d = f.readlines()
    print("check log :{}\n".format(d))
    f.close()


    for s in d:
        print ("{}\n".format(s))
    finish = list(filter(lambda x: len(x.split(" ")) >= 4, d))
    finish_run = list(map(lambda x: x.split(" ")[1], finish))
    need_run = list(filter(lambda x: x not in finish_run, run_list))
    print("finish: {}\nfinish_run: {}\nneed_run".format(finish,finish_run,need_run))
    print("finish length: {}\nfinish_run length: {}\nneed_run length: ".format(len(finish), len(finish_run), len(need_run)))
    print("Toal", len(need_run), "sra runs need to downlaod.")

    num = len(finish_run)
    for x in need_run:
        print("---------------------\n---------------------[ {} / {} ]---------------------\n".format(num, len(idlist)))
        num += 1
        print("x = {}".format(x))
        # outdir__ = os.path.join(output, "out")
        outdir__ = os.path.join(outdir, "Assembled")
        final_dir = os.path.join(outdir__, "{}_contig.fa".format(x))
        if os.path.isfile(final_dir):
            print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
        else:
            utils_.prefetch_sra(x, sra_dir)
            print("Download {}\n.".format(x))
            with open(check_log,"a+") as f:
                f.write("Run {} is ok.\n".format(num,x))





    print("Download all {}".format(date))

if __name__ == '__main__':
    main()