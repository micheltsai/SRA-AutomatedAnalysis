from __future__ import print_function
# -*- coding: utf-8 -*-
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
    #parser.add_argument("--sra_dir", required=True, help="Temp folder to save runs.")
    #parser.add_argument("--log", required=True, help="The file to recored runs which finished assembly each time.")
    parser.add_argument("--output", required=True, help="Folder to save Contigs.fa after assembly of eah runs.")
    #parser.add_argument("--assembly_dir", default=".", help="Temp folder to save assembly.")
    #parser.add_argument("--tmpdir", default="/tmp", help="Directory of temp folder default: '/tmp'")
    parser.add_argument("--gsize", default='', help="Estimated genome size(MB) eg. 3.2M. default: ''")
    parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    parser.add_argument("--n", default=3, help="count of download sra file eahc time", type=int)
    args = parser.parse_args()

    pattern = args.pattern
    # print(pattern)
    date = args.PDAT
    #sra_dir = args.sra_dir
    #assem_dir = args.assembly_dir
    #log_dir = args.log
    output = args.output
    #tmpdir = args.tmpdir
    threads = args.threads
    gsize = args.gsize
    n = args.n
    utils_.progress_bar("read arguments")
    utils_.mkdir_join(output)
    check_log = os.path.join(output,"check.log")
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
    idlist = utils_.IdList_esearch(pattern, 'sra', count)
    print(idlist)
    runinfo = utils_.Get_RunInfo(idlist)
    #progress_bar("get SRAfile name List stored in run_list")
    run_list = list(runinfo['Run']) #get SRAfile nameList stored in run_list
    print("runinfo: {}\n run_list: {}\n".format(runinfo, run_list))

    f = open(check_log, 'r+')
    d = f.readlines()
    print("chcke log :{}\n".format(d))
    f.close()

    for s in d:
        print ("{}\n".format(s))
    finish = list(filter(lambda x: len(x.split(" ")) >= 4, d))
    finish_run = list(map(lambda x: x.split(" ")[1], finish))
    need_run = list(filter(lambda x: x not in finish_run, run_list))
    print("finish: {}\nfinish_run: {}\nneed_run".format(finish,finish_run,need_run))
    print("finish length: {}\nfinish_run length: {}\nneed_run length: ".format(len(finish), len(finish_run), len(need_run)))
    print("Toal", len(need_run), "sra runs need to downlaod.")

    utils_.mkdir_join(output)

    sra_dir = os.path.join(output, "sra")  # .sra file
    utils_.mkdir_join(sra_dir)

    #k,每三個一輪迴
    k = list(range(0, len(need_run), n))
    print (k)

    for i in k:
        run_id = need_run[i:i+n]
        print("###### i = {}\n".format(i))
        print("run_id: {}\n".format(run_id))
        time.sleep(1)
        num=len(finish_run)
        for x in run_id:
            print ("---------------------\n---------------------[ {} / {} ]---------------------\n".format(num,len(idlist)))
            num+=1
            print ("x = {}".format(x))
            outdir__ = os.path.join(output, "out")
            final_dir = os.path.join(outdir__, "{}_contig.fa".format(x))
            if os.path.isfile(final_dir):
                print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
            else:
                utils_.prefetch_sra(x,sra_dir)
                utils_.run_for_114(x,sra_dir,output,threads,gsize,start,check_log)
                current_path = os.path.join(os.path.abspath(os.getcwd()), x)
                print("current_path: ", current_path, "\n")
                # print ("shutil.rmtree({})\n".format(current_path))
                utils_.run_cmd2("rm -rf {}".format(current_path))
                print ("remove {}\n".format(current_path))
        print("shutil.rmtree(sra_dir)\n")
        shutil.rmtree(sra_dir)
        utils_.mkdir_join(sra_dir)








    print('Done,total cost', time.time() - start, 'secs')


if __name__ == '__main__':
    main()
