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

    print("output: {}\n".format(output))

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
    c_e_=time.time()
    pattern, count = utils_.count_egquery(pattern, date, date)
    print ("pattern: {}\ncount: {}\n".format(pattern,count))
    with open("./ana_time.csv","a+")as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "count_egquery", "time": str(time.time()-c_e_)})
    i_e_=time.time()
    idlist = utils_.IdList_esearch(pattern, 'sra', count)
    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "IdList_esearch", "time": str(time.time() - i_e_)})
    print(idlist)
    g_r_=time.time()
    runinfo = utils_.Get_RunInfo(idlist)
    #progress_bar("get SRAfile name List stored in run_list")
    run_list = list(runinfo['Run']) #get SRAfile nameList stored in run_list
    print("runinfo: {}\n run_list: {}\n".format(runinfo, run_list))
    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "get_RunINfo", "time": str(time.time() - g_r_)})


    utils_.mkdir_join(output)
    fastq_dir = os.path.join(output, 'fastq')
    assemble_dir = os.path.join(output, "assembly_result")
    sra_dir = os.path.join(output, "sra")  # .sra file
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

    if len(need_run) == 0:
        #utils_.progress_bar("ALL is assembled.")
        shutil.rmtree(sra_dir)
        shutil.rmtree(fastq_dir)
        shutil.rmtree(assemble_dir)
        print("ALL is assembled.\n")

        print("********************  ASSEMBLED END  ************************\n\n")
        return 0




    #k,每三個一輪迴
    k = list(range(0, len(need_run), n))
    print (k)
    num = len(finish_run)

    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "read_check and get need_runList", "time": str(time.time() - read_log_)})

    for i in k:
        run_id = need_run[i:i+n]
        print("###### i = {}\n".format(i))
        print("run_id: {}\n".format(run_id))
        time.sleep(1)
        for x in run_id:
            print ("---------------------\n---------------------[ {} / {} ]---------------------\n".format(num,len(idlist)))
            num+=1
            print ("x = {}".format(x))
            #outdir__ = os.path.join(output, "out")
            outdir__ = os.path.join(output, "Assembled")

            final_dir = os.path.join(outdir__, "{}_contig.fa".format(x))
            if os.path.isfile(final_dir):
                print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
            else:
                #one_run_prefetch = time.time()
                #utils_.prefetch_sra(x,sra_dir)
                one_run_ass=time.time()
                sra_file=os.path.join(sra_dir,"{}/{}.sra".format(x,x))
                if os.path.isfile(sra_file):
                    print("have {}.sra file.\n".format(x))
                else:
                    utils_.prefetch_sra(x, sra_dir)
                    print("not found {}.sra, Download now\n".format(x))

                utils_.run_for_114(x,sra_dir,fastq_dir,assemble_dir,output,threads,gsize,start,check_log)


                current_path = os.path.join(os.path.abspath(os.getcwd()), x)
                print("current_path: ", current_path, "\n")
                three_run = time.time()
                # print ("shutil.rmtree({})\n".format(current_path))
                utils_.run_cmd2("rm -rf {}".format(current_path))
                print ("remove {}\n".format(current_path))

            with open("./ana_time.csv", "a+") as f:
                fieldnames = ["func", "time"]
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow({"func": "one file assembled", "time": str(time.time() - one_run_ass)})

            if num==5:
                print("break for loop\n")
                break
        print("shutil.rmtree(sra_dir)\n")
        shutil.rmtree(sra_dir)
        utils_.mkdir_join(sra_dir)
        with open("./ana_time.csv", "a+") as f:
            fieldnames = ["func", "time"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow({"func": "three so removed '/sra'", "time": str(time.time() - three_run)})

    if num == count:
        shutil.rmtree(sra_dir)
        shutil.rmtree(fastq_dir)
        shutil.rmtree(assemble_dir)
        with open(check_log,"a+") as f:
            f.write("ALL({}/{}) is ok.\n".format(num,count))

    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "all time", "time": str(time.time() - start)})

    print('Done,total cost', time.time() - start, 'secs')

    print("********************  ASSEMBLED END  ************************\n\n")
    return 0
if __name__ == '__main__':
    main()
