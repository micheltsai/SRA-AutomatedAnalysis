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
    pdat = date.replace("/", "")
    output = os.path.join(output, pdat)
    utils_.mkdir_join(output)
    print("output: {}\n".format(output))

    check_log = os.path.join(output, "Assemblecheck.log")
    run_txt=os.path.join(output,"checkDownload.log")
    fastq_dir = os.path.join(output, 'fastq')
    assemble_dir = os.path.join(output, "assembly_result")
    sra_dir = os.path.join(output, "sra")  # .sra file
    utils_.mkdir_join(sra_dir)

    read_log_ = time.time()
    myfile1 = Path(check_log)
    myfile1.touch(exist_ok=True)
    f = open(check_log, 'r')
    d = f.readlines()
    print("check log :{}\n".format(d))
    f.close()

    myfile = Path(run_txt)
    myfile.touch(exist_ok=True)

    with open(run_txt, "r")as f:
        run_list=f.readlines()
        print(run_list)



    for s in d:
        print ("{}\n".format(s))
    finish = list(filter(lambda x: len(x.split(" ")) >= 4, d))
    finish_run = list(map(lambda x: x.split(" ")[1], finish))
    need_run = list(filter(lambda x: x not in finish_run, run_list.split(" ")[1]))
    print("finish: {}\nfinish_run: {}\nneed_run".format(finish,finish_run,need_run))
    print("finish length: {}\nfinish_run length: {}\nneed_run length: ".format(len(finish), len(finish_run), len(need_run)))
    print("Toal", len(need_run), "sra runs need to downlaod.")

    if len(need_run) == 0:
        #utils_.progress_bar("ALL is assembled.")
        #shutil.rmtree(sra_dir)
        #shutil.rmtree(fastq_dir)
        #shutil.rmtree(assemble_dir)
        print("ALL is assembled.\n")

        print("********************  ASSEMBLED END  ************************\n\n")
        return 0

        # k,每三個一輪迴
        k = list(range(0, len(need_run), n))
        print(k)
        num = len(finish_run)
        for i in k:
            run_id = need_run[i:i + n]
            print("###### i = {}\n".format(i))
            print("run_id: {}\n".format(run_id))
            time.sleep(1)
            for x in run_id:
                print("---------------------\n---------------------[ {} / {} ]---------------------\n".format(num,len(idlist)))
                num += 1
                print("x = {}".format(x))
                # outdir__ = os.path.join(output, "out")
                outdir__ = os.path.join(output, "Assembled")
                final_dir = os.path.join(outdir__, "{}_contig.fa".format(x))
                if os.path.isfile(final_dir):
                    print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
                else:
                    utils_.run_for_114(x, sra_dir, fastq_dir, assemble_dir, output, threads, gsize, start, check_log)
                    current_path = os.path.join(os.path.abspath(os.getcwd()), x)
                    print("current_path: ", current_path, "\n")
                    # print ("shutil.rmtree({})\n".format(current_path))
                    # utils_.run_cmd2("rm -rf {}".format(current_path))
                    # print ("remove {}\n".format(current_path))
                sra_file=os.path.join(sra_dir, x)
                print("shutil.rmtree(sra_dir)\n")
                shutil.rmtree(sra_file)
                print("shutil.rmtree({}.sra)\n".format(x))
                #with open(run_txt, "a+") as f:
                #    f.write("Run {} is ok.\n".format(x))

        if num == count:
            shutil.rmtree(sra_dir)
            shutil.rmtree(fastq_dir)
            shutil.rmtree(assemble_dir)
            print("ALL({}/{}) is ok.\n".format(num, count))
            print("shutil.rmtree(sra_dir)\n")
            print("shutil.rmtree(fastq_dir)\n")
            print("shutil.rmtree(assemble_dir)\n")
            #with open(check_log, "a+") as f:
                #f.write("ALL({}/{}) is ok.\n".format(num, count))

        print('Done,total cost', time.time() - start, 'secs\n')

        print("********************  ASSEMBLED END  ************************\n\n")
        return 0
if __name__ == '__main__':
    main()