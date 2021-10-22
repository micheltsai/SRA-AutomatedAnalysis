from __future__ import print_function

import csv
import datetime
import multiprocessing
import os
import shlex
import subprocess
import sys
import time
import traceback
from pathlib import Path
import analysisv4
import utils_


def run_cmd(cmd):
    cmd=shlex.split(cmd)
    p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #print (cmd)
    print("--------------------------------------\nSubprogram output:\n")
    while p.poll() is None:
        #progress_bar("sub excuting")
        line = p.stdout.readline()
        line = line.strip()
        if line:
            line_=line.decode().split("\n")
            for s in line_:
                print (str("{}\n".format(s)))
            sys.stdout.flush()
            sys.stderr.flush()
    if p.returncode ==0:
        print ("Subprogram sucess")
    else:
        print ("Subprogram failed")

    print ("-------------------------\n")
    return p

def progress_bar(Category):
    for i in range(1, 101):
        print("\r{}: ".format(Category),end="")
        print("[{}] {}%".format("*" * (i // 2), i), end="")
        sys.stdout.flush()
        time.sleep(0.02)
    print ("\n")

def Download(x):
    one_ = time.time()
    #print(
    #   "---------------------\n---------------------[ {} / {} ]---------------------\n".format(num + 1,
    #                                                                                           len(idlist)))
    #num += 1
    print("x = {}".format(x))
    # outdir__ = os.path.join(output, "out")
    outdir__ = os.path.join(new_outdir, "Assembled")

    final_dir = os.path.join(outdir__, "{}_contig.fa".format(x))
    if os.path.isfile(final_dir):
        print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
    else:
        utils_.prefetch_sra(x, sra_dir)
        print("Download {}\n.".format(x))
        #with open(Downloadcheck_log, "a+") as f:
        #    f.write("{}\n".format(x))
    print('Done,total cost', time.time() - one_, 'secs')
    print("###########################################################")



def Assembled(x):
    final_dir = os.path.join(ass_dir, "{}_contig.fa".format(x))
    if os.path.isfile(final_dir):
        print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
    else:
        # one_run_prefetch = time.time()
        # utils_.prefetch_sra(x,sra_dir)
        one_run_ass = time.time()
        sra_file = os.path.join(sra_dir, "{}/{}.sra".format(x, x))
        print(sra_file)
        if os.path.isfile(sra_file):
            print("have {}.sra file.\n".format(x))
        else:
            utils_.prefetch_sra(x, sra_dir)
            print("not found {}.sra, Download now\n".format(x))

        utils_.run_for_114(x, sra_dir, fastq_dir, assemble_dir, new_outdir, thread, gsize, start, check_log)
        ### unnecessary ERR file
        #ERR_path = os.path.join(os.path.abspath(os.getcwd()), x)
        #print("suERR_path: ", ERR_path, "\n")
        # print ("shutil.rmtree({})\n".format(current_path))
        #utils_.run_cmd2("rm -rf {}".format(current_path))
        #print("remove {}\n".format(current_path))



def SRA_Analysis(x):
    SRA_start=time.time()
    Download(x)
    Assembled(x)
    #####
    with open("./threads_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "{}".format(x), "time": str(time.time() - SRA_start)})

    genome = os.path.join(ass_dir, "{}_contig.fa".format(x))
    qual_cmd = "python3 QualityCheckv3-124.py -r {} -g {} -db {} -m {} -o {}".format(ref_dir,genome , buscoDB, buscoMode,new_outdir)
    try:
    #    targetPath = run_cmd(qual_cmd)
    except Exception as e:
        error_class = e.__class__.__name__  # 取得錯誤類型
        detail = e.args[0]  # 取得詳細內容
        cl, exc, tb = sys.exc_info()  # 取得Call Stack
        lastCallStack = traceback.extract_tb(tb)[-1]  # 取得Call Stack的最後一筆資料
        fileName = lastCallStack[0]  # 取得發生的檔案名稱
        lineNum = lastCallStack[1]  # 取得發生的行號
        funcName = lastCallStack[2]  # 取得發生的函數名稱
        errMsg = "File \"{}\", line {}, in {}: [{}] {}".format(fileName, lineNum, funcName, error_class, detail)
        print(errMsg)
        sys.exit(e)
    #print("targetPAth = {}\n######\n".format(targetPath.encode("utf-8")))

    #######


    with open(check_log,"a+") as f:
        f.write("Run {} is ok.\n".format(x))
    return 0
if __name__ == '__main__':
    start=time.time()
    Month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")
    ########################
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

    thread = setList[1].strip("\n").split("=")[1]
    gsize = setList[4].strip("\n").split("=")[1]
    n = setList[7].strip("\n").split("=")[1]
    outdir = setList[10].strip("\n").split("=")[1]
    utils_.mkdir_join(outdir)
    ref_dir = setList[12].strip("\n").split("=")[1]
    buscoDB = setList[13].strip("\n").split("=")[1]
    buscoMode = setList[14].strip("\n").split("=")[1]
    mlstS = setList[16].strip("\n").split("=")[1]
    amrS = setList[17].strip("\n").split("=")[1]

    #####################
    for mon in range(7, 8):
        for d in range(1, Month[mon] + 1):
            pattern = "salmonella enterica[ORGN] AND illumina[PLAT] AND wgs[STRA] AND genomic[SRC] AND paired[LAY]"
            ds = time.time()

            date = datetime.date(2020, mon + 1, d).strftime("%Y/%m/%d")
            #temp="{}/{}/{}".format(str(2020),str(mon+1),str(d))
            ######
            pdat = date.replace("/", "")
            new_outdir = os.path.join(outdir, pdat)
            utils_.mkdir_join(new_outdir)
            print("output: {}\n".format(new_outdir))

            #Downloadcheck_log = os.path.join(new_outdir, "Downloadcheck.log")
            check_log =os.path.join(new_outdir,"check.log")
            # commit
            # run_cmd2("touch {}".format("check.log"))
            #myfile = Path(Downloadcheck_log)
            #myfile.touch(exist_ok=True)
            #f = open(Downloadcheck_log, 'a+')
            #t = str(datetime.datetime.now()).split(".")[0]
            #f.write(t)
            #f.write("\n")
            #f.close()

            myfile2 = Path(check_log)
            myfile2.touch(exist_ok=True)
            with open(check_log,"a+") as f:
                f.write(str(datetime.datetime.now()).split(".")[0])
                f.write("\n")

            # mkdir_join(log_dir)
            # check_log = os.path.join(log_dir, "check.log")

            # print(date)

            pattern, count = utils_.count_egquery(pattern, date, date)
            print("pattern: {}\ncount: {}\n".format(pattern, count))

            i_e_ = time.time()
            idlist = utils_.IdList_esearch(pattern, 'sra', count)

            print(idlist)

            runinfo = utils_.Get_RunInfo(idlist)
            # progress_bar("get SRAfile name List stored in run_list")
            run_list = list(runinfo['Run'])  # get SRAfile nameList stored in run_list
            print("runinfo: {}\n run_list: {}\n".format(runinfo, run_list))

            sra_dir = os.path.join(new_outdir, "sra")  # .sra file
            utils_.mkdir_join(sra_dir)
            ass_dir = os.path.join(new_outdir, "Assembled")
            utils_.mkdir_join(ass_dir)
            fastq_dir = os.path.join(new_outdir, 'fastq')
            assemble_dir = os.path.join(new_outdir, "assembly_result")

            read_log_ = time.time()
            f = open(check_log, 'r+')
            line = f.readlines()
            print("check log :{}\n".format(line))
            f.close()

            for s in line:
                print("{}\n".format(s))
            finish = list(filter(lambda x: len(x.split(" ")) >= 4, line))
            finish_run = list(map(lambda x: x.split(" ")[1], finish))
            need_run = list(filter(lambda x: x not in finish_run, run_list))
            print("finish: {}\nfinish_run: {}\nneed_run".format(finish, finish_run, need_run))
            print("finish length: {}\nfinish_run length: {}\nneed_run length: ".format(len(finish), len(finish_run),
                                                                                       len(need_run)))
            print("Toal", len(need_run), "sra runs need to downlaod.")

            num = len(finish_run)
            progress_list = []
            prog_num = 0
            for k in need_run:

                progress_list.append(multiprocessing.Process(target=SRA_Analysis, args=(k,)))
                progress_list[prog_num].start()
                prog_num += 1

            for i in range(prog_num):
                progress_list[i].join()


            with open("./Automate_check.log", "a+") as f:
                f.write("{}:{}:{}\n".format(date, time.time() - ds, time.time() - start))
            print("Download all {} ".format(date),'Done,total cost', time.time() - ds, 'secs')
            time.sleep(5)
    print('Done,total cost', time.time() - start, 'secs')
    ##########





