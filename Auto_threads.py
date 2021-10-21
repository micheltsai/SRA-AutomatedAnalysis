from __future__ import print_function

import datetime
import multiprocessing
import os
import shlex
import subprocess
import sys
import time
import traceback
from pathlib import Path

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

def Download():
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
        with open(check_log, "a+") as f:
            f.write("{}\n".format(x))
    print('Done,total cost', time.time() - one_, 'secs')


def Analysis(date):
    print("Analysis threadPID", os.getgid())
    print("Analysis 父進程編號", os.getppid())
    run_cmd("python3 SRA_Analysisv2.py --PDAT {}".format(date))

    return 0
if __name__ == '__main__':
    start=time.time()
    Month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")
    pattern = "salmonella enterica[ORGN] AND illumina[PLAT] AND wgs[STRA] AND genomic[SRC] AND paired[LAY]"
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
    try:
        for x in range(0, 12):
            for d in range(1, Month[x] + 1):
                ds = time.time()
                date = datetime.datetime(2020, x + 1, d).strftime("%Y/%m/%d")
                ######
                pdat = date.replace("/", "")
                new_outdir = os.path.join(outdir, pdat)
                utils_.mkdir_join(new_outdir)
                print("output: {}\n".format(new_outdir))

                check_log = os.path.join(new_outdir, "checkDownload.log")
                # commit
                # run_cmd2("touch {}".format("check.log"))
                myfile = Path(check_log)
                myfile.touch(exist_ok=True)
                f = open(check_log, 'a+')
                t = str(datetime.datetime.now()).split(".")[0]
                f.write(t)
                f.write("\n")
                f.close()

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

                sra_dir = os.path.join(outdir, "sra")  # .sra file
                utils_.mkdir_join(sra_dir)

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
                progress_list=[]
                prog_num=0
                for x in need_run:
                    progress_list.append(multiprocessing.Process(target=Download, args=(x,)))
                    progress_list[prog_num].start()
                    prog_num+=1

                for i in range(prog_num):
                    progress_list[i].join()
                print("Download all {}".format(date))
                print('Done,total cost', time.time() - start, 'secs')
                ##########


                with open("./Automate_check.log", "a+") as f:
                    f.write("{}:{}:{}\n".format(date, time.time() - ds, time.time() - start))
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

