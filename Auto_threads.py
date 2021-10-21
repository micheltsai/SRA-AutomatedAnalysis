from __future__ import print_function

import datetime
import multiprocessing
import os
import shlex
import subprocess
import sys
import time
import traceback

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

def Download(date):
    print("Dowload threadPID",os.getgid())
    print("Download 父進程編號",os.getppid())
    run_cmd("python3 DownloadFile.py --PDAT {}".format(date))
def Analysis(date):
    print("Analysis threadPID", os.getgid())
    print("Analysis 父進程編號", os.getppid())
    run_cmd("python3 SRA_Analysisv2.py --PDAT {}".format(date))

    return 0
if __name__ == '__main__':
    start=time.time()
    Month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    try:
        for x in range(0, 12):
            for d in range(1, Month[x] + 1):
                ds = time.time()
                c = datetime.datetime(2020, x + 1, d)
                tmp = c.strftime("%Y/%m/%d")
                download_prog=multiprocessing.Process(target=Download, args=(tmp,))

                analysis_prog=multiprocessing.Process(target=Analysis, args=(tmp,))

                download_prog.start()
                progress_bar("Download")
                analysis_prog.start()
                progress_bar("Analysis")

                with open("./Automate_check.log", "a+") as f:
                    f.write("{}:{}:{}\n".format(tmp, time.time() - ds, time.time() - start))
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

