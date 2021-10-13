import argparse
import csv
import datetime
import os
import shutil
import sys
import time
import traceback
from pathlib import Path

import utils_

def main():
    start = time.time()


    # read command arguments------
    # get ref_path, qen_path, and outdir
    utils_.progress_bar("read command")
    parser = argparse.ArgumentParser("python3 Analysis.py ")
    parser.add_argument("--pattern",
                        default="salmonella enterica[ORGN] AND illumina[PLAT] AND wgs[STRA] AND genomic[SRC] AND paired[LAY]",
                        help="Searching condition.")
    # PDAT格式：YYYY/MM/DD
    parser.add_argument("--PDAT", required=True, help="Publication Date[PDAT] of Runs.")
    # parser.add_argument("-i", "--input", required=True, help="genome")
    # parser.add_argument("-o", "--outdir", required=True, help="Output folder")
    # MLST -s
    # parser.add_argument("--s", "--organism", help="MLST need -s organism")
    # Inc type -p
    # parser.add_argument("-p", "--plasmidfinderDB", required=True, help="Path of plasmidfinder database")
    # Resistance genes & Mutations identification -d
    # parser.add_argument("-d", "--amrfinderDB", required=True, help="Path of amrifinder database")
    # parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    # serotype
    # parser.add_argument("-m", "--mode", default="geno", required=False, help="busco mode")
    args = parser.parse_args()
    # referencelist=args.ref
    # busco_db=args.busco_datatbase
    date = args.PDAT
    pattern=args.pattern
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

    thread = setList[1].strip("\n").split("=")[1]
    gsize = setList[4].strip("\n").split("=")[1]
    n = setList[7].strip("\n").split("=")[1]
    outdir = setList[10].strip("\n").split("=")[1]
    ref_dir = setList[12].strip("\n").split("=")[1]
    buscoDB = setList[13].strip("\n").split("=")[1]
    buscoMode = setList[14].strip("\n").split("=")[1]
    mlstS = setList[16].strip("\n").split("=")[1]
    amrS = setList[17].strip("\n").split("=")[1]

    gen_dir = os.path.join(outdir, "Assembled")
    print(
        "thread: {}\ngsize: {}\nn:{}\noutdir: {}\nref_dir: {}\nbuscoDB: {}\nbusco_mode: {}\nmlstS: {}\namrS: {}\n".format(
            thread, gsize, n, outdir, ref_dir, buscoDB, buscoMode, mlstS, amrS))

    utils_.mkdir_join(outdir)
    # ./YYYYMMDD
    pdat = date.replace("/", "")
    outdir = os.path.join(outdir, pdat)
    utils_.mkdir_join(outdir)

    pattern, count = utils_.count_egquery(pattern, date, date)
    print("pattern: {}\ncount: {}\n".format(pattern, count))
    i_e_ = time.time()
    idlist = utils_.IdList_esearch(pattern, 'sra', count)
    print(idlist)

    runinfo = utils_.Get_RunInfo(idlist)
    # progress_bar("get SRAfile name List stored in run_list")
    run_list = list(runinfo['Run'])  # get SRAfile nameList stored in run_list
    run_txt=os.path.join(outdir,"run_list.txt")
    with open(run_txt, "a+") as f:
        for run in run_list:
            f.write(run)
            f.write("\n")

    print("runinfo: {}\n run_list: {}\n".format(runinfo, run_list))

    utils_.mkdir_join(outdir)
    fastq_dir = os.path.join(outdir, 'fastq')
    assemble_dir = os.path.join(outdir, "assembly_result")
    sra_dir = os.path.join(outdir, "sra")  # .sra file
    utils_.mkdir_join(sra_dir)

    ####check_log
    check_log = os.path.join(outdir, "final_check.log")

    # commit
    # run_cmd2("touch {}".format("check.log"))
    myfile = Path(check_log)
    myfile.touch(exist_ok=True)
    f = open(check_log, 'a+')
    t = str(datetime.datetime.now()).split(".")[0]
    f.write(t)
    f.write("\n")
    f.close()

    read_log_ = time.time()
    f = open(check_log, 'r+')
    d = f.readlines()
    print("check log :{}\n".format(d))
    f.close()

    for s in d:
        print("{}\n".format(s))
    finish = list(filter(lambda x: len(x.split(" ")) >= 4, d))
    finish_run = list(map(lambda x: x.split(" ")[1], finish))
    need_run = list(filter(lambda x: x not in finish_run, run_list))
    print("finish: {}\nfinish_run: {}\nneed_run".format(finish, finish_run, need_run))
    print("finish length: {}\nfinish_run length: {}\nneed_run length: ".format(len(finish), len(finish_run),
                                                                               len(need_run)))
    print("Toal", len(need_run), "sra runs need to downlaod.")


    ##all assembled
    if len(need_run) == 0:
        # utils_.progress_bar("ALL is assembled.")
        shutil.rmtree(sra_dir)
        shutil.rmtree(fastq_dir)
        shutil.rmtree(assemble_dir)
        print("ALL is assembled.\n")


    # ./YYYYMMDD/Assembled
    # assembly_date.py
    assem_cmd = "python3 assembly_datev3.py --PDAT {} --output {}".format(date, outdir)
    print("\n\n\n")
    print(assem_cmd)
    n=3
    k = list(range(0, len(need_run), 3))
    print(k)
    num = len(finish_run)

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
            outdir__ = os.path.join(outdir, "Assembled")

            final_dir = os.path.join(outdir__, "{}_contig.fa".format(x))
            ##Download and Assembled
            if os.path.isfile(final_dir):
                print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
            else:
                #Download
                utils_.prefetch_sra(x,sra_dir)
                #Assembled
                utils_.run_for_114v2(x,sra_dir,fastq_dir,assemble_dir,outdir,thread,gsize,start,check_log)
                current_path = os.path.join(os.path.abspath(os.getcwd()), x)
                print("current_path: ", current_path, "\n")
                # print ("shutil.rmtree({})\n".format(current_path))
                utils_.run_cmd2("rm -rf {}".format(current_path))
                print ("remove {}\n".format(current_path))
            ########################
            ##QualityCheck
            print("**********************************  QUALITYCHECKED  **********************************\n")
            target_path = os.path.join(outdir, "target.txt")
            check_file = os.path.join(outdir, "QCcheck.log")
            gnum = 0
            genome=os.path.join(outdir__,"{}_cotig.fa".format(x))
            if os.path.isfile(check_file):
                with open(check_file, "r") as f:
                    finishList = f.readlines()
                print(finishList)
                gnum = len(finishList)
                print("gnum={}".format(gnum))
            if gnum != count:
                print("Now QualityCheck excutting------------\n")

            print("**********************************   {} / {}   **********************************\n".format(gnum + 1,
                                                                                                               count))
            print(x)
            genome_path=os.path.join(outdir__,"{}_contig.fa".format(x))
            # python3 QualityCheck.py -r /data/usrhome/LabSSLin/user30/Desktop/RefSeq/ -g /data/usrhome/LabSSLin/user30/Desktop/SRA/test0812/assembly_result/contigs.fa -db enterobacterales_odb10 -m geno -o /data/usrhome/LabSSLin/user30/Desktop/QualityCheck
            qual_cmd = "python3 QualityCheckv3-124.py -r {} -g {} -db {} -m {} -o {}".format(ref_dir, genome_path, buscoDB, buscoMode,
                                                                                         outdir)
            print("run cmd: {}\n".format(qual_cmd))

            try:
                targetPath = utils_.run_cmd3(qual_cmd)
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
            gnum+=1
            print("targetPAth = {}".format(targetPath))
            print("**********************************  QUALITYCHECKED END  **********************************\n")

            with open(target_path, "r") as f:
                tlines = f.readlines()
                print(tlines)

            ########### analysis
            print("**********************************  Analysis  **********************************\n")
            Anacheck = os.path.join(outdir, "Anacheck.log")
            anum = 0
            if os.path.isfile(Anacheck):
                with open(Anacheck, "r") as f:
                    acheck = f.readlines()
                    print(acheck)
                    anum = len(acheck) - 1

            for target in tlines[anum:]:
                print("**********************************   {} / {}   **********************************\n".format(
                    anum + 1,
                    len(tlines)))
                target = target.split(":")[0]
                target_ = target.replace(current_path, ".")
                ana_cmd = "python3 analysisv4.py -i {} -o {} -mlstS {} -amrS {}".format(target_, outdir, mlstS, amrS)
                print(ana_cmd)
                utils_.run_cmd3(ana_cmd)
                anum += 1
            print("**********************************  ANA  End**********************************\n")
            print("Analysis Done.\n")
            f = open(check_log, 'a')
            f.write("Run {} is ok\n".format(x))
            f.close()
            print("Run {} is ok\n".format(x))



        print("shutil.rmtree(sra_dir)\n")
        shutil.rmtree(sra_dir)
        utils_.mkdir_join(sra_dir)





    print("**********************************  ASSEMBLED  END  **********************************\n")
    print("Assemble Done.\n")
    print("********************************************************************\n")
    sys.exit()


if __name__ == '__main__':
    main()