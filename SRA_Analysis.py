from __future__ import print_function

import argparse
import os
import sys
import time
import traceback

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

    # ./YYYYMMDD/Assembled
    # assembly_date.py
    assem_cmd = "python3 assembly_datev3.py --PDAT {} --output {}".format(date, outdir)
    print("\n\n\n")
    print(assem_cmd)
    print("**********************************  ASSEMBLED  **********************************\n")

    try:
        utils_.run_cmd3(assem_cmd)
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
    print("**********************************  ASSEMBLED  END  **********************************\n")
    print("Assemble Done.\n")
    print("********************************************************************\n")

    ####################
    # QualityCheck
    assembled_path = os.path.join(outdir, "Assembled/")
    assembled_file = utils_.getGenomeListPath(assembled_path, outdir)
    with open(assembled_file, "r") as f:
        genomes = f.readlines()
    target_path = os.path.join(outdir, "target.txt")
    check_file=os.path.join(outdir, "QCcheck.log")
    gnum=0

    if os.path.isfile(check_file):
        with open(check_file,"r") as f:
            finishList=f.readlines()
        print(finishList)
        gnum=len(finishList)
        print("gnum={}".format(gnum))
    if gnum != len(genomes):
        print("Now QualityCheck excutting------------\n")

    #gnum=1
    print("**********************************  QUALITYCHECKED  **********************************\n")
    for g in genomes[gnum:]:
        g = g.strip("\n")
        print("**********************************   {} / {}   **********************************\n".format(gnum+1, len(genomes)))
        print(g)
        # python3 QualityCheck.py -r /data/usrhome/LabSSLin/user30/Desktop/RefSeq/ -g /data/usrhome/LabSSLin/user30/Desktop/SRA/test0812/assembly_result/contigs.fa -db enterobacterales_odb10 -m geno -o /data/usrhome/LabSSLin/user30/Desktop/QualityCheck
        qual_cmd = "python3 QualityCheckv3.py -r {} -g {} -db {} -m {} -o {}".format(ref_dir, g, buscoDB, buscoMode, outdir)
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
        gnum += 1
        print("targetPAth = {}".format(targetPath))
    print("**********************************  QUALITYCHECKED END  **********************************\n")
    print("QualityCheck Done.\n")
    print("********************************************************************\n")


    with open(target_path , "r") as f:
        tlines=f.readlines()
        print(tlines)

    ########### analysis
    print("**********************************  Analysis  **********************************\n")
    Anacheck = os.path.join(outdir, "Anacheck.log")
    anum=0
    if os.path.isfile(Anacheck):
        with open(Anacheck, "r") as f:
            acheck=f.readlines()
            print(acheck)
            anum=len(acheck)-1

    for target in tlines[anum:]:
        print("**********************************   {} / {}   **********************************\n".format(anum + 1,
                                                                                                           len(tlines)))
        target=target.split(":")[0]
        target_=target.replace(current_path,".")
        ana_cmd="python3 analysisv4.py -i {} -o {} -mlstS {} -amrS {}".format(target_,outdir,mlstS,amrS)
        print(ana_cmd)
        utils_.run_cmd3(ana_cmd)
        anum+=1
    print("**********************************  ANA  End**********************************\n")
    print("Analysis Done.\n")
    return 0


if __name__ == '__main__':
    main()