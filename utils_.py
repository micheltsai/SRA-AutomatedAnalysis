from __future__ import print_function
# -*- coding: utf-8 -*-
import csv
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
from Bio import Entrez
import datetime
import glob
import xml.etree.cElementTree as ET
from io import StringIO
from tempfile import TemporaryDirectory
Entrez.email = 'ann850324@gmail.com'

def progress_bar(Category):
    for i in range(1, 101):
        print("\r{}: ".format(Category),end="")
        print("[{}] {}%".format("*" * (i // 2), i), end="")
        sys.stdout.flush()
        time.sleep(0.02)
    print ("\n")

def mkdir_join(dir):
    try:
        os.makedirs(dir)
    except FileExistsError:
        print("folder is exist")
    except Exception as e:
        print(e)
    return dir

def prefetch_sra(sralist,outdir):
    #ss = " ".join(sralist)
    #progress_bar("prefetch")
    outdir
    try:
        cmd = "prefetch "+sralist+" --output-directory "+outdir
        print (cmd,"\n")
        run_cmd2(cmd)
        time.sleep(1)
    except Exception as e:
        print ("prefetch has problem:\n")
        error_class = e.__class__.__name__  # 取得錯誤類型
        detail = e.args[0]  # 取得詳細內容
        cl, exc, tb = sys.exc_info()  # 取得Call Stack
        lastCallStack = traceback.extract_tb(tb)[-1]  # 取得Call Stack的最後一筆資料
        fileName = lastCallStack[0]  # 取得發生的檔案名稱
        lineNum = lastCallStack[1]  # 取得發生的行號
        funcName = lastCallStack[2]  # 取得發生的函數名稱
        errMsg = "File \"{}\", line {}, in {}: [{}] {}".format(fileName, lineNum, funcName, error_class, detail)
        print(errMsg)

        #處理檔案已存在問題
        current_path = os.path.join(os.path.abspath(os.getcwd()), sralist)
        print(current_path, "\n")
        print("shutil.rmtree({})\n".format(current_path))
        shutil.rmtree(current_path)

        print("######run again\n")
        run_cmd2(cmd)
        time.sleep(1)
        pass
    print("now download", sralist, "runs.")

def run_cmd2(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True)
    return p

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

def run_cmd3(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    returncode = p.poll()
    while returncode is None:
        line = p.stdout.readline()
        returncode = p.poll()
        line = line.strip()
        line_ = line.decode().split("\n")
        #print (line.decode(),"\n\n")
        for s in line_:
            print(str("{}".format(s)))
            err=s
    return p,err

#assemble
###############################
# 獲得分布和質量並生成數據?
def bases_percentage(filepath, qscore=0):
    print("####################\nbases_percentage\n")
    cmd=f"seqtk fqchk -q {qscore} {filepath} | grep ALL | awk '{{print $NF}}'"
    print ("cmd: ",cmd)
    p = run_cmd2(cmd)
    print ("bases: ",p.stdout.decode())
    return float(p.stdout)


CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
ADAPTERS = os.path.join(CURRENT_DIR, 'trimmomatic.fa')
MIN_BQ = 3


# 把DNA序列兩端定序結果比較差的序列去掉
def crop_position(filepath, window_size=3, gap=10):
    print("#################### crop_postion ####################\n")
    p = run_cmd2("seqtk fqchk {}".format(filepath))
    #progress_bar("crop_postion")
    print("crop_postion\n")
    fq_check = p.stdout.decode().strip().split('\n')[3:]
    fq_check = (line.split()[2:6] for line in fq_check)
    content_gaps = []
    for line in fq_check:
        a, c, g, t = [float(value) for value in line]
        content_gaps.append(max(abs(a - t), abs(c - g)))
    # check from forward
    for start in range(len(content_gaps)):
        # 每四個做一次判斷
        # content_graps[start, end]
        end = start + window_size
        window = content_gaps[start: end]
        if max(window) < gap:
            headcrop = start
            break
        else:
            headcrop = 0
    # check from revers
    for start in range(len(content_gaps), 0, -1):
        end = start - window_size
        window = content_gaps[end:start]
        if max(window) < 10:
            crop = start
            break
        else:
            crop = len(content_gaps)
    print("#################### crop_postion END ####################\n")
    return crop, headcrop


def trimming(forward_reads, reverse_reads, outdir, threads):
    print("#################### trimming ####################\n")
    crop, headcrop = crop_position(forward_reads)
    opt = f"CROP:{crop} HEADCROP:{headcrop} ILLUMINACLIP:{ADAPTERS}:2:30:10 LEADING:{MIN_BQ} TRAILING:{MIN_BQ} " \
          f"SLIDINGWINDOW:4:20 MINLEN:36 TOPHRED33"
    paired_1 = os.path.join(outdir, 'R1.fq')
    paired_2 = os.path.join(outdir, 'R2.fq')
    cmd = f"java -jar ./trimmomatic-0.39.jar PE -threads {threads} {forward_reads} {reverse_reads} {paired_1} /dev/null" \
          f" {paired_2} /dev/null {opt}"
    run_cmd2(cmd)
    #progress_bar("trimming")
    print("#################### trimming END ####################\n")
    return paired_1, paired_2

def trimmingv2(forward_reads, reverse_reads, outdir, threads):
    print("#################### trimming ####################\n")
    crop, headcrop = crop_position(forward_reads)
    paired_1 = os.path.join(outdir, 'R1.fq')
    paired_2 = os.path.join(outdir, 'R2.fq')
    cmd=f"/data/usrhome/LabSSLin/user30/Desktop/fastp -i {forward_reads} -I {reverse_reads} -o {paired_1} -O {paired_2} --length_required 36 --detect_adapter_for_pe --cut_front 3 --cut_tail 3 -w {threads} -j /dev/null -h /dev/null"

    run_cmd2(cmd)
    #progress_bar("trimming")
    print("#################### trimming END ####################\n")
    return paired_1, paired_2


# sra轉換成fastq
def dump_fastq_from_sra(srafile, outdir):
    #progress_bar("dump_fastq_from_sra")
    print("#################### dump_fastq_from_sra ####################\n")
    run_cmd(f'fastq-dump --split-files --outdir {outdir} {srafile}')
    print("#################### dump_fastq_from_sra END ####################\n")


class SequenceReadArchivev2:
    #sra_file-->filepath
    def __init__(self,sraid):
       self._set_sraid(sraid)
       self._get_stat()

    #設置filepath
    def _set_sraid(self,sraid):
        self._sraid = sraid

    #發出指令
    def _get_stat(self):
        #sra-stat統計sra文件
        p = run_cmd2(f'sra-stat -x -s -b 1 -e 2 {self._sraid}')
        #xml: ElementTree, xml節點: Element
        #fromstring: 從xml_str構成Element賦予變數self._stat_tree
        self._stat_tree = ET.fromstring(p.stdout.decode())

    @property
    def sraid(self):
        return self._sraid

    @property
    def layout(self):
        #find()從節點的直接子節點中查詢,非遞回
        #.attrib: class dict
        return self._stat_tree.find('Statistics').attrib['nreads']

    def base_percentage(self):
        root = self._stat_tree.find('QualityCount')
        no30_count=0
        q30_count=0
        for child in root:
            #print(child.tag,":",child.attrib['value'].strip())
            q30=int(child.attrib['value'].strip())
            count=int(child.attrib['count'].strip())
            if q30>=30:
                q30_count+=count
            else:
                no30_count+=count
        print(q30_count , " : ", no30_count)
        return q30_count/(q30_count+no30_count)

class SequenceReadArchive:
    #sra_file-->filepath
    def __init__(self,filepath):
       self._set_filepath(filepath)
       self._get_stat()

    #設置filepath
    def _set_filepath(self,filepath):
        if os.access(filepath, os.F_OK) is False:
           raise FileNotFoundError("File not found.")
        with open(filepath,'rb') as handle:
           #if handle.read(8).decode() == 'NCBI.sra' is False:
           if handle.read(8).decode() != 'NCBI.sra':
              raise Exception(f"File format is not 'NCBI.sra'.")
        self._filepath = filepath

    #發出指令
    def _get_stat(self):
        #sra-stat統計sra文件
        p = run_cmd2(f'sra-stat -x -s -b 1 -e 2 {self._filepath}')
        #xml: ElementTree, xml節點: Element
        #fromstring: 從xml_str構成Element賦予變數self._stat_tree
        self._stat_tree = ET.fromstring(p.stdout.decode())

    @property
    def filepath(self):
        return self._filepath

    @property
    def layout(self):
        #find()從節點的直接子節點中查詢,非遞回
        #.attrib: class dict
        return self._stat_tree.find('Statistics').attrib['nreads']


def count_egquery(term, date_from , date_to, db = 'sra'):
    print("#################### count_egquery ####################\n")
    pattern = term+f" AND {date_from}[PDAT]:{date_to}[PDAT]"
    #progress_bar("count_egquery")
    print('Searching pattern:',pattern)
    #查詢PubMed中所有和pattern變數內容相關的文章
    handle = Entrez.egquery(term = pattern)
    d = Entrez.read(handle)
    #查詢d中的標籤`eGQueryResult`,  將DbName改成Db(`sra`)
    query = list(filter(lambda x: x['DbName'] == db, d['eGQueryResult']))[0]
    print('Total',query['Count'],'results in NCBI',db,'database.')
    print("#################### count_egquery END ####################\n")
    return pattern, query['Count']
    #term<---main.pattern
    #db='sra'
    #count=query['Count']

###### new add ######
def IdList_esearch(term, db, count):
    handle = Entrez.esearch(term = term, db = db, retmax = count)#不設定retmax的話只有20筆資料
    #progress_bar("read and stored IdList")
    d = Entrez.read(handle)
    return d['IdList']

def Get_RunInfo(idlist):
    print('Getting run_info table by idlist...')
    c = len(idlist)
    #progress_bar("read and stored RunInfo")
    if c >= 10000:
        print("over 10000 results")
        df_all = pd.DataFrame()
        #設定k為list,由0開始,最大值為c, 間隔10000
        k = list(range(0,c,10000))
        k.append(len(idlist))
        print(k)
        for i in range(1,len(k)):
            s = time.time()
            start = k[i-1]
            end = k[i]
            #下載GenBank records
            handle = Entrez.efetch(db = 'sra', id = idlist[start:end],rettype = 'runinfo',retmode = 'csv')
            #查看原始的Genbank文件
            d = handle.read()
            #讀檔
            df = pd.read_csv(StringIO(d))
            df = df[df['Run'] != 'Run']
            df_all = pd.concat([df_all,df])
            print("finish",i-1,"round,cost",time.time()-s,"secs")
            print("counts_of_run:",len(df))
            time.sleep(2)
    else:
        progress_bar("Entrez.efetch")
        handle = Entrez.efetch(db = 'sra', id = idlist,rettype = 'runinfo',retmode = 'csv')
        d = handle.read()
        df = pd.read_csv(StringIO(d))
        df_all = df[df['Run'] != 'Run']
    return df_all

def getlayout(path_):
    seq_readArchive = time.time()
    try:
        print("SequenceReadArchive\n")
        sra = SequenceReadArchive(path_)
        # print("layout:", sra.layout)
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
    if sra.layout != '2':
        sys.exit(f'File layout is not pair-end')
    print("layout=2\n")
    # if sra_layout==2 continue


def run_for_114v2(sra_id,sra_dir,fastq_dir,assemble_dir,outdir,threads,gsize,start,check_log):
    print ("sra_id = {}\nsra_dir = {}\noutdir= {}\n".format(sra_id,sra_dir,outdir))
    path_ = os.path.join(sra_dir,sra_id)
    path_=path_+"/"+sra_id+".sra"
    #outdir__=os.path.join(outdir, "Assembled")
    #mkdir_join(outdir__)
    #path_= os.path.join(path_,str("{}.sra".format(sra_id)))
    print ("srafile_path: {}\n".format(path_))
    # os.path.join(path, *paths)連接路徑

    #outdir = assem_dir + "/" + "".join(sra_id)
    print ("outdir = {}".format(outdir))
    #fastq_dir = os.path.join(outdir, 'fastq')
    fastq_dir_ = os.path.join(fastq_dir, sra_id)
    os.makedirs(fastq_dir_, exist_ok=True)
    print ("fastq_dir = {}".format(fastq_dir_))

    #assemble_dir = os.path.join(outdir, "assembly_result")
    assemble_dir_ = os.path.join(assemble_dir, sra_id)
    mkdir_join(assemble_dir_)
    contig_tmp = os.path.join(assemble_dir_, "contigs.fa")
    outdir__ = os.path.join(outdir, "Assembled")
    mkdir_join(outdir__)
    final_dir = os.path.join(outdir__, "{}_contig.fa".format(sra_id))
    #如果做過則下一個
    if os.path.isfile(final_dir):
        print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
        return 0

    dump_time=time.time()
    # 解壓縮成fastq
    print('Dump fastq.')
    # run_cmd
    dump_fastq_from_sra(path_, fastq_dir_)
    # os.listdir(fastq_dir) list files in dir
    print (fastq_dir_)

    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "dump_fastq_from_sra", "time": str(time.time() - dump_time)})

    reverse_time=time.time()
    try:
        forward_reads, reverse_reads = [os.path.join(fastq_dir_, fa) for fa in os.listdir(fastq_dir_)]
    except ValueError as e:
        if os.path.isfile(final_dir):
            print ("was ran assembly ,r1 and r2 is exist\n------------------------------\n\n")
            return 0
        else:
            #run_cmd("rm {}/R1.fq {}/R2.fq".format(fastq_dir,fastq_dir))
            run_cmd("rm -r {}".format(fastq_dir_))
            #agian
            run_for_114(sra_id,sra_dir,fastq_dir,assemble_dir,outdir,threads,gsize,start,check_log)
            #forward_reads, reverse_reads = [os.path.join(fastq_dir, fa) for fa in os.listdir(fastq_dir)]
            pass
    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "get forward and reverse reads", "time": str(time.time() - reverse_time)})

    ## up ok
    # 資料前處理：刪除爛的序列
    # Trimming sequence (trimmomatic)------- Q30 base >= 90% -----------> 預測基因組大小與定序深度(KMC & seqtk)

    #print('Trim sequences.')
    trim_time=time.time()
    r1, r2 = trimmingv2(forward_reads, reverse_reads, fastq_dir_, threads)
    print ("r1= {}, r2={}".format(r1,r2))
    # Q30>=90
    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "trimming", "time": str(time.time() - trim_time)})

    #bases_percentage_time=time.time()
    #if bases_percentage(r1, 30) < 90 and bases_percentage(r2, 30) < 90:
    #    #shutil.rmtree(outdir)
    #    sys.exit('Reads quality is too low.')

    #with open("./ana_time.csv", "a+") as f:
    #    fieldnames = ["func", "time"]
    #    writer = csv.DictWriter(f, fieldnames=fieldnames)
    #    writer.writeheader()
    #    writer.writerow({"func": "bases_percentage", "time": str(time.time() - bases_percentage_time)})

    # 預測基因組大小與定序深度(KMC & seqtk)--- depth>=80 ----> 抽樣(seqtk) -----------> SPAdes
    #                                 --- depth<80 ---------> SPAdes
    # de-novo assembly(SPAdes)------>Polish(pilon)---->Contings(最後成果檔案:conting.fa)
    #print("Run assembly pipline 'shovill'")
    progress_bar("Run assembly pipline 'shovill'")
    # depth >= 80

    shovill_time=time.time()
    #cmd = f"shovill --R1 {r1} --R2 {r2} --outdir {assemble_dir_} --depth 100 --tmpdir . --cpus {threads} --ram 3 --force"
    cmd = f"shovill --R1 {r1} --R2 {r2} --outdir {assemble_dir_} --depth 80 --tmpdir . --cpus {threads} --ram 3 --force"
    if gsize:
        cmd += f" --gsize {gsize}"
    print(cmd)
    run_cmd(cmd)

    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "shovill", "time": str(time.time() - shovill_time)})


    #cmd2 = "mv " + contig_tmp + " " + assemble_dir + "/" + sra_id + "_contig.fa && mv " + assemble_dir + "/" + sra_id + "_contig.fa " + outdir

    cmd2="cp {} {}".format(contig_tmp,final_dir)
    print("contig_tmp: {}\nfinal_dir: {}\ncmd2={}\n".format(contig_tmp,final_dir,cmd2))
    run_cmd(cmd2)
    #cmd3 = "rm -rf {}".format(assemble_dir)
    #run_cmd(cmd3)
    #f=open(check_log,"a")
    #f.write("Run {} is ok\n".format(sra_id))
    #f.close()


    shutil.rmtree(fastq_dir_)
    progress_bar("remove fastq dir")
    shutil.rmtree(assemble_dir_)
    progress_bar("remove assemble dir")


#run_for_114(x,sra_dir,output,threads,gsize,start,check_log)
def run_for_114(sra_id,sra_dir,fastq_dir,assemble_dir,outdir,threads,gsize,start,check_log):
    start_114=time.time()
    print ("sra_id = {}\nsra_dir = {}\noutdir= {}\n".format(sra_id,sra_dir,outdir))
    path_ = os.path.join(sra_dir,sra_id)
    path_=path_+"/"+sra_id+".sra"
    #outdir__=os.path.join(outdir, "Assembled")
    #mkdir_join(outdir__)
    #path_= os.path.join(path_,str("{}.sra".format(sra_id)))
    print ("srafile_path: {}\n".format(path_))
    seq_readArchive=time.time()
    try:
        print ("SequenceReadArchive\n")
        sra = SequenceReadArchive(path_)
        #print("layout:", sra.layout)
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
        run_cmd("rm -rf {}".format(fastq_dir))
        return run_for_114(sra_id,sra_dir,fastq_dir,assemble_dir,outdir,threads,gsize,start,check_log)
    if sra.layout != '2':
        sys.exit(f'File layout is not pair-end')
    print ("layout=2\n")
    # if sra_layout==2 continue
    # os.path.join(path, *paths)連接路徑

    #outdir = assem_dir + "/" + "".join(sra_id)
    print ("outdir = {}".format(outdir))
    #fastq_dir = os.path.join(outdir, 'fastq')
    fastq_dir_ = os.path.join(fastq_dir, sra_id)
    os.makedirs(fastq_dir_, exist_ok=True)
    print ("fastq_dir = {}".format(fastq_dir_))

    #assemble_dir = os.path.join(outdir, "assembly_result")
    assemble_dir_ = os.path.join(assemble_dir, sra_id)
    mkdir_join(assemble_dir_)
    contig_tmp = os.path.join(assemble_dir_, "contigs.fa")
    outdir__ = os.path.join(outdir, "Assembled")
    mkdir_join(outdir__)
    final_dir = os.path.join(outdir__, "{}_contig.fa".format(sra_id))
    #如果做過則下一個
    if os.path.isfile(final_dir):
        print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
        return 0

    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "SequenceReadArchive", "time": str(time.time() - seq_readArchive)})
    dump_time=time.time()
    # 解壓縮成fastq
    print('Dump fastq.')
    # run_cmd
    dump_fastq_from_sra(path_, fastq_dir_)
    # os.listdir(fastq_dir) list files in dir
    print (fastq_dir_)

    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "dump_fastq_from_sra", "time": str(time.time() - dump_time)})

    reverse_time=time.time()
    try:
        forward_reads, reverse_reads = [os.path.join(fastq_dir_, fa) for fa in os.listdir(fastq_dir_)]
    except ValueError as e:
        if os.path.isfile(final_dir):
            print ("was ran assembly ,r1 and r2 is exist\n------------------------------\n\n")
            return 0
        else:
            #run_cmd("rm {}/R1.fq {}/R2.fq".format(fastq_dir,fastq_dir))
            run_cmd("rm -r {}".format(fastq_dir_))
            #agian
            run_for_114(sra_id,sra_dir,fastq_dir,assemble_dir,outdir,threads,gsize,start,check_log)
            print("remove fastq\n resart run_for_114 to Assembled\ncommit to error.txt\n")
            with open("./SRA_run_error.txt", "a+") as f:
                f.write("{} :\n{}\n".format(sra_id,e))
            #forward_reads, reverse_reads = [os.path.join(fastq_dir, fa) for fa in os.listdir(fastq_dir)]
            exit()
    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "get forward and reverse reads", "time": str(time.time() - reverse_time)})

    ## up ok
    # 資料前處理：刪除爛的序列
    # Trimming sequence (trimmomatic)------- Q30 base >= 90% -----------> 預測基因組大小與定序深度(KMC & seqtk)

    #print('Trim sequences.')
    trim_time=time.time()
    r1, r2 = trimming(forward_reads, reverse_reads, fastq_dir_, threads)
    print ("r1= {}, r2={}".format(r1,r2))
    # Q30>=90
    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "trimming", "time": str(time.time() - trim_time)})

    bases_percentage_time=time.time()
    if bases_percentage(r1, 30) < 90 and bases_percentage(r2, 30) < 90:
        shutil.rmtree(outdir)
        sys.exit('Reads quality is too low.')

    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "bases_percentage", "time": str(time.time() - bases_percentage_time)})

    # 預測基因組大小與定序深度(KMC & seqtk)--- depth>=80 ----> 抽樣(seqtk) -----------> SPAdes
    #                                 --- depth<80 ---------> SPAdes
    # de-novo assembly(SPAdes)------>Polish(pilon)---->Contings(最後成果檔案:conting.fa)
    #print("Run assembly pipline 'shovill'")
    progress_bar("Run assembly pipline 'shovill'")
    # depth >= 80

    shovill_time=time.time()
    cmd = f"shovill --R1 {r1} --R2 {r2} --outdir {assemble_dir_} --depth 100 --tmpdir . --cpus {threads} --ram 3 --force"
    if gsize:
        cmd += f" --gsize {gsize}"
    print(cmd)
    run_cmd(cmd)

    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "shovill", "time": str(time.time() - shovill_time)})


    #cmd2 = "mv " + contig_tmp + " " + assemble_dir + "/" + sra_id + "_contig.fa && mv " + assemble_dir + "/" + sra_id + "_contig.fa " + outdir

    cmd2="cp {} {}".format(contig_tmp,final_dir)
    print("contig_tmp: {}\nfinal_dir: {}\ncmd2={}\n".format(contig_tmp,final_dir,cmd2))
    run_cmd(cmd2)
    #cmd3 = "rm -rf {}".format(assemble_dir)
    #run_cmd(cmd3)
    #f=open(check_log,"a")
    #f.write("Run {} is ok\n".format(sra_id))
    #f.close()
    f = open(check_log, 'a')
    f.write("Run {} is ok\n".format(sra_id))
    f.close()
    print ("Assembled Run {} is ok\n".format(sra_id))
    print('{} Done,total cost'.format(sra_id), time.time() - start_114, 'secs')

    shutil.rmtree(fastq_dir_)
    #progress_bar("remove fastq dir")
    print("shutil.rmtree(fastq_dir_)\n")
    shutil.rmtree(assemble_dir_)
    print("shutil.rmtree(assemble_dir_)\n")
    #progress_bar("remove assemble dir")

#Qualitycheck.py
##################################
#get reference list and path of reference list
def getRefListPath(refSeqPath,outdir):
    print("getRefListPath:\n")
    print("refSeqPath: " + refSeqPath + "\n")
    refListPath = os.path.join(outdir, 'ref.txt')
    if os.path.isfile(refListPath):
        run_cmd2("rm -f {}".format(refListPath))
        #shutil.rmtree(refListPath)
    run_cmd2("find {} -type f >{}".format(refSeqPath,refListPath))
    print("refListPath: "+refListPath+"\n")
    #run_cmd2("cat {}".format(refListPath))
    return refListPath

#get qenome list and path of qenome list
def getGenomeListPath(genome_Path,outdir):
    print("getGenomeListPath:\n")
    print("genSeqPath: " + genome_Path + "\n")
    genListPath = os.path.join(outdir, 'Assembled.txt')
    if os.path.isfile(genListPath):
        run_cmd2("rm {}".format(genListPath))
    run_cmd2("find {} -type f >{}".format(genome_Path,genListPath))
    print("genListPath: "+genListPath+"\n")
    return genListPath



