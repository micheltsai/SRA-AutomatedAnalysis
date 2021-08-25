from __future__ import print_function
# -*- coding: utf-8 -*-
import shlex
import time
import os
import re
import sys
import shutil
import argparse
import subprocess
import traceback

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
    return dir

def prefetch_sra(sralist,outdir):
    ss = " ".join(sralist)
    print("now download",ss,"runs.")
    cmd = "prefetch "+ss+" --output-directory "+outdir
    run_cmd2(cmd)

def run_cmd2(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True)
    return p

def run_cmd(cmd):
    cmd=shlex.split(cmd)
    p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #print (cmd)
    print("--------------------------------------\nSubprogram output:\n")
    while p.poll() is None:
        progress_bar("sub excuting")
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

# 獲得分布和質量並生成數據?
def bases_percentage(filepath, qscore=0):
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
    p = run_cmd2("seqtk fqchk {}".format(filepath))
    #
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
    return crop, headcrop


def trimming(forward_reads, reverse_reads, outdir, threads):
    crop, headcrop = crop_position(forward_reads)
    opt = f"CROP:{crop} HEADCROP:{headcrop} ILLUMINACLIP:{ADAPTERS}:2:30:10 LEADING:{MIN_BQ} TRAILING:{MIN_BQ} " \
          f"SLIDINGWINDOW:4:20 MINLEN:36 TOPHRED33"
    paired_1 = os.path.join(outdir, 'R1.fq')
    paired_2 = os.path.join(outdir, 'R2.fq')
    cmd = f"java -jar /data/usrhome/LabSSLin/user30/Desktop/SRA/trimmomatic-0.39.jar PE -threads {threads} {forward_reads} {reverse_reads} {paired_1} /dev/null" \
          f" {paired_2} /dev/null {opt}"
    run_cmd2(cmd)
    return paired_1, paired_2


# sra轉換成fastq
def dump_fastq_from_sra(srafile, outdir):
    run_cmd(f'fastq-dump --split-files --outdir {outdir} {srafile}')


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
           if handle.read(8).decode() == 'NCBI.sra' is False:
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
    pattern = term+f" AND {date_from}[PDAT]:{date_to}[PDAT]"
    print('Searching pattern:',pattern)
    #查詢PubMed中所有和pattern變數內容相關的文章
    handle = Entrez.egquery(term = pattern)
    d = Entrez.read(handle)
    #查詢d中的標籤`eGQueryResult`,  將DbName改成Db(`sra`)
    query = list(filter(lambda x: x['DbName'] == db, d['eGQueryResult']))[0]
    print('Total',query['Count'],'results in NCBI',db,'database.')
    return pattern, query['Count']
    #term<---main.pattern
    #db='sra'
    #count=query['Count']

###### new add ######
def IdList_esearch(term, db, count):
    handle = Entrez.esearch(term = term, db = db, retmax = count)#不設定retmax的話只有20筆資料
    progress_bar("read and stored IdList")
    d = Entrez.read(handle)
    return d['IdList']

def Get_RunInfo(idlist):
    print('Getting run_info table by idlist...')
    c = len(idlist)
    progress_bar("read and stored RunInfo")
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
        handle = Entrez.efetch(db = 'sra', id = idlist,rettype = 'runinfo',retmode = 'csv')
        d = handle.read()
        df = pd.read_csv(StringIO(d))
        df_all = df[df['Run'] != 'Run']
    return df_all



#run_for_114(x,sra_dir,output,threads,gsize,start,check_log)
def run_for_114(sra_id,sra_dir,outdir,threads,gsize,start,check_log):
    print ("sra_id = {}\nsra_dir = {}\noutdir= {}\n".format(sra_id,sra_dir,outdir))
    path_ = os.path.join(sra_dir,sra_id)
    path_=path_+"/"+sra_id+".sra"
    #path_= os.path.join(path_,str("{}.sra".format(sra_id)))
    print ("srafile_path: {}\n".format(path_))
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
        sys.exit(e)
    if sra.layout != '2':
        sys.exit(f'File layout is not pair-end')
    print ("layout=2\n")
    # if sra_layout==2 continue
    # os.path.join(path, *paths)連接路徑
    fastq_dir = os.path.join(outdir, 'fastq')
    os.makedirs(fastq_dir, exist_ok=True)
    print ("fastq_dir = {}".format(fastq_dir))
    # 解壓縮成fastq
    print('Dump fastq.')
    # run_cmd
    dump_fastq_from_sra(path_, fastq_dir)
    # os.listdir(fastq_dir) list files in dir
    forward_reads, reverse_reads = [os.path.join(fastq_dir, i) for i in os.listdir(fastq_dir)]
    ## up ok
    # 資料前處理：刪除爛的序列
    # Trimming sequence (trimmomatic)------- Q30 base >= 90% -----------> 預測基因組大小與定序深度(KMC & seqtk)

    print('Trim sequences.')
    r1, r2 = trimming(forward_reads, reverse_reads, fastq_dir, threads)
    print ("r1= {}, r2={}".format(r1,r2))
    # Q30>=90
    if bases_percentage(r1, 30) < 90 and bases_percentage(r2, 30) < 90:
        shutil.rmtree(outdir)
        sys.exit('Reads quality is too low.')

    # 預測基因組大小與定序深度(KMC & seqtk)--- depth>=80 ----> 抽樣(seqtk) -----------> SPAdes
    #                                 --- depth<80 ---------> SPAdes
    # de-novo assembly(SPAdes)------>Polish(pilon)---->Contings(最後成果檔案:conting.fa)
    print("Run assembly pipline 'shovill'")
    assemble_dir = os.path.join(outdir, 'assembly_result')
    # depth >= 80
    cmd = f"shovill --R1 {r1} --R2 {r2} --outdir {assemble_dir} --depth 100 --tmpdir . --cpus {threads} --ram 3 --force"
    if gsize:
        cmd += f" --gsize {gsize}"
    print(cmd)
    run_cmd(cmd)
    contig_tmp = os.path.join(assemble_dir,"/contig.fa")
    cmd2 = "mv " + contig_tmp + " " + assemble_dir + "/" + sra_id + "_contig.fa && mv " + assemble_dir + "/" + sra_id + "_contig.fa " + outdir
    run_cmd(cmd2)
    #f=open(check_log,"a")
    #f.write("Run {} is ok\n".format(sra_id))
    #f.close()
    print ("Run {} is ok\n".format(sra_id))

def main():
    # 計算時間
    start = time.time()
    progress_bar("read command")
    parser = argparse.ArgumentParser("Download_avaliable_runs_from_NCBI_sra_database_wit_PDAT.")
    parser.add_argument("--pattern",
                        default="salmonella enterica[ORGN] AND illumina[PLAT] AND wgs[STRA] AND genomic[SRC] AND paired[LAY]",
                        help="Searching condition.")
    # PDAT格式：YYYY/MM/DD
    parser.add_argument("--PDAT", required=True, help="Publication Date[PDAT] of Runs.")
    parser.add_argument("--sra_dir", required=True, help="Temp folder to save runs.")
    parser.add_argument("--log", required=True, help="The file to recored runs which finished assembly each time.")
    parser.add_argument("--output", required=True, help="Folder to save Contigs.fa after assembly of eah runs.")
    parser.add_argument("--assembly_dir", default=".", help="Temp folder to save assembly.")
    parser.add_argument("--tmpdir", default="/tmp", help="Directory of temp folder default: '/tmp'")
    parser.add_argument("--gsize", default='', help="Estimated genome size(MB) eg. 3.2M. default: ''")
    parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    parser.add_argument("--n", default=3, help="count of download sra file eahc time", type=int)
    args = parser.parse_args()

    pattern = args.pattern
    # print(pattern)
    date = args.PDAT
    sra_dir = args.sra_dir
    assem_dir = args.assembly_dir
    log_dir = args.log
    output = args.output
    tmpdir = args.tmpdir
    threads = args.threads
    gsize = args.gsize
    n = args.n
    progress_bar("read arguments")

    mkdir_join(log_dir)
    check_log = os.path.join(log_dir, "check.log")

    # print(date)
    pattern, count = count_egquery(pattern, date, date)
    print ("pattern: {}\ncount: {}\n".format(pattern,count))
    idlist = IdList_esearch(pattern, 'sra', count)
    print(idlist)
    runinfo = Get_RunInfo(idlist)
    progress_bar("get SRAfile name List stored in run_list")
    run_list = list(runinfo['Run']) #get SRAfile nameList stored in run_list
    print("runinfo: {}\n run_list: {}\n".format(runinfo, run_list))

    f = open(check_log, 'w+')
    d = f.read().split("\n")[1:-1]
    f.close()
    finish = list(filter(lambda x: len(x.split(" ")) >= 4, d))
    finish_run = list(map(lambda x: x.split(" ")[1], finish))
    need_run = list(filter(lambda x: x not in finish_run, run_list))
    print("finish: {}\nfinish_run: {}\nneed_run".format(finish,finish_run,need_run))
    print("finish length: {}\nfinish_run length: {}\nneed_run length: ".format(len(finish), len(finish_run), len(need_run)))

    print("Toal", len(need_run), "sra runs need to downlaod.")

    os.makedirs(sra_dir, exist_ok=True)
    os.makedirs(assem_dir, exist_ok=True)

    #f = open(check_log, 'a')
    #t = str(datetime.datetime.now()).split(".")[0]
    #f.write(t)
    #f.write("\n")
    #f.close()

    if os.path.exists(output) == False:
        os.makedirs(output)

    k = list(range(0, len(need_run), n))
    print (k)
    for i in k:
        run_id = need_run[i:i+n]
        prefetch_sra(run_id,sra_dir)
        print("###### i = {}\n".format(i))
        print("run_id: {}\n".format(run_id))
        time.sleep(1)
        for x in run_id:
            print ("x = {}".format(x))
            run_for_114(x,sra_dir,output,threads,gsize,start,check_log)





    print('Done,total cost', time.time() - start, 'secs')


if __name__ == '__main__':
    main()
