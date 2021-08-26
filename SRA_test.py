import time
import os
import re
import sys
import shutil
import argparse
import subprocess
import traceback
import xml.etree.cElementTree as ET
from os.path import isfile
from tempfile import TemporaryDirectory
import pandas as pd
import datetime
import glob
from Bio import Entrez
from io import StringIO
Entrez.email = 'ann850324@gmail.com'

def count_egquery(term, date_from , date_to, db = 'sra'):
    print ("---------------count_egquery------------------\n")

    pattern = term+f" AND {date_from}[PDAT]:{date_to}[PDAT]"
    print('Searching pattern:',pattern)
    handle = Entrez.egquery(term = pattern)
    d = Entrez.read(handle)
    query = list(filter(lambda x: x['DbName'] == db, d['eGQueryResult']))[0]
    print('Total',query['Count'],'results in NCBI',db,'database.')
    return pattern, query['Count']

def IdList_esearch(term, db, count):
    print ("---------------IdList_esearch------------------\n")
    handle = Entrez.esearch(term = term, db = db, retmax = count)#不設定retmax的話只有20筆資料
    d = Entrez.read(handle)
    return d['IdList']

def Get_RunInfo(idlist):
    print ("---------------Get_RunInfo------------------\n")
    print('Getting run_info table by idlist...')
    c = len(idlist)
    if c >= 10000:
        print("over 10000 results")
        df_all = pd.DataFrame()
        k = list(range(0,c,10000))
        k.append(len(idlist))
        print(k)
        for i in range(1,len(k)):
            s = time.time()
            start = k[i-1]
            end = k[i]
            handle = Entrez.efetch(db = 'sra', id = idlist[start:end],rettype = 'runinfo',retmode = 'csv')
            d = handle.read()
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

###
def run_cmd(cmd):
    print ("---------------run_cmd------------------\n")
    print ("{}\n".format(cmd))
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True)
    return p

def run_cmd2(cmd):
    print ("---------------run_cmd2------------------\n")
    print ("{}\n".format(cmd))
    #p = subprocess.run(cmd, stderr=subprocess.PIPE, shell=True, check=True)
    p = subprocess.run(cmd, shell=True, check=True)
    return p

###for_assembly
def bases_percentage(filepath, qscore=0):
    print ("---------------bases_percentage------------------\n")
    p = run_cmd(f"seqtk fqchk -q {qscore} {filepath} | grep ALL | awk '{{print $NF}}'")
    return float(p.stdout)

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
ADAPTERS = os.path.join(CURRENT_DIR, 'trimmomatic.fa')
MIN_BQ = 3

def crop_position(filepath, window_size=3, gap=10):
    print ("---------------crop_position------------------\n")
    p = run_cmd(f"seqtk fqchk {filepath}")
    fq_check = p.stdout.decode().strip().split('\n')[3:]
    fq_check = (line.split()[2:6] for line in fq_check)
    content_gaps = []
    for line in fq_check:
        a, c, g, t = [float(value) for value in line]
        content_gaps.append(max(abs(a - t), abs(c - g)))
    # check from forward
    for start in range(len(content_gaps)):
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
    print ("---------------trimming------------------\n")
    crop, headcrop = crop_position(forward_reads)
    opt = f"CROP:{crop} HEADCROP:{headcrop} ILLUMINACLIP:{ADAPTERS}:2:30:10 LEADING:{MIN_BQ} TRAILING:{MIN_BQ} " \
          f"SLIDINGWINDOW:4:20 MINLEN:36 TOPHRED33"
    paired_1 = os.path.join(outdir, 'R1.fq')
    paired_2 = os.path.join(outdir, 'R2.fq')
    cmd = f"java -jar trimmomatic-0.39.jar PE -threads {threads} {forward_reads} {reverse_reads} {paired_1} /dev/null" \
          f" {paired_2} /dev/null {opt}"
    run_cmd(cmd)
    return paired_1, paired_2

def dump_fastq_from_sra(srafile, outdir):
    print ("---------------dump_fastq_from_sra------------------\n")
    run_cmd(f'fastq-dump --split-files --outdir {outdir} {srafile}')


class SequenceReadArchive:
    def __init__(self,filepath):
       self._set_filepath(filepath)
       self._get_stat()

    def _set_filepath(self,filepath):
        if os.access(filepath, os.F_OK) is False:
           raise FileNotFoundError("File not found.")
        with open(filepath,'rb') as handle:
           if handle.read(8).decode() == 'NCBI.sra' is False:
              raise Exception(f"File format is not 'NCBI.sra'.")
        self._filepath = filepath

    def _get_stat(self):
        p = run_cmd(f'sra-stat -x -s -b 1 -e 2 {self._filepath}')
        self._stat_tree = ET.fromstring(p.stdout.decode())

    @property
    def filepath(self):
        return self._filepath

    @property
    def layout(self):
        return self._stat_tree.find('Statistics').attrib['nreads']
###

def prefetch_sra(sralist,outdir):
    print ("---------------prefetch_sra------------------\n")
    ss = " ".join(sralist)
    print("now download",ss,"runs.")
    cmd = "prefetch "+ss+" --output-directory "+outdir
    run_cmd2(cmd)

def run_dump_assembly(need_run,sra_dir,assem_dir,threads,gsize,output,check_log,n):
    print ("---------------run_dump_assembly------------------\n")
    k = list(range(0,len(need_run),n))
    for i in k:
        run_id = need_run[i:i+n]
        prefetch_sra(run_id,sra_dir)
        time.sleep(1)
        print ("run_id: {}\nsra_dir: {}\n".format(run_id,sra_dir))
        for x in run_id:
            start = time.time()

            while isfile(x):
                prefetch_sra(x, sra_dir)

            path = sra_dir+"/"+x+"/*.sra"
            re_path = "".join(glob.glob(path))
            sra_file = os.path.abspath(re_path)

            sra_file="{}/{}/{}.sra".format(sra_dir,x,x)
            print ("path: {}\nsra_file: {}\n".format(path,sra_file))
            try:
                sra = SequenceReadArchive(sra_file)
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
            print ("layout???---------------\n")
            if sra.layout != '2':
                f = open(check_log,"a")
                m = "Run "+x+" layout is not 2."
                f.write(m)
                f.close()
                #sys.exit(f'File layout is not pair-end')
                continue
            print ("layout=2---------------\n")
            print('Now assembly',sra_file,'...')
            outdir = assem_dir+"/"+"".join(x)
            fastq_dir = os.path.join(outdir,'fastq')
            print ("outdir: {}\nfastq_dir: {}\n".format(outdir,fastq_dir))
            os.makedirs(fastq_dir, exist_ok = True)
            #print('Dump fastq.')
            dump_fastq_from_sra(sra_file, fastq_dir)
            try:
                forward_reads, reverse_reads = [os.path.join(fastq_dir,i) for i in os.listdir(fastq_dir)]
            except Exception as e:
                continue
            print ("forward_reads: {}\nreverse_reads: {}\n".format(forward_reads,reverse_reads))
            #print('Trim sequences.')
            r1, r2 = trimming(forward_reads, reverse_reads, fastq_dir, threads)
            print ("r1: {}\nr2: {}\n".format(r1,r2))

            if bases_percentage(r1, 30) < 90 and bases_percentage(r2, 30) < 90:
                f = open(check_log,"a")
                m = "Run "+x+" reads quality is too low."
                f.write(m)
                f.close()
                #sys.exit('Reads quality is too low.')
                continue
            print("Run assembly pipline 'shovill'")
            assemble_dir = os.path.join(outdir,'assembly_result')
            print ("assemble_dir: {}\n".format(assemble_dir))
            cmd = f"shovill --R1 {r1} --R2 {r2} --outdir {assemble_dir} --depth 100 --tmpdir . --cpus {threads} --ram 3 --force"
            print ("cmd: {}\n".format(cmd))
            if gsize:
                cmd += f" --gsize {gsize}"
            print ("+gsize cmd: {}\n".format(cmd))
            run_cmd(cmd)
            print('Done.')
            contig_check = assemble_dir+"/contigs.fa"
            print ("contig_check: {}\n".format(contig_check))
            if os.path.exists(contig_check) == False:
                f = open(check_log,"a")
                m = "Run "+x+" have assembly error. Didn't have Contigs.fa file.\n"
                f.write(m)
                f.close()
                continue
            print ("os.path.exists(contig_check) == False\n")
            f = open(check_log,"a")
            m = "Run "+x+" is ok.\n"
            f.write(m)
            f.close()
            cmd = "mv "+contig_check+" "+assemble_dir+"/"+x+"_contig.fa && mv "+assemble_dir+"/"+x+"_contig.fa "+output
            print ("cmd: {}\n".format(cmd))
            run_cmd(cmd)
            print ("shutil.rmtree\n")
            shutil.rmtree(outdir)               
        shutil.rmtree(sra_dir)
        os.mkdir(sra_dir)

def main():
    parser = argparse.ArgumentParser("Download_avaliable_runs_from_NCBI_sra_database_wit_PDAT.") 
    parser.add_argument("--pattern",default = "salmonella enterica[ORGN] AND illumina[PLAT] AND wgs[STRA] AND genomic[SRC] AND paired[LAY]", help = "Searching condition.")
    parser.add_argument("--PDAT",required = True, help = "Publication Date[PDAT] of Runs.")
    parser.add_argument("--sra_dir",required = True, help ="Temp folder to save runs.")
    parser.add_argument("--log",required = True,help = "The file to recored runs which finished assembly each time.")
    parser.add_argument("--output",required = True, help = "Folder to save Contigs.fa after assembly of eah runs.")
    parser.add_argument("--assembly_dir", default = ".", help = "Temp folder to save assembly.")
    parser.add_argument("--tmpdir", default="/tmp", help="Directory of temp folder default: '/tmp'")
    parser.add_argument("--gsize", default='', help="Estimated genome size(MB) eg. 3.2M. default: ''")
    parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    parser.add_argument("--n", default = 3,help="count of download sra file eahc time", type=int)
    args = parser.parse_args()
    
    pattern = args.pattern
    #print(pattern)    
    date = args.PDAT
    print("date: {}\npattern: {}".format(date, pattern))
    #print(date)
    pattern, count = count_egquery(pattern,date,date)
    #print(pattern,count)
    idlist = IdList_esearch(pattern,'sra',count)
    print("after pattern: {}\ncount: {}\nidlist: {}".format(pattern,count,idlist))
    runinfo = Get_RunInfo(idlist)
    run_list = list(runinfo['Run'])
    sra_dir = args.sra_dir
    assem_dir = args.assembly_dir
    check_log = args.log
    output = args.output
    tmpdir = args.tmpdir
    threads = args.threads
    gsize = args.gsize
    print("runinfo : {}\nrun_list: {}\nsra_dir: {}\nassem_dir: {}\ncheck_log: {}\noutput: {}\ntmpdir: {}\nthreads: {}\ngsize: {}\n".format(runinfo, run_list, sra_dir,assem_dir,check_log,output,tmpdir,threads,gsize))
    n = args.n
    print ("n: {}\n".format(n))
    f = open(check_log)
    d = f.read().split("\n")[1:-1]
    f.close()
    print ("d: {}\n".format(d))
    finish = list(filter(lambda x:len(x.split(" ")) >= 4 ,d))
    finish_run = list(map(lambda x : x.split(" ")[1],finish))
    need_run = list(filter(lambda x : x not in finish_run, run_list))
    print("finish: {}\nfinish_run: {}\nneed_run".format(finish, finish_run, need_run))
    print("finish length: {}\nfinish_run length: {}\nneed_run length: ".format(len(finish), len(finish_run),
                                                                               len(need_run)))
    print("Toal",len(need_run),"sra runs need to downlaod.")
    os.makedirs(sra_dir, exist_ok = True)
    os.makedirs(assem_dir,exist_ok = True)
    f = open(check_log,'a')
    t = str(datetime.datetime.now()).split(".")[0]
    print ("t: {}\n".format(t))
    f.write(t)
    f.write("\n")
    f.close()
    if os.path.exists(output) == False:
        print ("os.path.exists(output) == False??")
        os.makedirs(output)
    print ("run dump assembly---------\n")
    run_dump_assembly(need_run,sra_dir,assem_dir,threads,gsize,output,check_log,n)
   
if __name__ == '__main__':
   main()
