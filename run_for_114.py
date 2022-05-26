#run_for_114.py
import time
import os
import re
import sys
import shutil
import argparse
import subprocess
import traceback
import xml.etree.cElementTree as ET
from tempfile import TemporaryDirectory

#於cmd執行程式,捕捉字串？
def run_cmd(cmd):
    print ("--------------------------------\ncmd: {} \n".format(cmd))
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True)
    return p

#獲得分布和質量並生成數據?
def bases_percentage(filepath, qscore=0):
    print ("---------------bases_prercentage------------------\n")
    p = run_cmd(f"seqtk fqchk -q {qscore} {filepath} | grep ALL | awk '{{print $NF}}'")
    return float(p.stdout)

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
ADAPTERS = os.path.join(CURRENT_DIR, 'trimmomatic.fa')
MIN_BQ = 3

#把DNA序列兩端定序結果比較差的序列去掉
def crop_position(filepath, window_size=3, gap=10):
    p = run_cmd(f"seqtk fqchk {filepath}")
    #
    fq_check = p.stdout.decode().strip().split('\n')[3:]
    print ("fq_check: {}\n".format(fq_check))
    fq_check = (line.split()[2:6] for line in fq_check)
    print ("fq_check: {}\n".format(fq_check))
    content_gaps = []
    for line in fq_check:
        a, c, g, t = [float(value) for value in line]
        content_gaps.append(max(abs(a - t), abs(c - g)))
    print ("content_gaps: {}\n".format(content_gaps))
    # check from forward
    for start in range(len(content_gaps)):
        #每四個做一次判斷
        #content_graps[start, end]
        end = start + window_size
        window = content_gaps[start: end]
        if max(window) < gap:
            headcrop = start
            break
        else:
            headcrop = 0
    print ("headcrop: {}\n".format(headcrop))
    # check from revers
    for start in range(len(content_gaps), 0, -1):
        end = start - window_size
        window = content_gaps[end:start]
        if max(window) < 10:
            crop = start
            break
        else:
            crop = len(content_gaps)
    print ("crop: {}\n".format(crop))
    return crop, headcrop

def trimming(forward_reads, reverse_reads, outdir, threads):
    crop, headcrop = crop_position(forward_reads)
    opt = f"CROP:{crop} HEADCROP:{headcrop} ILLUMINACLIP:{ADAPTERS}:2:30:10 LEADING:{MIN_BQ} TRAILING:{MIN_BQ} " \
          f"SLIDINGWINDOW:4:20 MINLEN:36 TOPHRED33"
    paired_1 = os.path.join(outdir, 'R1.fq')
    paired_2 = os.path.join(outdir, 'R2.fq')

    cmd = f"java -jar trimmomatic-0.39.jar PE -threads {threads} {forward_reads} {reverse_reads} {paired_1} /dev/null" \
          f" {paired_2} /dev/null {opt}"
    print ("crop: {}\nheadcrop: {}\nR1: {}\nR2: {}\ncmd: {}\n".format(cmd))
    run_cmd(cmd)
    return paired_1, paired_2

#sra轉換成fastq
def dump_fastq_from_sra(srafile, outdir):
    print ("dump_fastq_from_sra")
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
        p = run_cmd(f'sra-stat -x -s -b 1 -e 2 {self._filepath}')
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

def main():
    #計算時間
    start = time.time()
    #建立argparse.ArgumentParser()物件
    parser = argparse.ArgumentParser("NCBI web assembly pipeline.")
    #告知parser需要的指令引數（選項引數）
    parser.add_argument("--srafile", required=True, help="Path of web file")
    parser.add_argument("--outdir", required=True, help="Output folder")
    parser.add_argument("--tmpdir", default="/tmp", help="Directory of temp folder default: '/tmp'")
    parser.add_argument("--gsize", default='', help="Estimated genome size(MB) eg. 3.2M. default: ''")
    parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    #取得引數傳來的data, args：class `argparse.Namespace`
    args = parser.parse_args()

    #取得引數資料,   vars(args):class dict
    sra_file = args.srafile
    outdir = args.outdir
    tmpdir = args.tmpdir
    threads = args.threads
    gsize = args.gsize
    print(sra_file)

    try:
        sra = SequenceReadArchive(sra_file)
        print("layout:",sra.layout)
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
    print ("layout???\n")
    if sra.layout != '2':
        sys.exit(f'File layout is not pair-end')
    print ("layout=2\n")
    #if sra_layout==2 continue
    #os.path.join(path, *paths)連接路徑
    fastq_dir = os.path.join(outdir,'fastq')
    os.makedirs(fastq_dir, exist_ok = True)
    #解壓縮成fastq
    print('Dump fastq.')
    #run_cmd
    dump_fastq_from_sra(sra_file, fastq_dir)
    #os.listdir(fastq_dir) list files in dir
    forward_reads, reverse_reads = [os.path.join(fastq_dir,i) for i in os.listdir(fastq_dir)]

    # 資料前處理：刪除爛的序列
    # Trimming sequence (trimmomatic)------- Q30 base >= 90% -----------> 預測基因組大小與定序深度(KMC & seqtk)
    print('Trim sequences.')
    r1, r2 = trimming(forward_reads, reverse_reads, fastq_dir, threads)

    #Q30>=90
    if bases_percentage(r1, 30) < 90 and bases_percentage(r2, 30) < 90:
       shutil.rmtree(outdir)
       sys.exit('Reads quality is too low.')

    # 預測基因組大小與定序深度(KMC & seqtk)--- depth>=80 ----> 抽樣(seqtk) -----------> SPAdes
    #                                 --- depth<80 ---------> SPAdes
    #de-novo assembly(SPAdes)------>Polish(pilon)---->Contings(最後成果檔案:conting.fa)
    print("Run assembly pipline 'shovill'")
    assemble_dir = os.path.join(outdir,'assembly_result')
    #depth >= 80
    cmd = f"shovill --R1 {r1} --R2 {r2} --outdir {assemble_dir} --depth 100 --tmpdir . --cpus {threads} --ram 3 --force"
    if gsize:
       cmd += f" --gsize {gsize}"
    print(cmd)
    run_cmd(cmd)
    print('Done,total cost',time.time()-start,'secs')
 
if __name__ == '__main__':
    main()
