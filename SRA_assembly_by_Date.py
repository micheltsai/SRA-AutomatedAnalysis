import time
import os
import re
import sys
import shutil
import argparse
import subprocess
import xml.etree.cElementTree as ET
from tempfile import TemporaryDirectory
import pandas as pd
import datetime
import glob
from Bio import Entrez
from io import StringIO

# 調用Enterz時, 需使用email
Entrez.email = 'ann850324@gmail.com'


def count_egquery(term, date_from, date_to, db='sra'):
    pattern = term + f" AND {date_from}[PDAT]:{date_to}[PDAT]"
    print('Searching pattern:', pattern)
    # 查詢PubMed中所有和pattern變數內容相關的文章
    handle = Entrez.egquery(term=pattern)
    d = Entrez.read(handle)
    # 查詢d中的標籤`eGQueryResult`,  將DbName改成Db(`sra`)
    query = list(filter(lambda x: x['DbName'] == db, d['eGQueryResult']))[0]
    print('Total', query['Count'], 'results in NCBI', db, 'database.')
    return pattern, query['Count']


# term<---main.pattern
# db='sra'
# count=query['Count']
def IdList_esearch(term, db, count):
    handle = Entrez.esearch(term=term, db=db, retmax=count)  # 不設定retmax的話只有20筆資料
    d = Entrez.read(handle)
    return d['IdList']


def Get_RunInfo(idlist):
    print('Getting run_info table by idlist...')
    c = len(idlist)
    if c >= 10000:
        print("over 10000 results")
        df_all = pd.DataFrame()
        # 設定k為list,由0開始,最大值為c, 間隔10000
        k = list(range(0, c, 10000))
        k.append(len(idlist))
        print(k)
        for i in range(1, len(k)):
            s = time.time()
            start = k[i - 1]
            end = k[i]
            # 下載GenBank records
            handle = Entrez.efetch(db='sra', id=idlist[start:end], rettype='runinfo', retmode='csv')
            # 查看原始的Genbank文件
            d = handle.read()
            # 讀檔
            df = pd.read_csv(StringIO(d))
            df = df[df['Run'] != 'Run']
            df_all = pd.concat([df_all, df])
            print("finish", i - 1, "round,cost", time.time() - s, "secs")
            print("counts_of_run:", len(df))
            time.sleep(2)
    else:
        handle = Entrez.efetch(db='sra', id=idlist, rettype='runinfo', retmode='csv')
        d = handle.read()
        df = pd.read_csv(StringIO(d))
        df_all = df[df['Run'] != 'Run']
    return df_all


###
def run_cmd(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True)
    return p


def run_cmd2(cmd):
    # p = subprocess.run(cmd, stderr=subprocess.PIPE, shell=True, check=True)
    p = subprocess.run(cmd, shell=True, check=True)
    return p


###for_assembly
def bases_percentage(filepath, qscore=0):
    p = run_cmd(f"seqtk fqchk -q {qscore} {filepath} | grep ALL | awk '{{print $NF}}'")
    return float(p.stdout)


CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
ADAPTERS = os.path.join(CURRENT_DIR, 'trimmomatic.fa')
MIN_BQ = 3


def crop_position(filepath, window_size=3, gap=10):
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
    crop, headcrop = crop_position(forward_reads)
    opt = f"CROP:{crop} HEADCROP:{headcrop} ILLUMINACLIP:{ADAPTERS}:2:30:10 LEADING:{MIN_BQ} TRAILING:{MIN_BQ} " \
          f"SLIDINGWINDOW:4:20 MINLEN:36 TOPHRED33"
    paired_1 = os.path.join(outdir, 'R1.fq')
    paired_2 = os.path.join(outdir, 'R2.fq')
    cmd = f"java -jar trimmomatic-0.39.jar PE -threads {threads} {forward_reads} {reverse_reads} {paired_1} /dev/null" \
          f" {paired_2} /dev/null {opt}"
    run_cmd(cmd)
    return paired_1, paired_2


# Produces two fastq files (--split-files) containing ".1" and ".2" read suffices (-I) for paired-end data.
# 生成兩個fasta文件
def dump_fastq_from_sra(srafile, outdir):
    run_cmd(f'fastq-dump --split-files --outdir {outdir} {srafile}')


class SequenceReadArchive:
    def __init__(self, filepath):
        self._set_filepath(filepath)
        self._get_stat()

    def _set_filepath(self, filepath):
        if os.access(filepath, os.F_OK$) is False:
            raise FileNotFoundError("File not found.")
        with open(filepath, 'rb') as handle:
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

def prefetch_sra(sralist, outdir):
    ss = " ".join(sralist)
    print("now download", ss, "runs.")
    cmd = "prefetch " + ss + " --output-directory " + outdir
    run_cmd2(cmd)


def run_dump_assembly(need_run, sra_dir, assem_dir, threads, gsize, output, check_log, n):
    k = list(range(0, len(need_run), n))
    for i in k:
        run_id = need_run[i:i + n]
        prefetch_sra(run_id, sra_dir)
        time.sleep(1)
        for x in run_id:
            start = time.time()
            path = sra_dir + "/" + x + "/*.sra"
            re_path = "".join(glob.glob(path))
            sra_file = os.path.abspath(re_path)
            try:
                sra = SequenceReadArchive(sra_file)
            except Exception as e:
                sys.exit(e)
            if sra.layout != '2':
                f = open(check_log, "a")
                m = "Run " + x + " layout is not 2."
                f.write(m)
                f.close()
                # sys.exit(f'File layout is not pair-end')
                continue
            print('Now assembly', sra_file, '...')
            outdir = assem_dir + "/" + "".join(x)
            fastq_dir = os.path.join(outdir, 'fastq')
            os.makedirs(fastq_dir, exist_ok=True)
            # print('Dump fastq.')
            dump_fastq_from_sra(sra_file, fastq_dir)
            forward_reads, reverse_reads = [os.path.join(fastq_dir, i) for i in os.listdir(fastq_dir)]
            # print('Trim sequences.')
            r1, r2 = trimming(forward_reads, reverse_reads, fastq_dir, threads)
            if bases_percentage(r1, 30) < 90 and bases_percentage(r2, 30) < 90:
                f = open(check_log, "a")
                m = "Run " + x + " reads quality is too low."
                f.write(m)
                f.close()
                # sys.exit('Reads quality is too low.')
                continue
            print("Run assembly pipline 'shovill'")
            assemble_dir = os.path.join(outdir, 'assembly_result')
            cmd = f"shovill --R1 {r1} --R2 {r2} --outdir {assemble_dir} --depth 100 --tmpdir . --cpus {threads} --ram 3 --force"
            if gsize:
                cmd += f" --gsize {gsize}"
            run_cmd(cmd)
            print('Done.')
            contig_check = assemble_dir + "/contigs.fa"

            if os.path.exists(contig_check) == False:
                f = open(check_log, "a")
                m = "Run " + x + " have assembly error. Didn't have Contigs.fa file.\n"
                f.write(m)
                f.close()
                continue

            # check.log確定資料SRR******已經被抓取下載
            f = open(check_log, "a")
            m = "Run " + x + " is ok.\n"
            f.write(m)
            f.close()

            # 最後希望得到的檔案
            cmd = "mv " + contig_check + " " + assemble_dir + "/" + x + "_contig.fa && mv " + assemble_dir + "/" + x + "_contig.fa " + output
            run_cmd(cmd)
            shutil.rmtree(outdir)
        shutil.rmtree(sra_dir)
        os.mkdir(sra_dir)


def main():
    # 設定命令列參數
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
    # print(date)
    pattern, count = count_egquery(pattern, date, date)
    # print(pattern,count)
    # IdList_esearch(pattern,'sra',count) 回傳d[IdList]
    idlist = IdList_esearch(pattern, 'sra', count)
    print(idlist)

    runinfo = Get_RunInfo(idlist)
    run_list = list(runinfo['Run'])

    sra_dir = args.sra_dir
    assem_dir = args.assembly_dir
    check_log = args.log
    output = args.output
    tmpdir = args.tmpdir
    threads = args.threads
    gsize = args.gsize
    n = args.n

    f = open(check_log, 'a')
    d = f.read().split("\n")[1:-1]
    f.close()
    finish = list(filter(lambda x: len(x.split(" ")) >= 4, d))
    finish_run = list(map(lambda x: x.split(" ")[1], finish))
    need_run = list(filter(lambda x: x not in finish_run, run_list))
    print("Toal", len(need_run), "sra runs need to downlaod.")

    os.makedirs(sra_dir, exist_ok=True)
    os.makedirs(assem_dir, exist_ok=True)
    f = open(check_log, 'a')
    # 輸出當前的日期及時間
    t = str(datetime.datetime.now()).split(".")[0]
    f.write(t)
    f.write("\n")
    f.close()

    if os.path.exists(output) == False:
        os.makedirs(output)

    run_dump_assembly(need_run, sra_dir, assem_dir, threads, gsize, output, check_log, n)


if __name__ == '__main__':
    main()

