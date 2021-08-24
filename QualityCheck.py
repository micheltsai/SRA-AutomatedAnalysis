from __future__ import print_function
import argparse
import os
import shlex
import subprocess
import sys
import time
from datetime import datetime
# -*- coding: utf-8 -*-


#show program running
from builtins import FileExistsError


def progress_bar(Category):
    for i in range(1, 101):
        print("\r{}: ".format(Category),end="")
        print("[{}] {}%".format("*" * (i // 2), i), end="")
        sys.stdout.flush()
        time.sleep(0.02)
    print ("\n")

def mkdir_join(dir, add_dir):
    final=os.path.join(dir, add_dir)
    try:
            os.makedirs(final)
    except FileExistsError:
        print("folder is exist")
    return final

#excute subprogram
def run_cmd2(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True)
    return p


#show subprogram running status
def run_cmd(cmd):
    cmd=shlex.split(cmd)
    p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
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

#get reference list and path of reference list
def getRefListPath(refSeqPath,outdir):
    print("getRefListPath:\n")
    print("refSeqPath: " + refSeqPath + "\n")
    refListPath = os.path.join(outdir, 'ref.txt')
    run_cmd2("find {}>{}".format(refSeqPath,refListPath))
    print("refListPath: "+refListPath+"\n")
    #run_cmd2("cat {}".format(refListPath))
    return refListPath

#get qenome list and path of qenome list
def getQenomeListPath(genome_Path,outdir):
    print("getQenomeListPath:\n")
    print("refSeqPath: " + genome_Path + "\n")
    genListPath = os.path.join(outdir, 'qen.txt')
    run_cmd2("find {}>{}".format(genome_Path,genListPath))
    print("qenListPath: "+genListPath+"\n")
    return genListPath

def main():
    current_time = datetime.now().strftime("%Y%m%d-%I:%M:%S")
    print(current_time)

    start = time.time()

    #read command arguments------
    #get ref_path, qen_path, and outdir
    parser = argparse.ArgumentParser("python3 QualityCheck.py -r {path/to/referenceSequencefiles} -g {path/to/genomefiles} -db {database} -m {mode} -o {outduir}")
    parser.add_argument("-r","--ref", required=True, help="Path of reference Sequence files")
    parser.add_argument("-g","--genome", required=True, help="Path of qenome files")
    parser.add_argument("-db", "--database", required=True, help="busco database")
    parser.add_argument("-m","--mode", required=True, help="busco mode")
    parser.add_argument("-o","--outdir", required=True, help="Output folder")
    args = parser.parse_args()


    outdir = args.outdir
    outdir = mkdir_join(outdir, str(current_time))
    print("outdir: ",outdir)
    refPath = getRefListPath(args.ref, outdir)
    genome_Path = getQenomeListPath(args.genome, outdir)

    outdir_ani=mkdir_join(outdir, 'fastani')
    #outdir_ani=os.path.join(outdir, 'fastani')
    out_txt=os.path.join(outdir_ani, 'out.txt') #stroed fastANI output in out.txt

    db=args.database
    mode=args.mode
    progress_bar("load args")

    #fastANI-------
    print("-------------------------------fastANI start.-------------------------------")
    print ("reseq: {}\n, qen: {}\n, outdir: {}\n,out_txt: {}".format(refPath, genome_Path, outdir,out_txt))
    progress_bar("fastANI excuting")
    #fasani_=run_cmd("/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI -h")
    fastani_="/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI --rl {} --ql {} -o {}".format(refPath,genome_Path,out_txt)
    print(fastani_+"\n")
    run_cmd(fastani_)
    print("fastANI done.\n")

    # ANI>=95------
    print ("-------------------------------fastANI end.-------------------------------\ncompare and calculate ANI\nget ANI refseqPath:\n")
    #open fastANI output
    f = open(out_txt, 'r')
    AverageANI=0.0
    num=0   #quantity of ANI>=95
    not_num=0 #quantity of ANI<95
    ANI_total=0.0 #total of all ANI value
    ANI_=f.readlines() #read file stored in ANI_
    for x in ANI_:
        tmp=float(x.split("\t")[2])#temporarily ANI value
        if(tmp>=95.0):
            num+=1
            ANI_total+=tmp
        else:
            not_num+=1
    AverageANI=ANI_total/num #if ANI>=95 ,calulate average value of ANI
    print ("Average ANI: {}\ntotal number: {}\n>= quantity: {}\nmax ANI: {}\n".format(AverageANI, num+not_num, num, ANI_[0].split("\t")[2]))
    targetPath=ANI_[0].split("\t")[1]
    #get out.txt line 1 (max ANI)
    print("targetPath: {}\n".format(targetPath))
    f.close()

    #BUSCO------
    print("-------------------------------ANI>=95 continue, BUSCO start-------------------------------\n")
    #use conda enterring the busco VM(vm name is "busco")

    #db="enterobacterales_odb10"
    #mode="geno"
    outdir_bus = mkdir_join(outdir, 'busco')

    #outdir_bus = os.path.join(outdir, 'busco')
    busco_db = mkdir_join(args.outdir, 'busco_db')
    #busco_db= os.path.join(outdir, 'busco_db')
    # subprocess.run('bash -c "conda activate busco"', shell=True)
    #run_cmd('bash -c "source /data/usrhome/LabSSLin/user30/anaconda3/etc/profile.d/conda.sh && conda activate busco"')
    #-f overwrite
    cmd_bus = 'bash -c "source /data/usrhome/LabSSLin/user30/anaconda3/etc/profile.d/conda.sh && conda activate busco && busco -i {} -o bus_out --out_path {} -l {} -m {} --download_path {} -f"'.format(targetPath, outdir_bus,
                                                                                              db, mode, busco_db)
    #cmd_bus="busco -i {} -o bus_out --out_path {} -l {} -m {} --download_path {} -f".format(targetPath, outdir_bus, db, mode, busco_db)
    print (cmd_bus,"\n")
    run_cmd(cmd_bus)
    #subprocess.run('bash -c "conda deactivate"',shell=True)

    print('Done,total cost', time.time() - start, 'secs\n')

if __name__ == '__main__':
    main()
