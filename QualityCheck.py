from __future__ import print_function
import argparse
import os
import shlex
import subprocess
import sys
import time

# -*- coding: utf-8 -*-


#show program running
def progress_bar(Category):
    for i in range(1, 101):
        print("\r{}: ".format(Category),end="")
        print("[{}] {}%".format("*" * (i // 2), i), end="")
        sys.stdout.flush()
        time.sleep(0.02)
    print ("\n")

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
    start = time.time()
    #get ref_path, qen_path, and outdir
    parser = argparse.ArgumentParser("python3 QualityCheck.py -r {path/to/referenceSequencefiles} -g {path/to/genomefiles} -db {database} -m {mode} -o {outduir}")
    parser.add_argument("-r","--ref", required=True, help="Path of reference Sequence files")
    parser.add_argument("-g","--genome", required=True, help="Path of qenome files")
    parser.add_argument("-db", "--database", required=True, help="busco database")
    parser.add_argument("-m","--mode", required=True, help="busco mode")
    parser.add_argument("-o","--outdir", required=True, help="Output folder")
    args = parser.parse_args()

    refPath=getRefListPath(args.ref,args.outdir)
    genome_Path=getQenomeListPath(args.genome,args.outdir)
    outdir = args.outdir
    outdir_ani=out_txt=os.path.join(outdir, 'fastani')
    out_txt=os.path.join(outdir_ani, 'out.txt')
    db=args.database
    mode=args.mode
    progress_bar("load args")

    #fastANI
    print("-------------------------------fastANI start.-------------------------------")
    print ("reseq: {}\n, qen: {}\n, outdir: {}\n,out_txt: {}".format(refPath, genome_Path, outdir,out_txt))
    progress_bar("fastANI excuting")
    #fasani_=run_cmd("/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI -h")
    fastani_="/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI --rl {} --ql {} -o {}".format(refPath,genome_Path,out_txt)
    print(fastani_+"\n")
    run_cmd(fastani_)
    print("fastANI done.\n")

    # ANI>=95
    print ("-------------------------------fastANI end.-------------------------------\ncompare and calculate ANI\nget ANI refseqPath:\n")
    f = open(out_txt, 'r')
    AverageANI=0.0
    num=0
    not_num=0
    ANI_total=0.0
    ANI_=f.readlines()
    for x in ANI_:
        tmp=float(x.split("\t")[2])
        if(tmp>=95.0):
            num+=1
            ANI_total+=tmp
        else:
            not_num+=1
    AverageANI=ANI_total/num
    print ("Average ANI: {}\ntotal number: {}\n>= quantity: {}\nmax ANI: {}\n".format(AverageANI, num+not_num, num, ANI_[0].split("\t")[2]))
    targetPath=ANI_[0].split("\t")[1]
    #get out.txt line 1 (max ANI)
    print("targetPath: {}\n".format(targetPath))
    f.close()

    #BUSCO
    print("-------------------------------ANI>=95 continue, BUSCO start-------------------------------\n")
    #use conda enterring the busco VM(vm name is "busco")
    run_cmd2("conda activate busco")
    #db="enterobacterales_odb10"
    #mode="geno"
    outdir_bus = out_txt = os.path.join(outdir, 'busco')
    busco_db= os.path.join(outdir, 'busco_db')
    cmd_bus="busco -i {} -o bus_out --out_path {} -l  -m {} --download_path {}".format(targetPath, outdir_bus, db, mode, busco_db)
    print (cmd_bus,"\n")
    run_cmd(cmd_bus)
    run_cmd2("conda deactivate")

    print('Done,total cost', time.time() - start, 'secs\n')

if __name__ == '__main__':
    main()
