from __future__ import print_function
import argparse
import glob
import os
import re
import shlex
import subprocess
import sys
import time
from datetime import datetime
# -*- coding: utf-8 -*-
#python3 QualityCheck.py -r /data/usrhome/LabSSLin/user30/Desktop/RefSeq/ -g /data/usrhome/LabSSLin/user30/Desktop/SRA/test0812/assembly_result/contigs.fa -db enterobacterales_odb10 -m geno -o /data/usrhome/LabSSLin/user30/Desktop/QualityCheck

#show program running
from builtins import FileExistsError

import numpy as np

import utils_


###### outdir ./PADT/Qualitycheck/
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
    utils_.mkdir_join(outdir)
    refPath = utils_.getRefListPath(args.ref, outdir)
    Assem_path = os.path.join(outdir, "Assembled/")
    BUSCOresult = os.path.join(outdir,"BUSCOresult.txt")
    check = os.path.join(outdir, "QCcheck.log")
    outdir = os.path.join(outdir,"QualityCheck")
    utils_.mkdir_join(outdir)

    #outdir = utils_.mkdir_join(outdir, str(current_time))
    print("outdir: \n",outdir)
    print("check: \n",check)
    print("BUSCOresult= {}".format(BUSCOresult))



    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")
    #genome_Path = utils_.getGenomeListPath(args.genome, outdir)

    genome_Path=args.genome

    gID=genome_Path.replace(Assem_path,"")

    gID=gID.split(".")[0]
    print("gID: {}\n".format(gID))

    outdir_ani = os.path.join(outdir, 'fastani')
    utils_.mkdir_join(outdir_ani)
    #print("outdir_ani: {}\n".format(outdir_ani))
    # outdir_ani=os.path.join(outdir, 'fastani')

    outfile='{}_ani.txt'.format(gID)

    outfile = os.path.join(outdir_ani, outfile)  # stroed fastANI output in out.txt
    info_txt = os.path.join(outdir_ani, '{}_info.txt'.format(gID))  # stroed fastANI output in out.txt
    db=args.database
    mode=args.mode
    utils_.progress_bar("load args")

    #fastANI-------

    print("-------------------------------fastANI start.-------------------------------")
    print ("reseq: {}\n qen: {}\n outdir: {}\nout_txt: {}".format(refPath, genome_Path, outdir, outfile))
    utils_.progress_bar("fastANI excuting")
    #fasani_=run_cmd("/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI -h")
    #fastani_="/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI --rl {} --ql {} -o {}".format(refPath,genome_Path,out_txt)
    #fastani_ = "/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI --rl {} -q {} -o {}".format(refPath,
    fastani_ = "/data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/FastANI/fastANI --rl {} -q {} -o {}".format(refPath,
                                                                                                    genome_Path,
                                                                                                    os.path.join("./SRAtest/20200704/QualityCheck/fastani", '{}_ani.txt'.format(gID)))
    print(fastani_+"\n")
    utils_.run_cmd(fastani_)
    print("fastANI done.\n")

    # ANI>=95------
    print ("-------------------------------fastANI end.-------------------------------\ncompare and calculate ANI\nget ANIoutPath\n")
    #open fastANI output
    f = open(outfile, 'r+')
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
    # num =0
    AverageANI=ANI_total/num #if ANI>=95 ,calulate average value of ANI
    print ("Average ANI: {}\ntotal number: {}\n>= quantity: {}\nmax ANI: {}\n".format(AverageANI, num+not_num, num, ANI_[0].split("\t")[2]))
    targetPath=ANI_[0].split("\t")[1]

    with open(outfile,"a") as f:
        f.write("\ncalculate ANI:\n")
        f.write("Average ANI: {}\ntotal number: {}\n>= quantity: {}\nmax ANI: {}\n".format(AverageANI, num+not_num, num, ANI_[0].split("\t")[2]))
    #get out.txt line 1 (max ANI)
    print("targetPath: {}\n".format(targetPath))

    #save data info
    with open(info_txt,"w+") as f2:
        f2.write("Average ANI: {}\ntotal number: {}\n>= quantity: {}\nmax ANI: {}\n".format(AverageANI, num+not_num, num, ANI_[0].split("\t")[2]))
        f2.write("targetPath: {}\n".format(targetPath))

    f.close()

    if num == 0:
        print("All ANI value doesn't exceed 95, next genome run\n")
        with open(check, "a+") as f:
            f.write("{} is ANI<95.\n".format(gID))
        return 0


    #BUSCO------
    print("-------------------------------ANI>=95 continue, BUSCO start-------------------------------\n")
    #use conda enterring the busco VM(vm name is "busco")

    #db="enterobacterales_odb10"
    #mode="geno"

    outdir_bus = os.path.join(outdir, 'busco_db')
    busco_db = utils_.mkdir_join(outdir_bus)
    #busco_db= os.path.join(outdir, 'busco_db')
    # subprocess.run('bash -c "conda activate busco"', shell=True)
    #run_cmd('bash -c "source /data/usrhome/LabSSLin/user30/anaconda3/etc/profile.d/conda.sh && conda activate busco"')
    #-f overwrite


    #genome is "one excuting"
    #busco -i /data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/RefSeq/GCF_000335875.2.fa -o cofig --out_path /data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/SRA-AutomatedAnalysis/QualityCheck -l enterobacterales_odb10 -m geno --download_path /data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/SRA-AutomatedAnalysis/QualityCheck/QualityCheck/busco_db -f
    #cmd_bus = 'bash -c "source /data/usrhome/LabSSLin/user30/anaconda3/etc/profile.d/conda.sh && conda activate busco && busco -i {} -o {} --out_path {} -l {} -m {} --download_path {} -f"'.format(targetPath, gID, outdir, db, mode, busco_db)
    cmd_bus = 'busco -i {} -o {} --out_path {} -l {} -m {} --download_path {} -f"'.format(
        targetPath, gID, outdir, db, mode, busco_db)
    #cmd_bus="busco -i {} -o bus_out --out_path {} -l {} -m {} --download_path {} -f".format(targetPath, outdir_bus, db, mode, busco_db)
    print (cmd_bus,"\n")
    utils_.run_cmd(cmd_bus)
    #subprocess.run('bash -c "conda deactivate"',shell=True)

    #get BUSCO complete>=95 & duplicate>=3 ,or exit
    buscopath=os.path.join(outdir,"{}".format(gID))
    buscopath=os.path.join(buscopath,"run_{}".format(db))
    buscopath=glob.glob(buscopath+"/*.txt")
    print(buscopath)
    buscopath=os.path.abspath(buscopath[0])
    print(buscopath)
    with open(buscopath,"r")as f:
        b=f.readlines()
    print(b,"\n")
    b=b[8].strip("\t")
    print(b)
    print([float(s) for s in re.findall(r'-?\d+\.?\d*',b)])

    bC, bS, bD, bF, bM, bn=[float(s) for s in re.findall(r'-?\d+\.?\d*',b)]
    #bC, bS, bD, bF, bM =np.fromstring(b,dtype=float,sep=":")
    #bN=np.fromstring(b,dtype=int,sep="n:")

    print("c:{}%, d:{}%\n".format(bC,bD))

    with open(BUSCOresult,"a+")as f:
        print("stored BUSCO result\n")
        f.write("{}: {}\n".format(gID, [float(s) for s in re.findall(r'-?\d+\.?\d*',b)]))

    if bC<95.0 and bD > 3.0:
        with open(check, "a+") as f:
            f.write("{} is C<95 or D>3.\n".format(gID))
        return 0

    #continue
    targettxt=os.path.join(args.outdir, "target.txt")
    print("target.txt path: {}".format(targettxt))
    #continue need target path
    with open (targettxt, "a+") as f:
        f.write("{}:{}\n".format(genome_Path,targetPath))
    #check
    with open(check,"a+") as f:
        f.write("{} is ok.\n".format(gID))
        print("commit on check \n")
    print('Done,total cost', time.time() - start, 'secs\n')
    return targetPath
if __name__ == '__main__':
    main()
