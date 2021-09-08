from __future__ import print_function
import argparse
import os
import shlex
import subprocess
import sys
import time
from datetime import datetime
# -*- coding: utf-8 -*-
#python3 QualityCheck.py -r /data/usrhome/LabSSLin/user30/Desktop/RefSeq/ -g /data/usrhome/LabSSLin/user30/Desktop/SRA/test0812/assembly_result/contigs.fa -db enterobacterales_odb10 -m geno -o /data/usrhome/LabSSLin/user30/Desktop/QualityCheck

#show program running
from builtins import FileExistsError
import utils_

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
    outdir = utils_.mkdir_join(outdir, str(current_time))
    print("outdir: ",outdir)
    refPath = utils_.getRefListPath(args.ref, outdir)
    genome_Path = utils_.getQenomeListPath(args.genome, outdir)

    outdir_ani=utils_.mkdir_join(outdir, 'fastani')
    #outdir_ani=os.path.join(outdir, 'fastani')
    out_txt=os.path.join(outdir_ani, 'out.txt') #stroed fastANI output in out.txt
    info_txt = os.path.join(outdir_ani, 'ANIinfo.txt')  # stroed fastANI output in out.txt
    db=args.database
    mode=args.mode
    utils_.progress_bar("load args")

    #fastANI-------
    print("-------------------------------fastANI start.-------------------------------")
    print ("reseq: {}\n, qen: {}\n, outdir: {}\n,out_txt: {}".format(refPath, genome_Path, outdir,out_txt))
    utils_.progress_bar("fastANI excuting")
    #fasani_=run_cmd("/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI -h")
    fastani_="/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI --rl {} --ql {} -o {}".format(refPath,genome_Path,out_txt)
    print(fastani_+"\n")
    utils_.run_cmd(fastani_)
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

    #save data info
    with open(info_txt,"w+") as f2:
        f2.write("Average ANI: {}\ntotal number: {}\n>= quantity: {}\nmax ANI: {}\n".format(AverageANI, num+not_num, num, ANI_[0].split("\t")[2]))
        f2.write("targetPath: {}\n".format(targetPath))

    f.close()

    #BUSCO------
    print("-------------------------------ANI>=95 continue, BUSCO start-------------------------------\n")
    #use conda enterring the busco VM(vm name is "busco")

    #db="enterobacterales_odb10"
    #mode="geno"

    #outdir_bus = os.path.join(outdir, 'busco')
    busco_db = utils_.mkdir_join(args.outdir, 'busco_db')
    #busco_db= os.path.join(outdir, 'busco_db')
    # subprocess.run('bash -c "conda activate busco"', shell=True)
    #run_cmd('bash -c "source /data/usrhome/LabSSLin/user30/anaconda3/etc/profile.d/conda.sh && conda activate busco"')
    #-f overwrite
    cmd_bus = 'bash -c "source /data/usrhome/LabSSLin/user30/anaconda3/etc/profile.d/conda.sh && conda activate busco && busco -i {} -o busco --out_path {} -l {} -m {} --download_path {} -f"'.format(targetPath, outdir,
                                                                                              db, mode, busco_db)
    #cmd_bus="busco -i {} -o bus_out --out_path {} -l {} -m {} --download_path {} -f".format(targetPath, outdir_bus, db, mode, busco_db)
    print (cmd_bus,"\n")
    utils_.run_cmd(cmd_bus)
    #subprocess.run('bash -c "conda deactivate"',shell=True)

    print('Done,total cost', time.time() - start, 'secs\n')

if __name__ == '__main__':
    main()
