from __future__ import print_function

import argparse
import os
import time
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
    #parser.add_argument("-i", "--input", required=True, help="genome")
    #parser.add_argument("-o", "--outdir", required=True, help="Output folder")
    #MLST -s
    #parser.add_argument("--s", "--organism", help="MLST need -s organism")
    #Inc type -p
    #parser.add_argument("-p", "--plasmidfinderDB", required=True, help="Path of plasmidfinder database")
    #Resistance genes & Mutations identification -d
    #parser.add_argument("-d", "--amrfinderDB", required=True, help="Path of amrifinder database")
    #parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    #serotype
    #parser.add_argument("-m", "--mode", default="geno", required=False, help="busco mode")
    args = parser.parse_args()
    #referencelist=args.ref
    #busco_db=args.busco_datatbase
    date=args.PDAT

    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")

    ## read SRAsetting.txt
    utils_.progress_bar("read SRAsetting.txt")
    setting_path=os.path.join(current_path,"SRAsettings.txt")
    with open(setting_path,"r") as f:
        setList=f.readlines()

    print(setList)
    i=0
    for line in setList:
        line=line.strip("\n")
        line_=line.split("=")
        if line != "" and len(line_)==2:
            print(line_)
            print("line{}. {}:{}\n".format(i, line_[0], line_[1]))
        i+=1

    thread=setList[1].strip("\n").split("=")[1]
    gsize=setList[4].strip("\n").split("=")[1]
    n=setList[7].strip("\n").split("=")[1]
    outdir=setList[10].strip("\n").split("=")[1]
    ref_dir=setList[12].strip("\n").split("=")[1]
    buscoDB=setList[13].strip("\n").split("=")[1]
    busco_mode = setList[14].strip("\n").split("=")[1]
    mlstS = setList[16].strip("\n").split("=")[1]
    amrS = setList[17].strip("\n").split("=")[1]
    print("thread:{}\ngsize:{}\nn:{}\noutdir:{}\nref_dir:{}\nbuscoDB:{}\nbusco_mode:{}\nmlstS:{}\n,amrS:{}\n".format(
        thread,gsize,n,outdir,ref_dir,buscoDB,busco_mode,mlstS,amrS))

    utils_.mkdir_join(outdir)
    #./YYYYMMDD
    pdat=date.replace("/","")
    outdir=os.path.join(outdir,pdat)
    utils_.mkdir_join(outdir)

    #./YYYYMMDD/Assembled
    #assembly_date.py
    assem_cmd="python3 assembly_datev3.py --PDAT {} --output {}".format(date,outdir)
    print(assem_cmd)
    utils_.run_cmd(assem_cmd)

    assembled_path=os.path.join(outdir,"Assembled")
    assembled_file=utils_.getGenomeListPath(assembled_path,assembled_path)
    #QualityCheck
    print("Now QualityCheck excutting------------\n")
    #python3 QualityCheck.py -r /data/usrhome/LabSSLin/user30/Desktop/RefSeq/ -g /data/usrhome/LabSSLin/user30/Desktop/SRA/test0812/assembly_result/contigs.fa -db enterobacterales_odb10 -m geno -o /data/usrhome/LabSSLin/user30/Desktop/QualityCheck
    with open(assembled_file,"r") as f:
        line=f.readlines()
    print(line)



    qual_cmd="python3 QualityCheck.py -r {} -g {} -db {} -m geno -o {}".format(ref_dir, g, busco_db, ourdir)
    #print ("run cmd: {}\n".format(qual_cmd))

    #qual= utils_.run_cmd(qual_cmd)

    return 0

if __name__ == '__main__':
    main()