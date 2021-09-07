from __future__ import print_function

import argparse
import os
import time

import src_.myfunction_ as func_







def main():
    start = time.time()

    # read command arguments------
    # get ref_path, qen_path, and outdir
    func_.progress_bar("read command")
    parser = argparse.ArgumentParser("python3 Analysis.py ")
    parser.add_argument("--pattern",
                        default="salmonella enterica[ORGN] AND illumina[PLAT] AND wgs[STRA] AND genomic[SRC] AND paired[LAY]",
                        help="Searching condition.")
    # PDAT格式：YYYY/MM/DD
    parser.add_argument("-i", "--input", required=True, help="genome")
    parser.add_argument("-o", "--outdir", required=True, help="Output folder")
    #MLST -s
    parser.add_argument("-mlstS", "--mlstOrganism", help="MLST need -s organism....")
    #Inc type -p
    parser.add_argument("-p", "--plasmidfinderDB", help="Path of plasmidfinder database")
    #Resistance genes & Mutations identification -d
    parser.add_argument("-amrS", "--amrOrganism", help="Path of amrifinder organism")
    parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    #serotype
    parser.add_argument("-m", "--mode", default="geno", required=False, help="mode")
    args = parser.parse_args()

    input=args.input
    outdir=args.outdir
    mlst_organism=args.mlstOrganism
    plasmidfinderDB=args.plasmidfinderDB
    amr_organism=args.amrOrganism
    threads=args.threads
    mode=args.mode
    func_.mkdir_join(outdir)
    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")
    #run MLST
    print ("Now MLST analysis running--------")
    MLST_DB="/data/usrhome/LabSSLin/user30/Desktop/SRA_Analysis/mlst_db"
    mlst_outdir=os.path.join(outdir,"mlst")
    func_.mkdir_join(mlst_outdir)
    mlst_datajson=os.path.join(mlst_outdir,"data.json")
    f=open(mlst_datajson,"a+")
    f.close()
    mlst_cmd="docker run --rm -it \-v {}:/database \-v {}:/workdir \mlst-2 -i {} -o {} -s {}".format(MLST_DB,current_path,input,mlst_outdir,mlst_organism)
    print (mlst_cmd,"\n")
    func_.run_cmd(mlst_cmd)

    #run plasmidfinder
    print("Now plasmidfinder analysis running--------")
    PLASMID_DB="/data/usrhome/LabSSLin/user30/Desktop/SRA_Analysis/plasmidfinder_db"
    plas_outdir=os.path.join(outdir,"plasmidfinder")
    func_.mkdir_join(plas_outdir)
    plas_cmd="docker run --rm -it \-v {}:/database \-v {}:/workdir \plasmidfinder-2 -i {} -o {}".format(PLASMID_DB,current_path,input,plas_outdir)
    print (plas_cmd,"\n")
    func_.run_cmd(plas_cmd)

    #run amrfinder
    print("Now amrfinder analysis running--------")
    amr_outdir=os.path.join(outdir,"amrfinder")
    func_.mkdir_join(amr_outdir)
    amr_outdir = os.path.join(amr_outdir, "armout.tsv")
    amr_cmd="amrfinder -n {} -o {} -O {}".format(input,amr_outdir,amr_organism)
    print(amr_cmd, "\n")
    func_.run_cmd(amr_cmd)

    #run sistr
    print("Now sistr analysis running--------")
    sistr_outdir=os.path.join(outdir, "sistr")
    func_.mkdir_join(sistr_outdir)
    sistr_outdir=os.path.join(sistr_outdir,"sistr_out")
    input_list=input.split("/")
    input_name=input_list[len(input_list)-1]
    print ("name: ",input_name,"\n")
    sistr_cmd="sistr -i {} {} -f csv -o {} -m".format(input, input_name, sistr_outdir)
    print(sistr_cmd,"\n")
    func_.run_cmd(sistr_cmd)





    return 0

if __name__ == '__main__':
    main()