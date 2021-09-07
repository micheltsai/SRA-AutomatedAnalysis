from __future__ import print_function

import argparse
from datetime import time

import src_.myfunction_ as func_
import src_.myEntrez_ as myEz






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
    parser.add_argument("--s", "--organism", help="MLST need -s organism")
    #Inc type -p
    parser.add_argument("-p", "--plasmidfinderDB", required=True, help="Path of plasmidfinder database")
    #Resistance genes & Mutations identification -d
    parser.add_argument("-d", "--amrfinderDB", required=True, help="Path of amrifinder database")
    parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    #serotype
    parser.add_argument("-m", "--mode", default="geno", required=False, help="busco mode")
    args = parser.parse_args()
    referencelist=args.ref
    busco_db=args.busco_datatbase



    print("Download_avaliable_runs_from_NCBI_sra_database_wit_PDAT---------------\n")




    #assemble
    print("Now Assemble excutting------------\n")
    assem=func_.Assemble(PDAT, outdir, gsize, threads=8, n=3)

    #QualityCheck
    print("Now QualityCheck excutting------------\n")
    #python3 QualityCheck.py -r /data/usrhome/LabSSLin/user30/Desktop/RefSeq/ -g /data/usrhome/LabSSLin/user30/Desktop/SRA/test0812/assembly_result/contigs.fa -db enterobacterales_odb10 -m geno -o /data/usrhome/LabSSLin/user30/Desktop/QualityCheck
    qual_cmd="python3 QualityCheck.py -r {} -g {} -db {} -m geno -o {}".format(referencelist, genomelist, busco_db, ourdir)
    print ("run cmd: {}\n".format(qual_cmd))

    qual= func_.run_cmd(qual_cmd)

    return 0

if __name__ == '__main__':
    main()