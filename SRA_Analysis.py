from __future__ import print_function
import src_.myfunction_ as func_
import src_.myEntrez_ as myEz






def main():
    start = time.time()

    # read command arguments------
    # get ref_path, qen_path, and outdir
    progress_bar("read command")
    parser = argparse.ArgumentParser("python3 ")
    parser.add_argument("--pattern",
                        default="salmonella enterica[ORGN] AND illumina[PLAT] AND wgs[STRA] AND genomic[SRC] AND paired[LAY]",
                        help="Searching condition.")
    # PDAT格式：YYYY/MM/DD
    parser.add_argument("--PDAT", required=True, help="Publication Date[PDAT] of Runs.")
    parser.add_argument("--gsize", default='', help="Estimated genome size(MB) eg. 3.2M. default: ''")
    parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    parser.add_argument("--n", default=3, help="count of download sra file eahc time", type=int)
    parser.add_argument("-r", "--ref", required=True, help="Path of reference Sequence files")
    parser.add_argument("-db", "--busco_database", required=True, help="busco database")
    parser.add_argument("-m", "--mode", default="geno", required=False, help="busco mode")
    parser.add_argument("-o", "--outdir", required=True, help="Output folder")
    args = parser.parse_args()
    referencelist=args.ref
    busco_db=args.busco_datatbase



    print("Download_avaliable_runs_from_NCBI_sra_database_wit_PDAT---------------\n")




    #assemble
    print("Now Assemble excutting------------\n")
    assem=func_.Assemble(PDAT=, outdir=. gsize=, threads=, n=)

    #QualityCheck
    print("Now QualityCheck excutting------------\n")
    #python3 QualityCheck.py -r /data/usrhome/LabSSLin/user30/Desktop/RefSeq/ -g /data/usrhome/LabSSLin/user30/Desktop/SRA/test0812/assembly_result/contigs.fa -db enterobacterales_odb10 -m geno -o /data/usrhome/LabSSLin/user30/Desktop/QualityCheck
    qual_cmd="python3 QualityCheck.py -r {} -g {} -db {} -m geno -o {}".format(referencelist, genomelist, busco_db, ourdir)
    print ("run cmd: {}\n".format(qual_cmd))

    qual= func_.run_cmd(qual_cmd)

    return 0

if __name__ == '__main__':
    main()