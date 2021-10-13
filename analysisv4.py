from __future__ import print_function

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
import pandas as pd
import pymysql

import utils_

#python3 analysisv3.py -i ./contigs.fa -o /data/usrhome/LabSSLin/user30/Desktop/SRA_Analysis/analysis -mlstS senterica -amrS Salmonella
#python3 analysisv3.py -i ./SRAtest/20200704/SRR12144668_contig.fa -o /data/usrhome/LabSSLin/user30/Desktop/SRA_Analysis/analysis -mlstS senterica -amrS Salmonella


def run_cmd(cmd):
    cmd=shlex.split(cmd)
    p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #print (cmd)
    print("--------------------------------------\nSubprogram output:\n")
    utils_.progress_bar("sub excuting")
    while p.poll() is None:
        line = p.stdout.readline()
        line = line.strip()
        if line:
            line_=line.decode().split("\n")
            for s in line_:
                print (str("{}".format(s)))
            sys.stdout.flush()
            sys.stderr.flush()
            print ("\n")
    if p.returncode ==0:
        print ("Subprogram sucess")
    else:
        print ("Subprogram failed")
        print (p.stderr)
    print ("-------------------------\n")
    return p

def main():
    start = time.time()

    ###set Database setting
    db_settings = {
        "host": "127.0.0.1",
        "port": 3306,
        "user": "user",
        "password": "tumvgk01",
        "db": "SRA_Analysis",
        "charset": "utf8"
    }


    # read command arguments------
    # get ref_path, qen_path, and outdir
    utils_.progress_bar("read command")
    #python3 analysis.py -i ./contigs.fa -o /data/usrhome/LabSSLin/user30/Desktop/SRA_Analysis/0906test -mlstS senterica -amrS Salmonella
    parser = argparse.ArgumentParser("python3 analysis.py -i [./contigs.fa] -o [/data/usrhome/LabSSLin/user30/Desktop/SRA_Analysis/0906test] -mlstS [senterica] -amrS [Salmonella]")
    parser.add_argument("--pattern",
                        default="salmonella enterica[ORGN] AND illumina[PLAT] AND wgs[STRA] AND genomic[SRC] AND paired[LAY]",
                        help="Searching condition.")
    # PDAT格式：YYYY/MM/DD
    parser.add_argument("-i", "--input", required=True, help="genome")
    parser.add_argument("-o", "--outdir", required=True, help="Output folder[absolute path]")
    #MLST -s
    parser.add_argument("-mlstS", "--mlstOrganism", help="MLST need -s organism....")
    #Inc type -p
    #parser.add_argument("-p", "--plasmidfinderDB", help="Path of plasmidfinder database")
    #Resistance genes & Mutations identification -d
    parser.add_argument("-amrS", "--amrOrganism", help="Path of amrifinder organism")
    parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    #serotype
    parser.add_argument("-m", "--mode", default="geno", required=False, help="mode")
    args = parser.parse_args()

    input=args.input
    outdir=args.outdir
    mlst_organism=args.mlstOrganism
    #plasmidfinderDB=args.plasmidfinderDB
    amr_organism=args.amrOrganism
    threads=args.threads
    mode=args.mode
    utils_.mkdir_join(outdir)

    #get input id
    inlist=input.split("/")
    inId=inlist[len(inlist)-1]
    print("input Id: {}\n".format(inId))
    inId=inId.split(".")[0]
    print("input Id: {}\n".format(inId))

    #workdir
    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")



    origin_outdir=args.outdir
    allinfopath=os.path.join(origin_outdir,"info.txt")
    check=os.path.join(origin_outdir,"Anacheck.log")

    #add outpath "analysis"
    utils_.mkdir_join(outdir)
    outdir_=os.path.join(outdir,"analysis")
    utils_.mkdir_join(outdir_)
    print("analysis outdir: {}\n".format(outdir_))

    #set {genomoe}_log_output
    logpath = os.path.join(outdir_, inId)
    utils_.mkdir_join(logpath)
    logpath = os.path.join(logpath, "analysis_log.txt")


    #get relative output dir path
    outdir_list=outdir_.split("/")
    relative_path2=outdir_.replace(current_path,".")
    print ("relative2: {}\n".format(relative_path2))
    #relative_path="./"+outdir_list[len(outdir_list)-1]
    print ("relative_path: {}".format(relative_path2))

    relative_path2=os.path.join(relative_path2,inId)
    utils_.mkdir_join(relative_path2)
    print ("relative_path: {}".format(relative_path2))

    #load log.txt read running statedat
    step=0
    filename = Path(logpath)
    filename.touch(exist_ok=True)

    with open(logpath,"r") as f:
        line=f.readlines()
        print (line)
        for x in line:
            ana=x.split(" ")[0]
            if ana == "mlst":
                step=1
            elif ana == "plasmidfinder":
                step=2
            elif ana == "amr":
                step=3
            elif ana == "sistr":
                step=4
            print("ana: {}, step: {}\n".format(ana,step))


    #run MLST
    if step<1:
        print("STEP{}\n".format(step+1))
        print ("********** Now MLST analysis running. **********\n")
        MLST_DB="/data/usrhome/LabSSLin/user30/Desktop/SRA_Analysis/mlst_db"
        #mlst_outdir=os.path.join(outdir,"mlst")
        mlst_outdir = os.path.join(relative_path2, "mlst")
        utils_.mkdir_join(mlst_outdir)
        mlst_datajson=os.path.join(mlst_outdir,"data.json")
        f=open(mlst_datajson,"a+")
        f.close()
        #mlst_cmd="sudo docker run --rm -it \-v {}:/database \-v {}:/workdir \mlst -i {} -o {} -s {}".format(MLST_DB,current_path,input,mlst_outdir,mlst_organism)
        mlst_cmd = "sudo docker run --rm -it \-v {}:/databases \-v {}:/workdir \mlst -i {} -o {} -s {}".format(MLST_DB,
                                                                                                              current_path,
                                                                                                              input,
                                                                                                              mlst_outdir,
                                                                                                              mlst_organism)
        print (mlst_cmd,"\n")
        mlst,err=utils_.run_cmd3(mlst_cmd)
        with open(logpath, "a+") as f:
            if mlst.returncode != 0:
                print(mlst.stdout.readline())
                f.write(err+"\n")
                sys.exit()
            else:
                f.write("mlst is ok\n")
        step += 1
    else:
        print ("**********       mlst was running.      **********\n next step\n")



    #run plasmidfinder
    if step<2:
        print("STEP{}\n".format(step+1))
        print("********** Now plasmidfinder analysis running. **********\n")
        PLASMID_DB="/data/usrhome/LabSSLin/user30/Desktop/SRA_Analysis/plasmidfinder_db"
        #plas_outdir=os.path.join(outdir,"plasmidfinder")
        plas_outdir = os.path.join(relative_path2, "plasmidfinder")
        utils_.mkdir_join(plas_outdir)
        plas_cmd="sudo docker run --rm -it \-v {}:/databases \-v {}:/workdir \plasmidfinder -i {} -o {}".format(PLASMID_DB,current_path,input,plas_outdir)
        print (plas_cmd,"\n")
        plas=run_cmd(plas_cmd)
        with open(logpath, "a+") as f:
            if plas.returncode != 0:
                #print(mlst.stdout.readline())

                sys.exit()
            else:
                f.write("plasmidfinder is ok\n")
        step += 1
    else:
        print("********** plasmidfinder was running. **********\n next step\n")


    #run amrfinder
    if step < 3:

        print("STEP{}\n".format(step+1))
        print("********** Now amrfinder analysis running. **********\n")
        #amr_outdir=os.path.join(outdir,"amrfinder")
        amr_outdir = os.path.join(relative_path2, "amrfinder")
        utils_.mkdir_join(amr_outdir)
        amr_outdir = os.path.join(amr_outdir, "amrout.tsv")
        #amr_cmd="amrfinder -n {} -o {} -O {}".format(input,amr_outdir,amr_organism)
        amr_cmd = "/data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/amrfinder/amrfinder -n {} -o {} -O {}".format(
            input, amr_outdir, amr_organism)
        print(amr_cmd, "\n")
        amr=run_cmd(amr_cmd)
        with open(logpath,"a+") as f:
            if amr.returncode != 0:
                #print (mlst.stdout.readline())
                sys.exit()
            else:
                f.write("amr is ok\n")
        step += 1
    else:
        print("**********   amrfinderr was running.   **********\n next step\n")



    #run sistr
    if step < 4:
        print("STEP{}\n".format(step+1))
        print("********** Now sistr analysis running. **********")
        #sistr_outdir=os.path.join(outdir, "sistr")
        sistr_outdir = os.path.join(relative_path2, "sistr")
        utils_.mkdir_join(sistr_outdir)
        sistr_outdir=os.path.join(sistr_outdir,"sistr_out")
        input_list=input.split("/")
        input_name=input_list[len(input_list)-1]
        print ("name: ",input_name,"\n")
        sistr_cmd="sistr -i {} {} -f csv -o {} -m".format(input, input_name, sistr_outdir)
        print(sistr_cmd,"\n")

        sistr=run_cmd(sistr_cmd)
        with open(logpath, "a+") as f:
            if sistr.returncode != 0:
                #print(mlst.stdout.readline())
                sys.exit()
            else:
                f.write("sistr is ok\n")
        step += 1
    else:
        print("********** sistr was running. **********\n next step\n")


    ########################
    ########################
    #save data in analysis_final.csv and update DB

    #read mlst 'Sequence Type'
    #mlst_file=os.path.join(outdir, "mlst/results.txt")
    mlst_file = os.path.join(relative_path2, "mlst/results.txt")
    with open(mlst_file, "r") as f:
        data = f.readlines()
        print (data[6])
            #Sequence Type: 11
        sequenceType=data[6].split(" ")
        sequenceType=sequenceType[len(sequenceType)-1].strip("\n")

        ####update database table of "SRA" and table of "MLST"
        try:
            conn = pymysql.connect(**db_settings)

            with conn.cursor() as cursor:
                insertSRA = "INSERT INTO SRA(Genome) VALUES(%s);"
                insertMLST = "INSERT INTO MLST(Profile,Organism,SequenceType) VALUES(%s,%s,%s);"
                cursor.execute(
                    insertSRA, (inId))
                cursor.execute(
                    insertMLST,(data[2].split(": ")[1],data[4].split(": ")[1],sequenceType)
                )

        except Exception as e:
            print(e)
        print (sequenceType)

    #read plasmidfinder 'gene'
    plas_file = os.path.join(relative_path2, "plasmidfinder/results_tab.tsv")
    #plas_file=os.path.join(outdir,"plasmidfinder/results_tab.tsv")
    #plas_file = os.path.join(outdir, "./plastest/5524p/results_tab.tsv")
    pladf = pd.read_table(plas_file, sep='\t')
    #print(df)
    pladf=pd.DataFrame(pladf)
    print(pladf)
    print(pladf.columns)
    print (pladf.Plasmid)
    plist=list(pladf.Plasmid)
    plas_format=""

    for x in range(0,len(plist)-1):
        plas_format+=plist[x]
        if x < len(plist)-1:
            plas_format+=","
    print(plas_format)


    #read amrfinder 'Gene symbol', 'subtype'
    ##add amr "Point"
    #amr_file = os.path.join(outdir, "amrfinder/amrout.tsv")
    amr_file = os.path.join(relative_path2, "amrfinder/amrout.tsv")
    # plas_file = os.path.join(outdir, "./plastest/5524p/results_tab.tsv")
    amrdf = pd.read_table(amr_file, sep='\t')
    # print(df)
    amrdf = pd.DataFrame(amrdf)
    print(amrdf)
    print(amrdf.columns)
    #replace ' ' as '_'
    amrdf.columns=['Protein_identifier', 'Contig_id', 'Start', 'Stop', 'Strand',
       'Gene_symbol', 'Sequence_name', 'Scope', 'Element_type',
       'Element_subtype', 'Class', 'Subclass', 'Method', 'Target_length',
       'Reference_sequence_length', 'Coverage_of_reference_sequence',
       'Identity_to_reference_sequence', 'Alignment_length',
       'Accession_of_closest sequence', 'Name_of_closest sequence', 'HMM_id',
       'HMM_description']

    print(amrdf.columns)
    #if list is [], show NaN
    print(amrdf.Gene_symbol)
    print(amrdf.Element_subtype)

    ## arm format: [, , ,][][]
    amr_format=""
    amr_sym=""
    alist=list(amrdf.Gene_symbol)
    for y in range(0,len(alist)):
        amr_format+=alist[y]
        amr_sym+=alist[y]
        if y != len(alist)-1:
            amr_format+=","
            amr_sym+=","

    print(amr_format)
    amr_format += ","
    amr_sub=""
    aalist=list(amrdf.Element_subtype)
    for x in range(0,len(aalist)):
        amr_format+=aalist[x]
        amr_sub +=alist[x]
        if x != len(aalist)-1:
            amr_format += ","
            amr_sub+=","


    #amr_format+=amrdf.Element_subtype
    aaalist=list(amrdf.Method)
    amr_method=""
    for x in range(0,len(aaalist)):
        amr_format+=aaalist[x]
        amr_method+=aaalist[x]
        if x !=len(aaalist)-1:
            amr_format += ","
            amr_method+=","

    #amr_format += amrdf.Method
    try:
        conn = pymysql.connect(**db_settings)

        with conn.cursor() as cursor:
            insertAmr = "INSERT INTO Amrfinder(SRAID,GenSymbol,ElementSubtype,Method) VALUES(%s,%s,%s);"
            cursor.execute(
                insertMLST, (inId,amr_sym,amr_sub,amr_method)
            )

    except Exception as e:
        print(e)


    #read sistr 'serovar'
    #sistr_file = os.path.join(outdir, "sistr/sistr_out.csv")
    sistr_file = os.path.join(relative_path2, "sistr")
    utils_.mkdir_join(sistr_file)
    sistr_file = os.path.join(sistr_file, "sistr_out.csv")
    # plas_file = os.path.join(outdir, "./plastest/5524p/results_tab.tsv")

    sistrdf = pd.read_csv(sistr_file)
    # print(df)
    sistrdf = pd.DataFrame(sistrdf)
    print(sistrdf)
    print(sistrdf.columns)
    print(sistrdf.serovar)
    #dict={'Accession': pd.Series(input for a in range(0,3)),
    #      'mlst':sequenceType,
    #}
    in_abspath=input.replace(".",current_path)

    #dict = {'Accession': input,
    #        'mlst':sequenceType,
    #        'plasmidfinder':pladf.Plasmid,
    #        'amr_gane':amrdf.Gene_symbol,
    #        'amr_subtype':amrdf.Element_subtype,
    #        'sistr':sistrdf.serovar
    #        }

    dict = {'Accession': inId.split("_")[0],
            'MLST': sequenceType,
            'AMR': amr_format,
            #'Point':amrdf.Element_subtype,
            'Serotype': sistrdf.serovar,
            'Inc Type': plas_format
            }

    finaldf=pd.DataFrame(dict)
    print(finaldf)
    finalfile=os.path.join(origin_outdir,"analysis_final.csv")
    #if os.path.isfile(finalfile):
        #beforedf=pd.read_csv(finalfile)
        #beforedf=pd.DataFrame(beforedf)
        #finaldf=pd.concat([beforedf,finaldf])
        #print("merge df\n")
    finaldf.to_csv(finalfile,mode='a+',header=False))
    #finaldf.to_csv(finalfile)



    #after run all state, save ID in "Anackeck.log" and remove ./analysis
    with open(check,"a+") as f:
        f.write("Run {} is ok.\n".format(inId))
    shutil.rmtree(outdir_)
    print("remove ./analysis\n")

    ###info
    #with open(allinfopath,"a+") as f:
        #f.write("analysis outdir = {}\n".format(outdir_))
        #f.write("analysis file path = {}\n".format(finalfile))

if __name__ == '__main__':
    main()