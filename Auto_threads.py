from __future__ import print_function

import csv
import datetime
import glob
import multiprocessing
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path
import pandas as pd

import analysisv4
import utils_


def run_cmd(cmd):
    cmd=shlex.split(cmd)
    print(cmd)
    p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #print (cmd)
    print("--------------------------------------\nSubprogram output:\n")
    while p.poll() is None:
        #progress_bar("sub excuting")
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

def progress_bar(Category):
    for i in range(1, 101):
        print("\r{}: ".format(Category),end="")
        print("[{}] {}%".format("*" * (i // 2), i), end="")
        sys.stdout.flush()
        time.sleep(0.02)
    print ("\n")

def Download(x):
    one_ = time.time()
    #print(
    #   "---------------------\n---------------------[ {} / {} ]---------------------\n".format(num + 1,
    #                                                                                           len(idlist)))
    #num += 1
    print("x = {}".format(x))
    # outdir__ = os.path.join(output, "out")
    outdir__ = os.path.join(new_outdir, "Assembled")

    final_dir = os.path.join(outdir__, "{}_contig.fa".format(x))
    sra_file=os.path.join(sra_dir,"{}.sra".format(x))
    if os.path.isfile(final_dir):
        print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
    elif os.path.isfile(sra_file):
        print("was ran download ,sra is exist\n------------------------------\n\n")
    else:
        utils_.prefetch_sra(x, sra_dir)
        print("Download {}\n.".format(x))
        #with open(Downloadcheck_log, "a+") as f:
        #    f.write("{}\n".format(x))
    dltime=time.time() - one_
    print('Done,total cost',dltime, 'secs')
    print("###########################################################")



def Assembled(x):
    final_dir = os.path.join(ass_dir, "{}_contig.fa".format(x))
    if os.path.isfile(final_dir):
        print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
    else:
        # one_run_prefetch = time.time()
        # utils_.prefetch_sra(x,sra_dir)
        one_run_ass = time.time()
        sra_file = os.path.join(sra_dir, "{}/{}.sra".format(x, x))
        print(sra_file)
        if os.path.isfile(sra_file):
            print("have {}.sra file.\n".format(x))
        else:
            utils_.prefetch_sra(x, sra_dir)
            print("not found {}.sra, Download now\n".format(x))

        utils_.run_for_114(x, sra_dir, fastq_dir, assemble_dir, new_outdir, thread, gsize, start, check_log)
        ### unnecessary ERR file
        #ERR_path = os.path.join(os.path.abspath(os.getcwd()), x)
        #print("suERR_path: ", ERR_path, "\n")
        # print ("shutil.rmtree({})\n".format(current_path))
        #utils_.run_cmd2("rm -rf {}".format(current_path))
        #print("remove {}\n".format(current_path))
        rmsra_cmd="rm -rf {}".format(sra_file)
        print(rmsra_cmd)
        run_cmd(rmsra_cmd)


def QualityCheck(genome_Path):
    print("#####################  QualityCheck  #####################\n")
    refPath = utils_.getRefListPath(ref_dir, new_outdir)
    # refPath=args.ref
    Assem_path = os.path.join(new_outdir, "Assembled/")
    BUSCOresult = os.path.join(new_outdir, "BUSCOresult.txt")
    check = os.path.join(new_outdir, "QCcheck.log")
    outdir = os.path.join(new_outdir, "QualityCheck")
    utils_.mkdir_join(outdir)

    # outdir = utils_.mkdir_join(outdir, str(current_time))
    print("outdir: \n", new_outdir)
    print("check: \n", check)
    print("BUSCOresult= {}".format(BUSCOresult))

    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")
    # genome_Path = utils_.getGenomeListPath(args.genome, outdir)



    gID = genome_Path.replace(Assem_path, "")
    # gID=os.path.basename(genome_Path)
    gID = gID.split(".")[0]
    print("gID: {}\n".format(gID))

    # outdir_ani = os.path.join(outdir, 'fastani')
    outdir_ani = os.path.join(outdir, 'fastani')
    utils_.mkdir_join(outdir_ani)
    print("outdir_ani: {}\n".format(outdir_ani))
    # outdir_ani=os.path.join(outdir, 'fastani')

    outfile_ = '{}_ani.txt'.format(gID)

    outfile = os.path.join(outdir_ani, outfile_)  # stroed fastANI output in out.txt

    info_txt = os.path.join(outdir_ani, '{}_info.txt'.format(gID))  # stroed fastANI output in out.txt
    db = buscoDB
    mode = buscoMode
    utils_.progress_bar("load args")
    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "read args", "time": str(time.time() - start)})

    ##fastANI-------
    fastANI_time = time.time()
    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")
    replace_path = outdir.replace(current_path, ".")
    # fastani_outdir=os.path.join("./SRAtest/20200704/QualityCheck/fastani", '{}_ani.txt'.format(gID))
    fastani_outdir = os.path.join(replace_path, '{}_ani.txt'.format(gID))

    print("-------------------------------fastANI start.-------------------------------")
    print("reseq: {}\n qen: {}\n outdir: {}\nout_txt: {}\n{}\n".format(refPath, genome_Path, outdir, outfile,
                                                                       os.path.join(outdir_ani, outfile_)))
    utils_.progress_bar("fastANI excuting")
    # fasani_=run_cmd("/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI -h")
    # fastani_="/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI --rl {} --ql {} -o {}".format(refPath,genome_Path,out_txt)
    # fastani_ = "sudo /data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/FastANI/fastANI --rl {} -q {} -o {}".format(refPath,
    fastani_ = "/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI --rl {} -q {} -o {}".format(refPath, genome_Path, outfile)
    print(fastani_ + "\n")
    if os.path.isfile(outfile):
        print(outfile, " is exist.\n")
        print("fastANI was done.\n")
    else:
        utils_.run_cmd(fastani_)
        print("fastANI done.\n")

    # ANI>=95------
    print(
        "-------------------------------fastANI end.-------------------------------\ncompare and calculate ANI\nget ANIoutPath\n")
    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "{} fastANI".format(gID), "time": str(time.time() - start)})

    # open fastANI output
    f = open(outfile, 'r')
    AverageANI = 0.0
    num = 0  # quantity of ANI>=95
    not_num = 0  # quantity of ANI<95
    ANI_total = 0.0  # total of all ANI value
    ANI_ = f.readlines()  # read file stored in ANI_
    for x in ANI_:
        tmp = float(x.split("\t")[2])  # temporarily ANI value
        if (tmp >= 95.0):
            num += 1
            ANI_total += tmp
        else:
            not_num += 1
    # num =0
    AverageANI = ANI_total / num  # if ANI>=95 ,calulate average value of ANI
    print("Average ANI: {}\ntotal number: {}\n>= quantity: {}\nmax ANI: {}\n".format(AverageANI, num + not_num, num,
                                                                                     ANI_[0].split("\t")[2]))
    targetPath = ANI_[0].split("\t")[1]

    # save data info
    with open(info_txt, "w+") as f2:
        f2.write(
            "Average ANI: {}\ntotal number: {}\n>= quantity: {}\nmax ANI: {}\n".format(AverageANI, num + not_num, num,
                                                                                       ANI_[0].split("\t")[2]))
        f2.write("targetPath: {}\n".format(targetPath))

    f.close()

    if num == 0:
        print("All ANI value doesn't exceed 95, next genome run\n")
        with open(check, "a+") as f:
            f.write("{} is ANI<95.\n".format(gID))
        return 0

    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "{} fastANI".format(gID), "time": str(time.time() - fastANI_time)})

    # BUSCO------
    print("-------------------------------ANI>=95 continue, BUSCO start-------------------------------\n")
    # use conda enterring the busco VM(vm name is "busco")

    # db="enterobacterales_odb10"
    # mode="geno"
    busco_time = time.time()
    outdir_bus = os.path.join(outdir, 'busco_db')
    busco_db = utils_.mkdir_join(outdir_bus)
    # busco_db= os.path.join(outdir, 'busco_db')
    # subprocess.run('bash -c "conda activate busco"', shell=True)
    # run_cmd('bash -c "source /data/usrhome/LabSSLin/user30/anaconda3/etc/profile.d/conda.sh && conda activate busco"')
    # -f overwrite

    # genome is "one excuting"
    # busco -i /data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/RefSeq/GCF_000335875.2.fa -o cofig --out_path /data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/SRA-AutomatedAnalysis/QualityCheck -l enterobacterales_odb10 -m geno --download_path /data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/SRA-AutomatedAnalysis/QualityCheck/QualityCheck/busco_db -f
    cmd_bus = 'bash -c "source /data/usrhome/LabSSLin/user30/anaconda3/etc/profile.d/conda.sh && conda activate busco && busco -i {} -o {} --out_path {} -l {} -m {} --download_path {} -f"'.format(
        targetPath, gID, outdir, db, mode, busco_db)
    # cmd_bus = 'busco -i {} -o {} --out_path {} -l {} -m {} --download_path {} -f"'.format(
    #   targetPath, gID, outdir, db, mode, busco_db)
    # cmd_bus="busco -i {} -o bus_out --out_path {} -l {} -m {} --download_path {} -f".format(targetPath, outdir_bus, db, mode, busco_db)
    print(cmd_bus, "\n")
    utils_.run_cmd(cmd_bus)
    # subprocess.run('bash -c "conda deactivate"',shell=True)

    # get BUSCO complete>=95 & duplicate>=3 ,or exit
    buscopath = os.path.join(outdir, "{}".format(gID))
    buscopath = os.path.join(buscopath, "run_{}".format(db))
    buscopath = glob.glob(buscopath + "/*.txt")
    print(buscopath)
    buscopath = os.path.abspath(buscopath[0])
    print(buscopath)
    with open(buscopath, "r") as f:
        b = f.readlines()
    print(b, "\n")
    b = b[8].strip("\t")
    print(b)
    print([float(s) for s in re.findall(r'-?\d+\.?\d*', b)])

    bC, bS, bD, bF, bM, bn = [float(s) for s in re.findall(r'-?\d+\.?\d*', b)]
    # bC, bS, bD, bF, bM =np.fromstring(b,dtype=float,sep=":")
    # bN=np.fromstring(b,dtype=int,sep="n:")

    print("c:{}%, d:{}%\n".format(bC, bD))

    with open(BUSCOresult, "a+") as f:
        print("stored BUSCO result\n")
        f.write("{}: {}\n".format(gID, [float(s) for s in re.findall(r'-?\d+\.?\d*', b)]))

    if bC < 95.0 and bD > 3.0:
        with open(check, "a+") as f:
            f.write("{} is C<95 or D>3.\n".format(gID))
        return 0

    # continue
    targettxt = os.path.join(new_outdir, "target.txt")
    print("target.txt path: {}".format(targettxt))
    # continue need target path
    with open(targettxt, "a+") as f:
        f.write("{}:{}\n".format(genome_Path, targetPath))
    # check
    with open(check, "a+") as f:
        f.write("{} is ok.\n".format(gID))
        print("commit on check \n")

    with open("./ana_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "{} busco".format(gID), "time": str(time.time() - busco_time)})
    print('Done,total cost', time.time() - start, 'secs\n')
    return targetPath

def Analysis(input,anoutdir):
    print("#####################  Analysis  #####################\n")
    anoutdir = new_outdir
    mlst_organism = mlstS
    # plasmidfinderDB=args.plasmidfinderDB
    amr_organism = amrS
    utils_.mkdir_join(anoutdir)

    # get input id
    inlist = input.split("/")
    inId = inlist[len(inlist) - 1]
    print("input Id: {}\n".format(inId))
    inId = inId.split(".")[0]
    print("input Id: {}\n".format(inId))

    # workdir
    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")

    origin_outdir = new_outdir
    allinfopath = os.path.join(origin_outdir, "info.txt")
    check = os.path.join(origin_outdir, "Anacheck.log")

    # add outpath "analysis"
    utils_.mkdir_join(outdir)
    anoutdir_ = os.path.join(outdir, "analysis")
    utils_.mkdir_join(anoutdir_)
    print("analysis outdir: {}\n".format(anoutdir_))

    # set {genomoe}_log_output
    logpath = os.path.join(anoutdir_, inId)
    utils_.mkdir_join(logpath)
    logpath = os.path.join(logpath, "analysis_log.txt")

    # get relative output dir path
    outdir_list = anoutdir_.split("/")
    relative_path2 = anoutdir_.replace(current_path, ".")
    print("relative2: {}\n".format(relative_path2))
    # relative_path="./"+outdir_list[len(outdir_list)-1]
    print("relative_path: {}".format(relative_path2))

    relative_path2 = os.path.join(relative_path2, inId)
    utils_.mkdir_join(relative_path2)
    print("relative_path: {}".format(relative_path2))

    # load log.txt read running statedat
    step = 0
    filename = Path(logpath)
    filename.touch(exist_ok=True)

    with open(logpath, "r") as f:
        line = f.readlines()
        print(line)
        for x in line:
            ana = x.split(" ")[0]
            if ana == "mlst":
                step = 1
            elif ana == "plasmidfinder":
                step = 2
            elif ana == "amr":
                step = 3
            elif ana == "sistr":
                step = 4
            print("ana: {}, step: {}\n".format(ana, step))

    # run MLST
    if step < 1:
        step1_time = time.time()
        print("STEP{}\n".format(step + 1))
        print("********** Now MLST analysis running. **********\n")
        MLST_DB = "/data/usrhome/LabSSLin/user30/Desktop/SRA_Analysis/mlst_db"
        # mlst_outdir=os.path.join(outdir,"mlst")
        mlst_outdir = os.path.join(relative_path2, "mlst")
        utils_.mkdir_join(mlst_outdir)
        mlst_datajson = os.path.join(mlst_outdir, "data.json")
        f = open(mlst_datajson, "a+")
        f.close()
        # mlst_cmd="sudo docker run --rm -it \-v {}:/database \-v {}:/workdir \mlst -i {} -o {} -s {}".format(MLST_DB,current_path,input,mlst_outdir,mlst_organism)
        mlst_cmd = "docker run --rm -it \-v {}:/databases \-v {}:/workdir \mlst -i {} -o {} -s {}".format(MLST_DB,
                                                                                                          current_path,
                                                                                                          input,
                                                                                                          mlst_outdir,
                                                                                                          mlst_organism)
        print(mlst_cmd, "\n")
        mlst, err = utils_.run_cmd3(mlst_cmd)
        with open(logpath, "a+") as f:
            if mlst.returncode != 0:
                print(mlst.stdout.readline())
                f.write(err + "\n")
                sys.exit()
            else:
                f.write("mlst is ok\n")
        step += 1
        # time
        with open("./ana_time.csv", "a+") as f:
            fieldnames = ["func", "time"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow({"func": "{} mlst".format(inId), "time": str(time.time() - step1_time)})
    else:
        print("**********       mlst was running.      **********\n next step\n")

    # run plasmidfinder
    if step < 2:
        step2_time = time.time()
        print("STEP{}\n".format(step + 1))
        print("********** Now plasmidfinder analysis running. **********\n")
        PLASMID_DB = "/data/usrhome/LabSSLin/user30/Desktop/SRA_Analysis/plasmidfinder_db"
        # plas_outdir=os.path.join(outdir,"plasmidfinder")
        plas_outdir = os.path.join(relative_path2, "plasmidfinder")
        utils_.mkdir_join(plas_outdir)
        # plas_cmd="sudo docker run --rm -it \-v {}:/databases \-v {}:/workdir \plasmidfinder -i {} -o {}".format(PLASMID_DB,current_path,input,plas_outdir)
        plas_cmd = "docker run --rm -it \-v {}:/databases \-v {}:/workdir \plasmidfinder -i {} -o {}".format(
            PLASMID_DB, current_path, input, plas_outdir)
        print(plas_cmd, "\n")
        plas = run_cmd(plas_cmd)
        with open(logpath, "a+") as f:
            if plas.returncode != 0:
                # print(mlst.stdout.readline())

                sys.exit()
            else:
                f.write("plasmidfinder is ok\n")
        step += 1
        # time
        with open("./ana_time.csv", "a+") as f:
            fieldnames = ["func", "time"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow({"func": "{} plasmidfinder".format(inId), "time": str(time.time() - step2_time)})
    else:
        print("********** plasmidfinder was running. **********\n next step\n")

    # run amrfinder
    if step < 3:
        step3_time = time.time()
        print("STEP{}\n".format(step + 1))
        print("********** Now amrfinder analysis running. **********\n")
        # amr_outdir=os.path.join(outdir,"amrfinder")
        amr_outdir = os.path.join(relative_path2, "amrfinder")
        utils_.mkdir_join(amr_outdir)
        amr_outdir = os.path.join(amr_outdir, "amrout.tsv")
        # amr_cmd="amrfinder -n {} -o {} -O {}".format(input,amr_outdir,amr_organism)
        # amr_cmd = "/data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/amrfinder/amrfinder -n {} -o {} -O {}".format(input, amr_outdir, amr_organism)
        amr_cmd = "amrfinder -n {} -o {} -O {}".format(input, amr_outdir, amr_organism)
        print(amr_cmd, "\n")
        amr = run_cmd(amr_cmd)
        with open(logpath, "a+") as f:
            if amr.returncode != 0:
                # print (mlst.stdout.readline())
                sys.exit()
            else:
                f.write("amr is ok\n")
        step += 1
        # time
        with open("./ana_time.csv", "a+") as f:
            fieldnames = ["func", "time"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow({"func": "{} amr".format(inId), "time": str(time.time() - step3_time)})
    else:
        print("**********   amrfinderr was running.   **********\n next step\n")

    # run sistr
    if step < 4:
        step4_time = time.time()
        print("STEP{}\n".format(step + 1))
        print("********** Now sistr analysis running. **********")
        # sistr_outdir=os.path.join(outdir, "sistr")
        sistr_outdir = os.path.join(relative_path2, "sistr")
        utils_.mkdir_join(sistr_outdir)
        sistr_outdir = os.path.join(sistr_outdir, "sistr_out")
        input_list = input.split("/")
        input_name = input_list[len(input_list) - 1]
        print("name: ", input_name, "\n")
        sistr_cmd = "sistr -i {} {} -f csv -o {} -m".format(input, input_name, sistr_outdir)
        print(sistr_cmd, "\n")

        sistr = run_cmd(sistr_cmd)
        with open(logpath, "a+") as f:
            if sistr.returncode != 0:
                # print(mlst.stdout.readline())
                sys.exit()
            else:
                f.write("sistr is ok\n")
        step += 1
        # time
        with open("./ana_time.csv", "a+") as f:
            fieldnames = ["func", "time"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow({"func": "{} sistr".format(inId), "time": str(time.time() - step4_time)})
    else:
        print("********** sistr was running. **********\n next step\n")

    ########################
    ########################
    # save data in analysis_final.csv and update DB

    # read mlst 'Sequence Type'
    # mlst_file=os.path.join(outdir, "mlst/results.txt")
    mlst_file = os.path.join(relative_path2, "mlst/results.txt")
    with open(mlst_file, "r") as f:
        data = f.readlines()
        print(data[6])
        # Sequence Type: 11
        sequenceType = data[6].split(" ")
        sequenceType = sequenceType[len(sequenceType) - 1].strip("\n")

        print(sequenceType)

    # read plasmidfinder 'gene'
    plas_file = os.path.join(relative_path2, "plasmidfinder/results_tab.tsv")
    # plas_file=os.path.join(outdir,"plasmidfinder/results_tab.tsv")
    # plas_file = os.path.join(outdir, "./plastest/5524p/results_tab.tsv")
    pladf = pd.read_table(plas_file, sep='\t')
    # print(df)
    pladf = pd.DataFrame(pladf)
    print(pladf)
    print(pladf.columns)
    print(pladf.Plasmid)
    plist = list(pladf.Plasmid)
    plas_format = ""

    for x in range(0, len(plist)):
        plas_format += plist[x]
        if x < len(plist) - 1:
            plas_format += ","
    print(plas_format)

    # read amrfinder 'Gene symbol', 'subtype'
    ##add amr "Point"
    # amr_file = os.path.join(outdir, "amrfinder/amrout.tsv")
    amr_file = os.path.join(relative_path2, "amrfinder/amrout.tsv")
    # plas_file = os.path.join(outdir, "./plastest/5524p/results_tab.tsv")
    amrdf = pd.read_table(amr_file, sep='\t')
    # print(df)
    amrdf = pd.DataFrame(amrdf)
    print(amrdf)
    print(amrdf.columns)
    # replace ' ' as '_'
    amrdf.columns = ['Protein_identifier', 'Contig_id', 'Start', 'Stop', 'Strand',
                     'Gene_symbol', 'Sequence_name', 'Scope', 'Element_type',
                     'Element_subtype', 'Class', 'Subclass', 'Method', 'Target_length',
                     'Reference_sequence_length', 'Coverage_of_reference_sequence',
                     'Identity_to_reference_sequence', 'Alignment_length',
                     'Accession_of_closest sequence', 'Name_of_closest sequence', 'HMM_id',
                     'HMM_description']

    print(amrdf.columns)
    # if list is [], show NaN
    print(amrdf.Gene_symbol)
    print(amrdf.Element_subtype)

    amr_format = ""
    point_format = ""
    sym_list = list(amrdf.Gene_symbol)
    sub_list = list(amrdf.Element_subtype)
    for i in range(0, len(sym_list)):
        if sub_list[i] == "AMR":
            if len(point_format) > 0:
                amr_format += ","
            amr_format += sym_list[i]
        elif sub_list[i] == "POINT":
            if len(point_format) > 0:
                point_format += ","
            point_format += sym_list[i]

    print(amr_format)
    print(point_format)

    # if len(aalist)!=0:
    #    amr_format += ","
    # for x in range(0,len(aalist)):
    #    amr_format+=aalist[x]
    #    amr_sub +=alist[x]
    #    if x != len(aalist)-1:
    #        amr_format += ","
    #        amr_sub+=","

    # read sistr 'serovar'
    # sistr_file = os.path.join(outdir, "sistr/sistr_out.csv")
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
    # dict={'Accession': pd.Series(input for a in range(0,3)),
    #      'mlst':sequenceType,
    # }
    in_abspath = input.replace(".", current_path)

    # dict = {'Accession': input,
    #        'mlst':sequenceType,
    #        'plasmidfinder':pladf.Plasmid,
    #        'amr_gane':amrdf.Gene_symbol,
    #        'amr_subtype':amrdf.Element_subtype,
    #        'sistr':sistrdf.serovar
    #        }

    dict = {'Accession': inId.split("_")[0],
            'MLST': sequenceType,
            'AMR': amr_format,
            'Point': point_format,
            'Serotype': sistrdf.serovar,
            'Inc Type': plas_format
            }

    finaldf = pd.DataFrame(dict)
    print(finaldf)
    finalfile = os.path.join(origin_outdir, "analysis_final.csv")
    # if os.path.isfile(finalfile):
    # beforedf=pd.read_csv(finalfile)
    # beforedf=pd.DataFrame(beforedf)
    # finaldf=pd.concat([beforedf,finaldf])
    # print("merge df\n")
    finaldf.to_csv(finalfile, mode='a+', header=False)
    # finaldf.to_csv(finalfile)

    # after run all state, save ID in "Anackeck.log" and remove ./analysis
    with open(check, "a+") as f:
        f.write("Run {} is ok.\n".format(inId))


def SRA_Analysis(x):
    SRA_start=time.time()
    try:
        Download(x)
        Assembled(x)
        #####
        genome = os.path.join(ass_dir, "{}_contig.fa".format(x))
        targetPath=QualityCheck(genome)
        print("targetPAth = {}\n######\n".format(targetPath.encode("utf-8").decode()))
        target_ = targetPath.replace(current_path, ".")
        Analysis(target_,new_outdir)
        print("Run {} is ok\n".format(x))
    except Exception as e:
        error_class = e.__class__.__name__  # 取得錯誤類型
        detail = e.args[0]  # 取得詳細內容
        cl, exc, tb = sys.exc_info()  # 取得Call Stack
        lastCallStack = traceback.extract_tb(tb)[-1]  # 取得Call Stack的最後一筆資料
        fileName = lastCallStack[0]  # 取得發生的檔案名稱
        lineNum = lastCallStack[1]  # 取得發生的行號
        funcName = lastCallStack[2]  # 取得發生的函數名稱
        errMsg = "File \"{}\", line {}, in {}: [{}] {}".format(fileName, lineNum, funcName, error_class, detail)
        print(errMsg)
        with open("./SRA_run_error.txt", "a+") as f:
            f.write("{} :\n{}\n".format(x, errMsg))
        sys.exit(e)
    with open("./threads_time.csv", "a+") as f:
        fieldnames = ["func", "time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func": "{}".format(x), "time": str(time.time() - SRA_start)})
    #######
    with open(check_log,"a+") as f:
        f.write("Run {} is ok.\n".format(x))
    return 0
if __name__ == '__main__':
    start=time.time()
    Month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")
    ########################
    ## read SRAsetting.txt
    utils_.progress_bar("read SRAsetting.txt")
    setting_path = os.path.join(current_path, "SRAsettings.txt")
    with open(setting_path, "r") as f:
        setList = f.readlines()

    print(setList)
    i = 0
    for line in setList:
        line = line.strip("\n")
        line_ = line.split("=")
        if line != "" and len(line_) == 2:
            print(line_)
            print("line{}. {}:{}\n".format(i, line_[0], line_[1]))
        i += 1

    thread = setList[1].strip("\n").split("=")[1]
    gsize = setList[4].strip("\n").split("=")[1]
    n = setList[7].strip("\n").split("=")[1]
    outdir = setList[10].strip("\n").split("=")[1]
    utils_.mkdir_join(outdir)
    ref_dir = setList[12].strip("\n").split("=")[1]
    buscoDB = setList[13].strip("\n").split("=")[1]
    buscoMode = setList[14].strip("\n").split("=")[1]
    mlstS = setList[16].strip("\n").split("=")[1]
    amrS = setList[17].strip("\n").split("=")[1]
    thread=4

    #####################
    for mon in range(7, 8):
        for d in range(1, 4):
            pattern = "salmonella enterica[ORGN] AND illumina[PLAT] AND wgs[STRA] AND genomic[SRC] AND paired[LAY]"
            ds = time.time()

            date = datetime.date(2020, mon + 1, d).strftime("%Y/%m/%d")
            #temp="{}/{}/{}".format(str(2020),str(mon+1),str(d))
            ######
            pdat = date.replace("/", "")
            new_outdir = os.path.join(outdir, pdat)
            utils_.mkdir_join(new_outdir)
            print("output: {}\n".format(new_outdir))

            #Downloadcheck_log = os.path.join(new_outdir, "Downloadcheck.log")
            check_log =os.path.join(new_outdir,"Analysischeck.log")
            # commit
            # run_cmd2("touch {}".format("check.log"))
            #myfile = Path(Downloadcheck_log)
            #myfile.touch(exist_ok=True)
            #f = open(Downloadcheck_log, 'a+')
            #t = str(datetime.datetime.now()).split(".")[0]
            #f.write(t)
            #f.write("\n")
            #f.close()

            myfile2 = Path(check_log)
            myfile2.touch(exist_ok=True)
            with open(check_log,"a+") as f:
                f.write(str(datetime.datetime.now()).split(".")[0])
                f.write("\n")

            # mkdir_join(log_dir)
            # check_log = os.path.join(log_dir, "check.log")

            # print(date)

            pattern, count = utils_.count_egquery(pattern, date, date)
            print("pattern: {}\ncount: {}\n".format(pattern, count))

            i_e_ = time.time()
            idlist = utils_.IdList_esearch(pattern, 'sra', count)

            print(idlist)

            runinfo = utils_.Get_RunInfo(idlist)
            # progress_bar("get SRAfile name List stored in run_list")
            run_list = list(runinfo['Run'])  # get SRAfile nameList stored in run_list
            print("runinfo: {}\n run_list: {}\n".format(runinfo, run_list))

            sra_dir = os.path.join(new_outdir, "sra")  # .sra file
            utils_.mkdir_join(sra_dir)
            ass_dir = os.path.join(new_outdir, "Assembled")
            utils_.mkdir_join(ass_dir)
            fastq_dir = os.path.join(new_outdir, 'fastq')
            assemble_dir = os.path.join(new_outdir, "assembly_result")

            read_log_ = time.time()
            f = open(check_log, 'r+')
            line = f.readlines()
            print("check log :{}\n".format(line))
            f.close()

            for s in line:
                print("{}\n".format(s))
            finish = list(filter(lambda x: len(x.split(" ")) >= 4, line))
            finish_run = list(map(lambda x: x.split(" ")[1], finish))
            need_run = list(filter(lambda x: x not in finish_run, run_list))
            print("finish: {}\nfinish_run: {}\nneed_run".format(finish, finish_run, need_run))
            print("finish length: {}\nfinish_run length: {}\nneed_run length: ".format(len(finish), len(finish_run),
                                                                                       len(need_run)))
            print("Toal", len(need_run), "sra runs need to downlaod.")

            num = len(finish_run)
            progress_list = []
            prog_num = 0
            finish_num=0
            finish_num=len(finish_run)
            pool=multiprocessing.Pool(processes=1)
            for k in need_run:
                print("########### hello %d ############\n"%prog_num)
                print("########## {}/{} ###########".format(finish_num,count))
                pool.apply_async(SRA_Analysis, (k,))
                progress_list.append(multiprocessing.Process(target=SRA_Analysis, args=(k,)))
                #progress_list[prog_num].start()
                prog_num += 1
                finish_num+=1
            pool.close()
            pool.join()
            #for i in range(prog_num):
            #    progress_list[i].join()


            with open("./Automate_check.log", "a+") as f:
                f.write("{}:{}:{}\n".format(date, time.time() - ds, time.time() - start))
            print("Download all {} ".format(date),'Done,total cost', time.time() - ds, 'secs')
            time.sleep(5)
    print('Done,total cost', time.time() - start, 'secs')
    ##########





