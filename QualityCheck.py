from __future__ import print_function
import argparse
import os
import shlex
import subprocess
import sys
import time

# -*- coding: utf-8 -*-
def progress_bar(Category):
    for i in range(1, 101):
        print("\r{}: ".format(Category),end="")
        print("[{}] {}%".format("*" * (i // 2), i), end="")
        sys.stdout.flush()
        time.sleep(0.02)
    print ("\n")


def run_cmd2(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True)
    return p


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

def getRefListPath(refSeqPath,outdir):
    print("getRefListPath:\n")
    print("refSeqPath: " + refSeqPath + "\n")
    refListPath = os.path.join(outdir, 'ref.txt')
    run_cmd2("find {}>{}".format(refSeqPath,refListPath))
    print("refListPath: "+refListPath+"\n")
    #run_cmd2("cat {}".format(refListPath))
    return refListPath

def getQenomeListPath(qenome_Path,outdir):
    print("getQenomeListPath:\n")
    print("refSeqPath: " + qenome_Path + "\n")
    qenListPath = os.path.join(outdir, 'qen.txt')
    run_cmd2("find {}>{}".format(qenome_Path,qenListPath))
    print("qenListPath: "+qenListPath+"\n")
    return qenListPath

def main():
    start = time.time()
    #get ref_path, qen_path, and outdir
    parser = argparse.ArgumentParser("Analytics")
    parser.add_argument("-r","--ref", required=True, help="Path of reference Sequence files")
    parser.add_argument("-g","--qenome", required=True, help="Path of qenome files")
    parser.add_argument("-o","--outdir", required=True, help="Output folder")

    args = parser.parse_args()

    refPath=getRefListPath(args.ref,args.outdir)
    qenome_Path=getQenomeListPath(args.qenome,args.outdir)
    outdir = args.outdir
    out_txt=os.path.join(outdir, 'out.txt')
    progress_bar("load args")

    #fastANI
    print ("reseq: {}\n, qen: {}\n, outdir: {}\n,out_txt: {}".format(refPath, qenome_Path, outdir,out_txt))
    print("fastANI start.")
    progress_bar("fastANI excuting")
    #fasani_=run_cmd("/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI -h")
    fastani_="/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI --rl {} --ql {} -o {}".format(refPath,qenome_Path,out_txt)
    print(fastani_+"\n")
    run_cmd(fastani_)
    print("fastANI done.\n")

    #ANI>=95, get out.txt line 1
    print ("----------------fastANI end----------------\nget maxANI refseqPath:")
    f = open(out_txt, 'r')
    targetPath=f.readline().split("\t")[1]
    print(targetPath,"\n")
    f.close()






    print('Done,total cost', time.time() - start, 'secs\n')

if __name__ == '__main__':
    main()
