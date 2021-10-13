import subprocess
import time
def run_cmd3(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    returncode = p.poll()
    while returncode is None:
        line = p.stdout.readline()
        returncode = p.poll()
        line = line.strip()
        line_ = line.decode().split("\n")
        #print (line.decode(),"\n\n")
        for s in line_:
            print(str("{}".format(s)))
            err=s
    return p,err

def main():
    start=time.time()
    time1=time.time()

    run_cmd3("/data1/usrhome/LabSSLin/linss01/anaconda3/bin/shovill --R1 /data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/SRA-AutomatedAnalysis/SRAtest/20200704/fastq/SRR12145778/R1.fq --R2 /data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/SRA-AutomatedAnalysis/SRAtest/20200704/fastq/SRR12145778/R2.fq --outdir /data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/SRA-AutomatedAnalysis/SRAtest/20200704/assembly_result/SRR12145778 --depth 100 --tmpdir . --cpus 8 --ram 3 --force")
    with open("testCMD.txt","a+")as f:
        f.write("time1: {}\n".format(time.time()-time1))

    time2=time.time()
    run_cmd3("seqtk fqchk -q3 \/data1\/usrhome\/LabSSLin\/linss01\/Desktop\/SRA\-AutoAnalysis\/SRA\-AutomatedAnalysis\/SRAtest\/20200704\/fastq\/SRR12145778\/R1\.fq >/tmp/gkbYYK8NN8 2>&1 | sed 's/^/[seqtk] /' | tee -a shovill.log")
    with open("testCMD.txt","a+")as f:
        f.write("time2: {}\n".format(time.time()-time2))

    time3=time.time()
    run_cmd3("kmc -sm -m1 -t8 -k21 -ci10 \/data1\/usrhome\/LabSSLin\/linss01\/Desktop\/SRA\-AutoAnalysis\/SRA\-AutomatedAnalysis\/SRAtest\/20200704\/fastq\/SRR12145778\/R1\.fq /tmp/ok6feIchfJ/kmc /tmp/ok6feIchfJ 2>&1 | sed 's/^/[kmc] /' | tee -a shovill.log")
    with open("testCMD.txt","a+")as f:
        f.write("time3: {}\n".format(time.time()-time3))

    time4=time.time()
    run_cmd3("ln -sf \/data1\/usrhome\/LabSSLin\/linss01\/Desktop\/SRA\-AutoAnalysis\/SRA\-AutomatedAnalysis\/SRAtest\/20200704\/fastq\/SRR12145778\/R1\.fq R1.fq.gz 2>&1 | sed 's/^/[ln] /' | tee -a shovill.log")
    with open("testCMD.txt","a+")as f:
        f.write("time4: {}\n".format(time.time()-time4))

    time5=time.time()
    run_cmd3("ln -sf \/data1\/usrhome\/LabSSLin\/linss01\/Desktop\/SRA\-AutoAnalysis\/SRA\-AutomatedAnalysis\/SRAtest\/20200704\/fastq\/SRR12145778\/R2\.fq R2.fq.gz 2>&1 | sed 's/^/[ln] /' | tee -a shovill.log")
    with open("testCMD.txt","a+")as f:
        f.write("time5: {}\n".format(time.time()-time5))


    time6=time.time()
    run_cmd3("lighter -od . -r R1.fq.gz -r R2.fq.gz -K 32 4810242 -t 8 -maxcor 1 2>&1 | sed 's/^/[lighter] /' | tee -a shovill.log")
    with open("testCMD.txt","a+")as f:
        f.write("time6: {}\n".format(time.time()-time6))

    time7=time.time()
    run_cmd3("flash -m 20 -M 236 -d . -o flash -z -t 8 R1.cor.fq.gz R2.cor.fq.gz 2>&1 | sed 's/^//' | tee -a shovill.log")
    with open("testCMD.txt","a+")as f:
        f.write("time7: {}\n".format(time.time()-time7))

    time8=time.time()
    run_cmd3("spades.py -1 flash.notCombined_1.fastq.gz -2 flash.notCombined_2.fastq.gz --isolate --threads 8 --memory 3 -o spades --tmp-dir . -k 31,55,79,103,127  --merged flash.extendedFrags.fastq.gz 2>&1 | sed 's/^/[spades] /' | tee -a shovill.log")
    with open("testCMD.txt","a+")as f:
        f.write("time8: {}\n".format(time.time()-time8))

    time9=time.time()
    run_cmd3("bwa index spades.fasta 2>&1 | sed 's/^/[bwa-index] /' | tee -a shovill.log")
    with open("testCMD.txt","a+")as f:
        f.write("time9: {}\n".format(time.time()-time9))

    time10=time.time()
    run_cmd3("samtools faidx spades.fasta 2>&1 | sed 's/^/[faidx] /' | tee -a shovill.log")
    with open("testCMD.txt","a+")as f:
        f.write("time10: {}\n".format(time.time()-time10))

    time11=time.time()
    run_cmd3("(bwa mem -v 3 -x intractg -t 8 spades.fasta R1.fq.gz R2.fq.gz | samclip --ref spades.fasta.fai | samtools sort --threads 2 -m 512m --reference spades.fasta -T . -o shovill.bam) 2>&1 | sed 's/^/[bwa+samtools-sort] /' | tee -a shovill.log")
    with open("testCMD.txt","a+")as f:
        f.write("time11: {}\n".format(time.time()-time11))

    time12=time.time()
    run_cmd3("samtools index shovill.bam 2>&1 | sed 's/^/[samtools-index] /' | tee -a shovill.log")
    with open("testCMD.txt","a+")as f:
        f.write("time12: {}\n".format(time.time()-time12))

    time13=time.time()
    run_cmd3("pilon --genome spades.fasta --frags shovill.bam --minmq 60 --minqual 3 --fix bases --output pilon --threads 8 --changes --mindepth 0.25 2>&1 | sed 's/^/[pilon] /' | tee -a shovill.log")
    with open("testCMD.txt","a+")as f:
        f.write("time13: {}\n".format(time.time()-time13))


    print('Done,total cost', time.time() - start, 'secs')
if __name__ == '__main__':
    main()