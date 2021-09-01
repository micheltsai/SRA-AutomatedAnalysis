from __future__ import print_function
import myfunction_ as func_
import pandas as pd
from Bio import Entrez
import datetime
import glob
import xml.etree.cElementTree as ET
from io import StringIO
from tempfile import TemporaryDirectory
Entrez.email = 'ann850324@gmail.com'


###### new add ######
def IdList_esearch(term, db, count):
    handle = Entrez.esearch(term = term, db = db, retmax = count)#不設定retmax的話只有20筆資料
    func_.progress_bar("read and stored IdList")
    d = Entrez.read(handle)
    return d['IdList']

def Get_RunInfo(idlist):
    c = len(idlist)
    print('Getting run_info table by idlist...')
    func_.progress_bar("read and stored RunInfo")
    if c >= 10000:
        print("over 10000 results")
        df_all = pd.DataFrame()
        #設定k為list,由0開始,最大值為c, 間隔10000
        k = list(range(0,c,10000))
        k.append(len(idlist))
        print(k)
        for i in range(1,len(k)):
            s = datetime.time.time()
            start = k[i-1]
            end = k[i]
            #下載GenBank records
            handle = Entrez.efetch(db = 'sra', id = idlist[start:end],rettype = 'runinfo',retmode = 'csv')
            #查看原始的Genbank文件
            d = handle.read()
            #讀檔
            df = pd.read_csv(StringIO(d))
            df = df[df['Run'] != 'Run']
            df_all = pd.concat([df_all,df])
            print("finish", i - 1,"round,cost", datetime.time.time() - s, "secs")
            print("counts_of_run:",len(df))
            datetime.time.sleep(2)
    else:
        handle = Entrez.efetch(db = 'sra', id = idlist,rettype = 'runinfo',retmode = 'csv')
        d = handle.read()
        df = pd.read_csv(StringIO(d))
        df_all = df[df['Run'] != 'Run']
    return df_all
