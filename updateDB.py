import pandas as pd
import pymysql
def main():
    print("update data to Databases\n")
    db_settings={
        "host":"127.0.0.1",
        "port":3306,
        "user":"root",
        "password":"tumvgk01",
        "db":"SRA_Analysis",
        "charset":"utf8"
    }
    file="/data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/SRA-AutomatedAnalysis/finaltest.csv"
    df = pd.read_csv(file)
    # print(df)
    df = pd.DataFrame(df)
    print(df)
    df=df.fillna(value="NAN")
    print(df)
    #for i in range(len(df)):
    #    print(str(df.loc[i, "Accession"]) + " " + str(df.loc[i, "mlst"]))


        #conn = pymysql.connect(**db_settings)
    conn=pymysql.connect(host="127.0.0.1",user="root",password="tumvgk01",database="SRA_Analysis",port=3306)
    cursor=conn.cursor()
    insert = "INSERT INTO `Final`(`Accession`,`MLST`,`AMR`,`Serotype`,`Inc_Type`) VALUES({},{},{},{},{})".format(
        "22",222,"222",'jjj',"aaa")
    try:
        cursor.execute(insert)
        conn.commit()
    except Exception as e:
        print("ffff")
        conn.rollback()
        print(e)
    #insertSRA = "INSERT INTO SRA(Genome) VALUES(%s);"
    #insert = "INSERT INTO Final(Accession,MLST,AMR,Serotype,Inc_Type) VALUES(%s,%s,%s,%s,%s);"
    #for i in range(0,len(df)-1):
        #print(str(df.loc[i,"Accession"])+" "+str(df.loc[i,"mlst"])+" "+str(df.loc[i,"plasmidfinder"])+" "+str(df.loc[i,"amr_gane"])+" "+str(df.loc[i,"sistr"]))
    #    insert = "INSERT INTO `Final`(`Accession`,`MLST`,`AMR`,`Point`,`Serotype`,`Inc_Type`) VALUES ({},{},{},NAN,{},0)".format(
    #        str(df.loc[i, "Accession"]), int(df.loc[i, "mlst"]), str(df.loc[i, "amr_gane"]).replace(",","|"), str(df.loc[i,"sistr"]).replace(",","|")
     #       )
    #    try:
    #        cursor.execute(insert)
    #        conn.commit()
   #         #cursor.execute(insert, (str(df.loc[i,"Accession"]),str(df.loc[i,"mlst"]),str(df.loc[i,"amr_gane"]),str(df.loc["sistr"]),str(df.loc[i,"plasmidfinder"])))
   #     except Exception as e:
    #        print("ffff")
    #        conn.rollback()
    #        print(e)
    conn.close()






    return 0

if __name__ == '__main__':
    main()