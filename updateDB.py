import pandas as pd
import pymysql
def main():
    print("update data to Databases\n")
    db_settings={
        "host":"127.0.0.1",
        "port":3306,
        "user":"user",
        "password":"tumvgk01",
        "db":"SRA_Analysis",
        "charset":"utf8"
    }
    file="/data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/SRA-AutomatedAnalysis/finaltest.csv"
    df = pd.read_csv(file)
    # print(df)
    df = pd.DataFrame("df")
    print(df)



    try:
        conn = pymysql.connect(**db_settings)

        with conn.cursor() as cursor:
            #insertSRA = "INSERT INTO SRA(Genome) VALUES(%s);"
            insert = "INSERT INTO Final(Accession,MLST,AMR,Serotype,Inc_Type) VALUES(%s,%s,%s,%s,%s);"
            for i in range(len(df)):
                print(df.loc[i,"Accession"]+" "+df.loc[i,"MLST"])
                #cursor.execute(insert, ())


    except Exception as e:
        print(e)


    return 0

if __name__ == '__main__':
    main()