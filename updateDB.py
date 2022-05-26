import pandas as pd
import sys
import traceback
import pymysql
from pymysql.converters import escape_string


def main():
    print("update data to Databases\n")

    file = "./analysis_final08.csv"
    df = pd.read_csv(file)
    # print(df)
    df = pd.DataFrame(df)
    # df.columns = ["","Accession", "mlst", "AMR", "Point", "amr_gane", "sistr"]
    print(df)
    df = df.fillna(value="NAN")
    print(df)
    # for i in range(len(df)):
    #    print(str(df.loc[i, "Accession"]) + " " + str(df.loc[i, "mlst"]))

    # conn = pymysql.connect(**db_settings)
    conn = pymysql.connect(host="127.0.0.1", user="root", password="server@ntutrag", database="SRA_analysis", port=3306)
    cursor = conn.cursor()
    # insert = "INSERT INTO `Final`(`Accession`, `MLST`, `AMR`, `Point`, `Serotype`, `Inc_Type`) VALUES ({},{},{},{},{},{}) ".format("12",111,"111","111","111","1111")
    # try:
    #    cursor.execute(insert)
    #    conn.commit()
    # except Exception as e:
    #    print("ffff")
    #   conn.rollback()
    #   print(e)
    # insertSRA = "INSERT INTO web(Genome) VALUES(%s);"
    # insert = "INSERT INTO Final(Accession,MLST,AMR,Serotype,Inc_Type) VALUES(%s,%s,%s,%s,%s);"
    for i in range(0, len(df)):
        ## modify seq
        ST = escape_string(str(df.loc[i, "MLST"]))
        serotype_ = escape_string(str(df.loc[i, "Serotype"]))
        inctype_ = escape_string(str(df.loc[i, "IncType"]))
        AMR_ = escape_string(str(df.loc[i, "AMR"]))
        point_ = escape_string(str(df.loc[i, "Point"]))

        print(str(df.loc[i, "Accession"]) + " " + str(df.loc[i, "MLST"]) + " " + str(df.loc[i, "AMR"]) + str(
            df.loc[i, "Point"]) + " " + str(df.loc[i, "Serotype"]) + " " + str(df.loc[i, "IncType"]))
        print(str(df.loc[i, "Accession"]) + " " + ST + " " + serotype_ + " " + inctype_ + " " + AMR_ + " " + point_)
        insert = "INSERT INTO `Final`(`Accession`, `ST`, `Serotype`, `IncType`, `AMR`, `Point`) VALUES ('{}','{}','{}','{}','{}','{}') ".format(
            str(df.loc[i, "Accession"]), str(ST), serotype_, inctype_, AMR_, point_)
        print(insert)
        try:
            cursor.execute(insert)
            conn.commit()
            # cursor.execute(insert, (str(df.loc[i,"Accession"]),str(df.loc[i,"mlst"]),str(df.loc[i,"amr_gane"]),str(df.loc["sistr"]),str(df.loc[i,"plasmidfinder"])))
        except Exception as e:
            print("Failed\n")
            conn.rollback()
            error_class = e.__class__.__name__  # 取得錯誤類型
            detail = e.args[0]  # 取得詳細內容
            cl, exc, tb = sys.exc_info()  # 取得Call Stack
            lastCallStack = traceback.extract_tb(tb)[-1]  # 取得Call Stack的最後一筆資料
            fileName = lastCallStack[0]  # 取得發生的檔案名稱
            lineNum = lastCallStack[1]  # 取得發生的行號
            funcName = lastCallStack[2]  # 取得發生的函數名稱
            errMsg = "File \"{}\", line {}, in {}: [{}] {}".format(fileName, lineNum, funcName, error_class,
                                                                   detail)
            if "IntegrityError" not in error_class:
                print(errMsg)

                sys.exit(e)
    conn.close()
    print("Done\n")

    return 0


if __name__ == '__main__':
    main()