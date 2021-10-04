import datetime
import sys
import time
import traceback

import utils_
def main():
    start = time.time()
    cmd="python3 SRA_Analysis.py --PDAT "
    Month=[31,28,31,30,31,30,31,31,30,31,30,31]
    try:
        for x in range(0, 12):
            for d in range(1, Month[x] + 1):
                ds = time.time()
                c = datetime.datetime(2020, x + 1, d)
                tmp = c.strftime("%Y/%m/%d")
                cmd2 = cmd + tmp
                print(cmd2)
                utils_.run_cmd3(cmd2)
                with open("./Automate_check.log", "a+") as f:
                    f.write("{}:{}:{}\n".format(tmp, time.time() - ds, time.time() - start))
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
        sys.exit(e)

    print('Done,total cost', time.time() - start, 'secs')
    return 0


if __name__ == '__main__':
    main()