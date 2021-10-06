import csv


def main():

    with open("./test1006.csv","a+")as f:
        fieldnames=["func","time"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({"func":"ffff","time":str(123)})

    with open("./test1006.csv","r")as f:
        k=f.readlines()
        print(k)

if __name__ == '__main__':
    main()