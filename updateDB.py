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

    try:
        conn=pymysql.connect(**db_settings)

        with conn.cursor() as cursor:
            insertSRA="INSERT INTO SRA(Genome) VALUES({})".format()

    except Exception as e:
        print(e)


    return 0

if __name__ == '__main__':
    main()