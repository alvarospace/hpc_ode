import sqlite3
import argparse
from pathlib import Path

from constants import INPUT_DATA

DB_RELATIVE_LOCATION = "/app/database.py"

class DataBase:
    def __init__(self, path: str):
        self._path = path
        self._connection = sqlite3.connect(self._path, detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES)
        self._connection.execute("PRAGMA foreign_keys = ON;")
        self._cursor = self._connection.cursor()

    def query(self, sql: str) -> list:
        self._cursor.execute(sql)
        return self._cursor.fetchall()

    # TODO: insert to database
    def insert(self):
        pass

    def get_all_from_table(self, table: str) -> list | None:
        self._cursor.execute("SELECT * FROM :table", {"table": table})
        return self._cursor.fetchall()

    def close(self):
        self._connection.close()
        self._connection = None
        self._cursor = None

    def create_tables(self):
        self.create_input_table()
        self.create_execution_table()
        self.create_integrator_config_table()
        self.create_omp_table()
        self.create_log_files_table()

    def create_input_table(self) -> None:
        self._cursor.execute("""CREATE TABLE input (
            input_id         TEXT       NOT NULL    PRIMARY KEY, -- path to the file
            mechanism        TEXT       NOT NULL, -- example: gri30 (see resources for more examples)
            nsp              INTEGER    NOT NULL,
            systems          INTEGER    NOT NULL
        );
        """)

    def create_execution_table(self) -> None:
        self._cursor.execute("""CREATE TABLE execution (
            execution_id     INTEGER    NOT NULL    PRIMARY KEY,
            input_id         TEXT       NOT NULL,
            reader           TEXT       NOT NULL,
            writer           TEXT       NOT NULL,
            integrator       TEXT       NOT NULL,
            logger           TEXT       NOT NULL,
            output           TEXT       NOT NULL, -- path to the output file
            read_time        REAL       NOT NULL,
            integration_time REAL       NOT NULL,
            write_time       REAL       NOT NULL,
            total_time       REAL       NOT NULL,
            date             TIMESTAMP  NOT NULL,
            FOREIGN KEY (input_id) REFERENCES input( input_id )
        );
        """)

    def create_integrator_config_table(self) -> None:
        self._cursor.execute("""CREATE TABLE integrator_config (
            integrator_config_id  INTEGER    NOT NULL    PRIMARY KEY,
            execution_id          INTEGER    NOT NULL,
            reltol                REAL       NOT NULL,
            abstol                REAL       NOT NULL,
            pressure              REAL       NOT NULL,
            dt                    REAL       NOT NULL,
            FOREIGN KEY (execution_id) REFERENCES execution( execution_id )
        );
        """)

    def create_omp_table(self) -> None:
        self._cursor.execute("""CREATE TABLE omp (
            omp_id           INTEGER    NOT NULL    PRIMARY KEY,
            execution_id     INTEGER    NOT NULL,
            cpus             INTEGER    NOT NULL,
            schedule         TEXT       NOT NULL,
            chunk            INTEGER    NOT NULL,
            FOREIGN KEY (execution_id) REFERENCES execution( execution_id )
        );
        """)

    def create_log_files_table(self) -> None:
        self._cursor.execute("""CREATE TABLE log_files (
            integrator_config_id  INTEGER    NOT NULL    PRIMARY KEY,
            execution_id          INTEGER    NOT NULL,
            log_file              BLOB       NOT NULL,
            FOREIGN KEY (execution_id) REFERENCES execution( execution_id )
        );
        """)

    def init_input_table(self) -> None:
        for entry in INPUT_DATA:
            self._cursor.execute("""INSERT INTO input (input_id, mechanism, nsp, systems)
            VALUES (:input_id, :mechanism, :nsp, :systems);
            """, entry)
        self._connection.commit()


def main(args: argparse.Namespace):
    # Add the parent directory of the repository to the 
    # input file of each entry
    repo_directory = __file__.replace(DB_RELATIVE_LOCATION, "")
    global INPUT_DATA
    for entry in INPUT_DATA:
        entry["input_id"] = repo_directory + entry["input_id"]

    path = Path(args.path)
    if path.is_dir() is False:
        path.mkdir(parents=True)
    path = path.joinpath("ODEIntegratorDB.db")
    
    if args.create:
        db = DataBase(str(path))
        db.create_tables()
        db.init_input_table()
        db.close()

    if args.delete:
        path.unlink()

if __name__ == "__main__":
    # Arguments parser
    parser = argparse.ArgumentParser("Database manager")
    parser.add_argument("-p", "--path", required=True, type=str,
                        help="Host path where database is or where it will be created")
    parser.add_argument("-c", "--create", action="store_true", default=False,
                        help= "Flag to create new database or create tables in a empty database")
    parser.add_argument("-d", "--delete", action="store_true", default=False,
                        help="Flag to delete database")

    args = parser.parse_args()
    main(args)

    
