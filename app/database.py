import sqlite3
import argparse
from pathlib import Path

from common import get_repo_directory
from constants import INPUT_DATA

class DataBase:
    """Class to connect to the SQLite3 database
    """
    def __init__(self, path: str):
        """Constructor

        Args:
            path (str): absolute or relative path to the sqlite3 database
        """
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

    def add_item_into_input(self, id, mechanism, nsp, systems) -> None:
        row_dict = {
            "input_id": id,
            "mechanism": mechanism,
            "nsp": nsp,
            "systems": systems
        }

        # Check that the id does not already exists
        existing_items = self.query("SELECT * FROM Input")
        for item in existing_items:
            current_id: str = item[0]
            if current_id == id:
                raise Exception("Trying to add into input a item with an input_id that already exists")

        self._cursor.execute("""INSERT INTO Input (input_id, mechanism, nsp, systems)
        VALUES (:input_id, :mechanism, :nsp, :systems)
        """, row_dict)
        self._connection.commit()


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
        self._cursor.execute("""CREATE TABLE Input (
            input_id         TEXT       NOT NULL    PRIMARY KEY, -- path to the file
            mechanism        TEXT       NOT NULL, -- example: gri30 (see resources for more examples)
            nsp              INTEGER    NOT NULL,
            systems          INTEGER    NOT NULL
        );
        """)

    def create_execution_table(self) -> None:
        self._cursor.execute("""CREATE TABLE Execution (
            execution_id     INTEGER    NOT NULL    PRIMARY KEY,
            input_id         TEXT       NOT NULL,
            reader           TEXT       NOT NULL,
            writer           TEXT       NOT NULL,
            integrator       TEXT       NOT NULL,
            logger           TEXT       NOT NULL,
            output           BLOB       NOT NULL,
            read_time        REAL       NOT NULL,
            integration_time REAL       NOT NULL,
            write_time       REAL       NOT NULL,
            total_time       REAL       NOT NULL,
            date             TIMESTAMP  NOT NULL,
            FOREIGN KEY (input_id) REFERENCES Input( input_id )
        );
        """)

    def create_integrator_config_table(self) -> None:
        self._cursor.execute("""CREATE TABLE Integrator_Config (
            integrator_config_id  INTEGER    NOT NULL    PRIMARY KEY,
            execution_id          INTEGER    NOT NULL,
            reltol                REAL       NOT NULL,
            abstol                REAL       NOT NULL,
            pressure              REAL       NOT NULL,
            dt                    REAL       NOT NULL,
            FOREIGN KEY (execution_id) REFERENCES Execution( execution_id )
        );
        """)

    def create_omp_table(self) -> None:
        self._cursor.execute("""CREATE TABLE Omp (
            omp_id           INTEGER    NOT NULL    PRIMARY KEY,
            execution_id     INTEGER    NOT NULL,
            cpus             INTEGER    NOT NULL,
            schedule         TEXT       NOT NULL,
            chunk            INTEGER    NOT NULL,
            FOREIGN KEY (execution_id) REFERENCES Execution( execution_id )
        );
        """)

    def create_log_files_table(self) -> None:
        self._cursor.execute("""CREATE TABLE Log_Files (
            integrator_config_id  INTEGER    NOT NULL    PRIMARY KEY,
            execution_id          INTEGER    NOT NULL,
            log_file              BLOB       NOT NULL,
            FOREIGN KEY (execution_id) REFERENCES Execution( execution_id )
        );
        """)

    def init_input_table(self, input_data) -> None:
        for entry in input_data:
            self._cursor.execute("""INSERT INTO Input (input_id, mechanism, nsp, systems)
            VALUES (:input_id, :mechanism, :nsp, :systems);
            """, entry)
        self._connection.commit()


def main(args: argparse.Namespace):
    # Add the parent directory of the repository to the 
    # input file of each entry
    input_data = INPUT_DATA
    repo_directory = get_repo_directory()
    for entry in input_data:
        entry["input_id"] = repo_directory + entry["input_id"]

    path = Path(args.path)
    if path.is_dir() is False:
        path.mkdir(parents=True)
    path = path.joinpath("ODEIntegratorDB.db")
    
    if args.create:
        db = DataBase(str(path))
        db.create_tables()
        db.init_input_table(input_data)
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

    
