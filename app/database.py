import sqlite3

INPUT_TABLE = "input"
EXECUTION_TABLE = "execution"
OMP_TABLE = "omp"
INTEGRATOR_CONFIG_TABLE = "integrator_config"
LOG_FILES_TABLE = "log_files"

INPUT_DATA = [
    {
        "input_id": "/home/almousa/TFM/hpc_cvode/ref_data/res_gri3.0.csv",
        "mechanism": "gri30",
        "nsp": 53,
        "systems": 55076
    },
    {
        "input_id": "/home/almousa/TFM/hpc_cvode/ref_data/res_gri_32.csv",
        "mechanism": "gri30",
        "nsp": 53,
        "systems": 32
    },
]

class DataBase:
    def __init__(self, path: str):
        self._path = path
        self._connection = sqlite3.connect(self._path, detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES)
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
        self._cursor.execute("""CREATE TABLE :input_table (
            input_id         TEXT       NOT NULL    PRIMARY KEY, -- path to the file
            mechanism        TEXT       NOT NULL, -- example: gri30 (see resources for more examples)
            nsp              INTEGER    NOT NULL,
            systems          INTEGER    NOT NULL,
        );
        """, {"input_table": INPUT_TABLE})

    def create_execution_table(self) -> None:
        self._cursor.execute("""CREATE TABLE :execution_table (
            execution_id     INTEGER    NOT NULL    PRIMARY KEY,
            input_id         TEXT       NOT NULL    REFERENCES :input_table( input_id ),
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
        );
        """, {"execution_table": EXECUTION_TABLE, "input_table": INPUT_TABLE})

    def create_integrator_config_table(self) -> None:
        self._cursor.execute("""CREATE TABLE :integrator_config_table (
            integrator_config_id  INTEGER    NOT NULL    PRIMARY KEY,
            execution_id          INTEGER    NOT NULL    REFERENCES :execution_table( execution_id ),
            reltol                REAL       NOT NULL,
            abstol                REAL       NOT NULL,
            pressure              REAL       NOT NULL,
            dt                    REAL       NOT NULL
        );
        """, {"integrator_config_table": INTEGRATOR_CONFIG_TABLE, "execution_table": EXECUTION_TABLE})

    def create_omp_table(self) -> None:
        self._cursor.execute("""CREATE TABLE :omp_table (
            omp_id           INTEGER    NOT NULL    PRIMARY KEY,
            execution_id     INTEGER    NOT NULL    REFERENCES :execution_table( execution_id ),
            cpus             INTEGER    NOT NULL,
            schedule         TEXT       NOT NULL,
            chunk            INTEGER    NOT NULL
        );
        """, {"omp_table": OMP_TABLE, "execution_table": EXECUTION_TABLE})

    def create_log_files_table(self) -> None:
        self._cursor.execute("""CREATE TABLE :log_files_table (
            integrator_config_id  INTEGER    NOT NULL    PRIMARY KEY,
            execution_id          INTEGER    NOT NULL    REFERENCES :execution_table( execution_id ),
            log_file              BLOB       NOT NULL
        );
        """, {"log_files_table": LOG_FILES_TABLE, "execution_table": EXECUTION_TABLE})
    
            
            
if __name__ == "__main__":
    db = DataBase("../results/database.db")
    db.create_tables()
    db.close()
