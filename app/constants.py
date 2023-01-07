# SQL Tables
INPUT_TABLE = "Input"
EXECUTION_TABLE = "Execution"
OMP_TABLE = "Omp"
INTEGRATOR_CONFIG_TABLE = "Integrator_Config"
LOG_FILES_TABLE = "Log_Files"

INPUT_DATA = [
    {
        "input_id": "/ref_data/res_gri3.0.csv",
        "mechanism": "gri30",
        "nsp": 53,
        "systems": 55076
    },
    {
        "input_id": "/ref_data/res_gri_32.csv",
        "mechanism": "gri30",
        "nsp": 53,
        "systems": 32
    },
]

# Components and configuration
READERS = ["csvReader", "testReader"]
WRITERS = ["csvWriter", "testWriter"]
LOGGERS = ["ConsoleLogger", "FileLogger"]
LOGGER_LEVELS = ["DEBUG", "INFO"]
INTEGRATORS = [
    "CanteraIntegrator",
    "CanteraIntegratorOMP",
    "CVodeIntegrator",
    "CVodeIntegratorOMP",
    "CVodeIntegratorGPU"
]
OMP_SCHEDULES = [
    "static",
    "dynamic",
    "guided",
    "auto"
]
