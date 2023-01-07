from typing import Callable

from tkinter import ttk
from tkinter import *
from tkinter import filedialog

import yaml

from constants import (
    READERS,
    WRITERS,
    LOGGERS,
    LOGGER_LEVELS,
    INTEGRATORS,
    OMP_SCHEDULES,
)
from common import get_mechanism_list

PADDING = 20

class ConfigYamlPage:
    def __init__(self, root: Tk):
        self.root = root
        self.root.title("YAML Configuration")
        self.root.geometry("1280x720")
        self.yaml_dict = {}


        self.main_frame = ttk.Frame(self.root, padding=PADDING)
        self.main_frame.bind("<Configure>", self.print)
        self.left_column = ttk.Frame(self.main_frame, padding=PADDING, borderwidth=2, relief="raised")
        self.right_column = ttk.Frame(self.main_frame, padding=PADDING, borderwidth=2, relief="raised")
        
        self.main_frame.columnconfigure(0, weight=1)
        self.main_frame.columnconfigure(1, weight=1)
        self.main_frame.rowconfigure(0, weight=1)
        self.main_frame.rowconfigure(1, weight=1)
        self.left_column.columnconfigure(0, weight=1)
        self.left_column.rowconfigure(0, weight=1)
        self.right_column.columnconfigure(0, weight=1)
        self.right_column.rowconfigure(0, weight=1)

        self.main_frame.grid(column=0, row=0, sticky=(N, S, E, W))
        self.left_column.grid(column=0, row=0, sticky=(N, S, E, W))
        self.right_column.grid(column=1, row=0, sticky=(N, S, E, W))

        self.init_reader(self.left_column)
        self.init_writer(self.left_column)
        self.init_logger(self.left_column)
        self.init_outfolder(self.main_frame)
        self.init_integrator(self.right_column)

    def print(self, e: Event):
        print(e.widget)

    ######### READER LOGIC ###########
    def init_reader(self, frame: ttk.Frame) -> None:
        self.reader_frame = ttk.Labelframe(frame, padding=PADDING/2, text="Reader")
        ###### Important Var ######
        self.reader_type = StringVar()
        self.reader_combobox = ttk.Combobox(self.reader_frame, textvariable=self.reader_type, values=READERS)
        self.reader_combobox.bind("<<ComboboxSelected>>", self.reader_selected)
        self.csvReader_flag = False
        self.csvReader_options = False

        self.reader_frame.grid(column=0, row=0, sticky=(N, S, E, W))
        self.reader_combobox.grid(column=0, row=0, sticky=(E, W))
        self.reader_frame.columnconfigure(0, weight=1)
        self.reader_frame.rowconfigure(0, weight=1)

    def reader_selected(self, event: Event) -> None:
        combo: ttk.Combobox = event.widget
        value = combo.get()
        self.csvReader_flag = True if value == "csvReader" else False

        # Add csvReader options
        if self.csvReader_flag is True and self.csvReader_options is False:
            self.csvReader_options = True
            self.csvReader_frame = ttk.Frame(self.reader_frame, padding="0 20 0 0")
            self.csvReader_entry_frame = ttk.Frame(self.csvReader_frame, padding="0 0 10 15")
            self.csvReader_entry_label = ttk.Label(self.csvReader_entry_frame, text="Filename:")
            ###### Important Var ######
            self.csvReader_filename = StringVar()
            self.csvReader_entry = ttk.Entry(self.csvReader_entry_frame, textvariable=self.csvReader_filename)
            self.csvReader_browse_button = ttk.Button(self.csvReader_frame, text="Browse", command=self.csvReader_openfile_callback)

            # Display
            self.csvReader_frame.grid(column=0, row=1, sticky=(N, S, E, W))
            self.csvReader_entry_frame.grid(column=0, row=0, sticky=(N, S, E, W))
            self.csvReader_entry_label.grid(column=0, row=0, sticky=W)
            self.csvReader_entry.grid(column=0, row=1)
            self.csvReader_browse_button.grid(column=1, row=0)
        
        # Destroy all
        if self.csvReader_flag is False and self.csvReader_options is True:
            self.csvReader_frame.destroy()
            self.csvReader_options = False
            
    def csvReader_openfile_callback(self) -> None:
        filename = filedialog.askopenfilename()
        self.csvReader_filename.set(filename)
    #################################################

    ######### WRITER LOGIC ###########
    def init_writer(self, frame: ttk.Frame) -> None:
        self.writer_frame = ttk.Labelframe(frame, padding=PADDING/2, text="Writer")
        ###### Important Var ######
        self.writer_type = StringVar()
        self.writer_combobox = ttk.Combobox(self.writer_frame, textvariable=self.writer_type, values=WRITERS)
        self.writer_combobox.bind("<<ComboboxSelected>>", self.writer_selected)
        self.csvWriter_flag = False
        self.csvWriter_options = False

        self.writer_frame.grid(column=0, row=1, sticky=(N, S, E, W))
        self.writer_combobox.grid(column=0, row=0, sticky=(E, W))
        self.writer_frame.columnconfigure(0, weight=1)

    def writer_selected(self, event: Event) -> None:
        combo: ttk.Combobox = event.widget
        value = combo.get()
        self.csvWriter_flag = True if value == "csvWriter" else False

        # Add csvReader options
        if self.csvWriter_flag is True and self.csvWriter_options is False:
            self.csvWriter_options = True
            self.csvWriter_frame = ttk.Frame(self.writer_frame, padding="0 20 0 0")
            self.csvWriter_entry_label = ttk.Label(self.csvWriter_frame, text="Result filename:")
            ###### Important Var ######
            self.csvWriter_filename = StringVar(value="result.csv")
            self.csvWriter_entry = ttk.Entry(self.csvWriter_frame, textvariable=self.csvWriter_filename)

            # Display
            self.csvWriter_frame.grid(column=0, row=1, sticky=(N, S, E, W))
            self.csvWriter_entry_label.grid(column=0, row=0, sticky=W)
            self.csvWriter_entry.grid(column=0, row=1)
        
        # Destroy all
        if self.csvWriter_flag is False and self.csvWriter_options is True:
            self.csvWriter_frame.destroy()
            self.csvWriter_options = False
    
    #################################################

    ######### LOGGER LOGIC ###########

    def init_logger(self, frame: ttk.Frame) -> None:
        self.logger_frame = ttk.Labelframe(frame, padding=PADDING/2, text="Logger")
        ###### Important Var ######
        self.logger_type = StringVar()
        self.logger_combobox = ttk.Combobox(self.logger_frame, textvariable=self.logger_type, values=LOGGERS)

        self.logger_level_frame = ttk.Frame(self.logger_frame, padding="0 10 0 0")
        self.logger_level_label = ttk.Label(self.logger_level_frame, text="Log level:")
        self.logger_level = StringVar()
        self.logger_level_combobox = ttk.Combobox(self.logger_level_frame, textvariable=self.logger_level, values=LOGGER_LEVELS)

        self.logger_frame.grid(column=0, row=2, sticky=(N, S, E, W))
        self.logger_combobox.grid(column=0, row=0, sticky=(E, W))
        self.logger_level_frame.grid(column=0, row=1, sticky=(N, S, E, W))
        self.logger_level_label.grid(column=0, row=0, sticky=W)
        self.logger_level_combobox.grid(column=0, row=1, sticky=(E, W))

        self.logger_frame.columnconfigure(0, weight=1)
        self.logger_level_frame.columnconfigure(0, weight=1)

    ######### OUTFOLDER LOGIC ###########
    def init_outfolder(self, frame: ttk.Frame) -> None:
        self.outfolder_frame = ttk.Labelframe(frame, padding=PADDING, text="Out Execution Directory")
        info_message = "Default value \"null\" save the execution results in the \"./out\" directory (relative to the compiled binary)"
        self.outfolder_info_label = ttk.Label(self.outfolder_frame, text=info_message, padding="0 20 0 0")
        self.outfolder_browse_button = ttk.Button(self.outfolder_frame, text="Browse", command=self.outfolder_openfolder_callback)

        self.outfolder_entry_frame = ttk.Frame(self.outfolder_frame, padding="0 5 0 0")
        self.outfolder_value = StringVar(value="null")
        self.outfolder_entry = ttk.Entry(self.outfolder_entry_frame, textvariable=self.outfolder_value, width=70)

        self.outfolder_frame.grid(column=0, columnspan=2, row=1, sticky=(N, S, E, W))
        self.outfolder_entry_frame.grid(column=0, row=0, sticky=(N, S, E, W))
        self.outfolder_info_label.grid(column=0, row=1, sticky=W)
        self.outfolder_browse_button.grid(column=1, row=0)
        self.outfolder_entry.grid(column=0, row=0, sticky=(E, W))
    
    def outfolder_openfolder_callback(self) -> None:
        filename = filedialog.askdirectory()
        self.outfolder_value.set(filename)

    #################################################

    ######### INTEGRATOR LOGIC ###########
    def init_integrator(self, frame: ttk.Frame) -> None:
        self.integrator_frame = ttk.Labelframe(frame, padding=PADDING, text="Integrator")
        self.integrator_type = StringVar()
        self.integrator_combobox = ttk.Combobox(self.integrator_frame, textvariable=self.integrator_type, values=INTEGRATORS)
        self.integrator_combobox.bind("<<ComboboxSelected>>", self.integrator_combobox_selected)

        self.integrator_mechanism_frame = ttk.Frame(self.integrator_frame, padding="0 10 0 0")
        self.integrator_mechanism_value = StringVar()
        self.integrator_mechanism_label = ttk.Label(self.integrator_mechanism_frame, text="Mechanism:")
        self.integrator_mechanism_combobox = ttk.Combobox(self.integrator_mechanism_frame, textvariable=self.integrator_mechanism_value, values=get_mechanism_list())

        self.integrator_reltol_frame = ttk.Frame(self.integrator_frame, padding="0 10 0 0")
        self.integrator_reltol_label = ttk.Label(self.integrator_reltol_frame, text="Relative tolerance:")
        self.integrator_reltol_value = DoubleVar(value=1.0e-6)
        self.integrator_reltol_entry = ttk.Entry(self.integrator_reltol_frame, textvariable=self.integrator_reltol_value)

        self.integrator_abstol_frame = ttk.Frame(self.integrator_frame, padding="0 10 0 0")
        self.integrator_abstol_label = ttk.Label(self.integrator_abstol_frame, text="Absolute tolerance:")
        self.integrator_abstol_value = DoubleVar(value=1.0e-10)
        self.integrator_abstol_entry = ttk.Entry(self.integrator_abstol_frame, textvariable=self.integrator_abstol_value)

        self.integrator_pressure_frame = ttk.Frame(self.integrator_frame, padding="0 10 0 0")
        self.integrator_pressure_label = ttk.Label(self.integrator_pressure_frame, text="Pressure:")
        self.integrator_pressure_value = DoubleVar(value=101325.15)
        self.integrator_pressure_entry = ttk.Entry(self.integrator_pressure_frame, textvariable=self.integrator_pressure_value)

        self.integrator_time_frame = ttk.Frame(self.integrator_frame, padding="0 10 0 0")
        self.integrator_time_label = ttk.Label(self.integrator_time_frame, text="Integration time:")
        self.integrator_time_value = DoubleVar(value=1.0e-3)
        self.integrator_time_entry = ttk.Entry(self.integrator_time_frame, textvariable=self.integrator_time_value)

        # OpenMP
        self.openmp_frame = ttk.Labelframe(self.integrator_frame, padding="0 10 0 0", text="OpenMP")
        
        self.openmp_cpus_frame = ttk.Frame(self.openmp_frame, padding="0 10 0 0")
        self.openmp_cpus_label = ttk.Label(self.openmp_cpus_frame, text="Number of CPUs:")
        self.openmp_cpus_value = IntVar()
        self.openmp_cpus_entry = ttk.Entry(self.openmp_cpus_frame, textvariable=self.openmp_cpus_value)

        self.openmp_schedule_frame = ttk.Frame(self.openmp_frame, padding="0 10 0 0")
        self.openmp_schedule_label = ttk.Label(self.openmp_schedule_frame, text="For loop Schedule:")
        self.openmp_schedule_value = StringVar()
        self.openmp_schedule_entry = ttk.Combobox(self.openmp_schedule_frame, textvariable=self.openmp_schedule_value, values=OMP_SCHEDULES)

        self.openmp_chunk_frame = ttk.Frame(self.openmp_frame, padding="0 10 0 0")
        self.openmp_chunk_label = ttk.Label(self.openmp_chunk_frame, text="For loop Chunk:")
        self.openmp_chunk_value = IntVar()
        self.openmp_chunk_entry = ttk.Entry(self.openmp_chunk_frame, textvariable=self.openmp_chunk_value)

        # Display
        self.integrator_frame.grid(column=0, row=0, sticky=(N, S, E, W))
        self.integrator_frame.columnconfigure(0, weight=1)
        self.integrator_combobox.grid(column=0, row=0, sticky=(E, W))
        self.integrator_mechanism_frame.grid(column=0, row=1, sticky=(N, S, E, W))

        self.integrator_mechanism_label.grid(column=0, row=0, sticky=W)
        self.integrator_mechanism_combobox.grid(column=0, row=1, sticky=(E, W))

        self.integrator_reltol_frame.grid(column=0, row=2, sticky=(N, S, E, W))
        self.integrator_reltol_label.grid(column=0, row=0, sticky=W)
        self.integrator_reltol_entry.grid(column=0, row=1, sticky=(E, W))

        self.integrator_abstol_frame.grid(column=0, row=3, sticky=(N, S, E, W))
        self.integrator_abstol_label.grid(column=0, row=0, sticky=W)
        self.integrator_abstol_entry.grid(column=0, row=1, sticky=(E, W))





    def integrator_combobox_selected(self, event: Event) -> None:
        print(event.widget)

if __name__ == "__main__":
    window = Tk()
    configPage = ConfigYamlPage(window)
    window.mainloop()