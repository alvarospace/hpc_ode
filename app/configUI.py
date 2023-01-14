from pathlib import Path
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
LABELFRAME_PADDING = 10
INPUT_PADDING = 2.5

class InputBox:
    def __init__(self, parent: ttk.Frame, label: str, variable: Variable) -> None:
        self.frame = ttk.Frame(parent, padding=INPUT_PADDING)
        self.label = ttk.Label(self.frame, text=label)
        self.variable = variable

    def draw(self, column: int, row: int) -> None:
        self.frame.grid(column=column, row=row, sticky=(N, S, E, W))
        self.frame.columnconfigure(0, weight=1)
        # Draw inside frame
        self.label.grid(column=0, row=0, sticky=W)

class ComboBox(InputBox):
    def __init__(self, parent: ttk.Frame, label: str, variable: Variable, values: list[str]) -> None:
        super().__init__(parent, label, variable)
        self.combobox = ttk.Combobox(self.frame, textvariable=self.variable, values=values)
    
    def draw(self, column: int, row: int) -> None:
        super().draw(column, row)
        self.combobox.grid(column=0, row=1, sticky=(W, E))

    def add_callback(self, callback: Callable[[Event], None]) -> None:
        self.combobox.bind("<<ComboboxSelected>>", callback)

class EntryBox(InputBox):
    def __init__(self, parent: ttk.Frame, label: str, variable: Variable) -> None:
        super().__init__(parent, label, variable)
        self.entry = ttk.Entry(self.frame, textvariable=self.variable)

    def draw(self, column: int, row: int) -> None:
        super().draw(column, row)
        self.entry.grid(column=0, row=1, sticky=(W, E))

class BrowseBox(InputBox):
    def __init__(self, parent: ttk.Frame, label: str, variable: Variable, dir: bool=False) -> None:
        super().__init__(parent, label, variable)
        self.sub_frame = ttk.Frame(self.frame)
        self.entry = ttk.Entry(self.sub_frame, textvariable=self.variable)
        self.button_frame = ttk.Frame(self.sub_frame, padding="5 0 0 0")
        command: Callable[[], None] = self.browse_dir if dir is True else self.browse_file
        self.button = ttk.Button(self.button_frame, text="Browse", command=command)

    def draw(self, column: int, row: int) -> None:
        super().draw(column, row)
        self.sub_frame.grid(column=0, row=1, sticky=(N, S, E, W))
        self.sub_frame.columnconfigure(0, weight=1)

        self.entry.grid(column=0, row=0, sticky=(W, E))
        self.button_frame.grid(column=1, row=0)
        self.button.grid(column=0, row=0)

    def browse_file(self) -> None:
        filename = filedialog.askopenfilename()
        self.variable.set(filename)

    def browse_dir(self) -> None:
        dirname = filedialog.askdirectory()
        self.variable.set(dirname)

class ConfigYamlPage:
    def __init__(self, root: Tk):
        self.root = root
        self.root.title("YAML Configuration")
        self.root.attributes('-fullscreen', True)
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        # Variables
        self.yaml_dict = {}
        self.reader_type = StringVar()
        self.csvreader_filename = StringVar()
        self.writer_type = StringVar()
        self.csvwriter_filename = StringVar(value="result.csv")
        self.logger_type = StringVar()
        self.logger_level = StringVar()
        self.outfolder = StringVar(value=None)
        self.integrator_type = StringVar()
        self.integrator_mechanism = StringVar()
        self.integrator_reltol = DoubleVar(value=1.0e-6)
        self.integrator_abstol = DoubleVar(value=1.0e-10)
        self.integrator_pressure = DoubleVar(value=101325.15)
        self.integrator_time = DoubleVar(value=1.0e-3)
        self.openmp_cpus = IntVar()
        self.openmp_schedule = StringVar()
        self.openmp_chunk = IntVar()
        
        self.message = StringVar()

        # Grid structure
        self.main_frame = ttk.Frame(self.root, padding=PADDING, borderwidth=3, relief="raised")
        self.config_frame = ttk.Frame(self.main_frame, borderwidth=2, relief="raised")
        self.left_column = ttk.Frame(self.config_frame, padding=f"{PADDING} {PADDING} {PADDING} 0")
        self.right_column = ttk.Frame(self.config_frame, padding=f"0 {PADDING} {PADDING} {PADDING}")
        self.yaml_frame = ttk.Frame(self.main_frame, padding=PADDING, borderwidth=2, relief="raised")
        self.bottom = ttk.Frame(self.main_frame, padding=PADDING, borderwidth=2, relief="raised")
        
        self.main_frame.grid(column=0, row=0, sticky=(N, S, E, W))
        self.main_frame.columnconfigure(0, weight=1)
        self.main_frame.columnconfigure(1, weight=1)
        self.main_frame.rowconfigure(0, weight=1)

        # TODO: TEXT widget to the right to show the yaml
        # TODO: Message and button widgets

        self.config_frame.grid(column=0, row=0, sticky=(N, S, E, W))
        self.config_frame.columnconfigure(0, weight=1)
        self.config_frame.columnconfigure(1, weight=1)
        self.config_frame.rowconfigure(0, weight=1)

        self.left_column.grid(column=0, row=0, sticky=(N, S, E, W))
        self.left_column.columnconfigure(0, weight=1)

        self.right_column.grid(column=1, row=0, sticky=(N, S, E, W))
        self.right_column.columnconfigure(0, weight=1)

        self.yaml_frame.grid(column=1, row=0, rowspan=2, sticky=(N, S, E, W))
        self.yaml_frame.columnconfigure(0, weight=1)
        self.yaml_frame.rowconfigure(0, weight=1)
        
        self.bottom.grid(column=0, row=1, sticky=(N, S, E, W))
        self.bottom.columnconfigure(0, weight=1)
        self.bottom.rowconfigure(0, weight=1)

        # Create components
        self.init_reader(self.left_column)
        self.init_writer(self.left_column)
        self.init_logger(self.left_column)
        self.init_outfolder(self.left_column)
        self.init_integrator(self.right_column)
        self.init_yaml(self.yaml_frame)
        self.init_bottom(self.bottom)

    ######### YAML LOGIC ########
    def init_yaml(self, frame: ttk.Frame) -> None:
        self.yaml_text = Text(frame, wrap="none")
        self.vert_scroll = ttk.Scrollbar(frame, orient=VERTICAL, command=self.yaml_text.yview)
        self.horiz_scroll = ttk.Scrollbar(frame, orient=HORIZONTAL, command=self.yaml_text.xview)
        self.yaml_text.configure(yscrollcommand=self.vert_scroll.set, xscrollcommand=self.horiz_scroll.set)

        self.yaml_text.grid(column=0, row=0, sticky=(N, S, E, W))
        self.vert_scroll.grid(column=1, row=0, sticky=(N, S))
        self.horiz_scroll.grid(column=0, row=1, sticky=(W, E))

    def is_input_valid(self) -> bool:
        is_valid = True
        variables_minimum_list: list[Variable] = [
            self.reader_type,
            self.writer_type,
            self.logger_type,
            self.logger_level,
            self.outfolder,
            self.integrator_type,
            self.integrator_mechanism,
            self.integrator_reltol,
            self.integrator_abstol,
            self.integrator_pressure,
            self.integrator_time,
        ]
        openmp_list: list[Variable] = [
            self.openmp_cpus,
            self.openmp_schedule,
            self.openmp_chunk
        ]
        for min_var in variables_minimum_list:
            if not min_var.get():
                is_valid = False

        if not self.is_dir(self.outfolder.get()):
            is_valid = False

        if self.reader_type.get() == "csvReader":
            if not self.csvreader_filename.get():
                is_valid = False
            else:
                is_valid = self.is_csv(self.csvreader_filename.get())

        if self.writer_type.get() == "csvWriter":
            if not self.csvwriter_filename.get():
                is_valid = False
            else:
                is_valid = self.is_csv(self.csvwriter_filename.get())
        
        if self.integrator_type.get().endswith("OMP"):
            for omp_var in openmp_list:
                if not omp_var.get():
                    is_valid = False

        return is_valid

    def is_csv(self, filename: str) -> bool:
        return True if filename.endswith(".csv") else False

    def is_dir(self, dirname: str) -> bool:
        return Path(dirname).is_dir()

    ######## GENERATE LOGIC ########
    def init_bottom(self, frame: ttk.Frame) -> None:
        self.message_label = ttk.Label(frame, textvariable=self.message, padding=PADDING, anchor=CENTER)
        self.buttons_frame = ttk.Frame(frame, padding=LABELFRAME_PADDING)
        self.generate_button = ttk.Button(self.buttons_frame, text="Generate", command=self.generate_yaml)
        self.saveas_button = ttk.Button(self.buttons_frame, text="Save as", state="disabled")

        self.message_label.grid(column=0, row=0, sticky=(N, S, E, W))
        self.buttons_frame.grid(column=0, row=1, sticky=(N, S, E, W))
        self.buttons_frame.columnconfigure(0, weight=1)
        self.buttons_frame.columnconfigure(1, weight=1)
        self.buttons_frame.rowconfigure(0, weight=1)
        self.generate_button.grid(column=0, row=0, sticky=(N, S, E, W))
        self.saveas_button.grid(column=1, row=0, sticky=(N, S, E, W))

    def generate_yaml(self):
        self.message.set("")
        if self.is_input_valid():
            # Init yaml dict fields
            for field in ["outFileService", "reader", "writer", "integrator", "logger"]:
                self.yaml_dict[field] = {}

            # Out file service
            self.yaml_dict["outFileService"]["outFolder"] = self.outfolder.get()

            # Reader
            self.yaml_dict["reader"]["type"] = self.reader_type.get()
            if self.reader_type.get() == "csvReader":
                self.yaml_dict["reader"]["filename"] = self.csvreader_filename.get()
            
            # Writer
            self.yaml_dict["writer"]["type"] = self.writer_type.get()
            if self.writer_type.get() == "csvWriter":
                self.yaml_dict["writer"]["filename"] = self.csvwriter_filename.get()

            # Logger
            self.yaml_dict["logger"]["type"] = self.logger_type.get()
            self.yaml_dict["logger"]["logLevel"] = self.logger_level.get()

            # Integrator
            self.yaml_dict["integrator"]["type"] = self.integrator_type.get()
            self.yaml_dict["integrator"]["mechanism"] = self.integrator_mechanism.get()
            self.yaml_dict["integrator"]["reltol"] = self.integrator_reltol.get()
            self.yaml_dict["integrator"]["abstol"] = self.integrator_abstol.get()
            self.yaml_dict["integrator"]["pressure"] = self.integrator_pressure.get()
            self.yaml_dict["integrator"]["dt"] = self.integrator_time.get()
            if self.integrator_type.get().endswith("OMP"):
                self.yaml_dict["integrator"]["omp"] = {}
                self.yaml_dict["integrator"]["omp"]["schedule"] = {}
                self.yaml_dict["integrator"]["omp"]["cpus"] = self.openmp_cpus.get()
                self.yaml_dict["integrator"]["omp"]["schedule"]["type"] = self.openmp_schedule.get()
                self.yaml_dict["integrator"]["omp"]["schedule"]["chunk"] = self.openmp_chunk.get()

            self.yaml_text.delete("1.0", "end")
            self.yaml_text.insert("1.0", yaml.dump(self.yaml_dict))
        else:
            self.message.set("Some missing values, please fill in all the fields")
        

    ######### READER LOGIC ###########
    def init_reader(self, frame: ttk.Frame) -> None:
        self.reader_frame = ttk.Labelframe(frame, padding=LABELFRAME_PADDING, text="Reader")
        self.reader_combo = ComboBox(self.reader_frame, "Type:", self.reader_type, READERS)
        self.reader_combo.add_callback(self.reader_selected)
        self.reader_combo.draw(0,0)

        self.csvReader_flag = False
        self.csvReader_options = False

        self.reader_frame.grid(column=0, row=0, sticky=(N, S, E, W))
        self.reader_frame.columnconfigure(0, weight=1)

    def reader_selected(self, event: Event) -> None:
        combo: ttk.Combobox = event.widget
        value = combo.get()
        self.csvReader_flag = True if value == "csvReader" else False

        # Add csvReader options
        if self.csvReader_flag is True and self.csvReader_options is False:
            self.csvReader_options = True
            self.csvReader_browsebox = BrowseBox(self.reader_frame, "Filename:", self.csvreader_filename)
            self.csvReader_browsebox.draw(0, 1)
        
        # Destroy all
        if self.csvReader_flag is False and self.csvReader_options is True:
            self.csvReader_browsebox.frame.destroy()
            self.csvReader_options = False

    ######### WRITER LOGIC ###########
    def init_writer(self, frame: ttk.Frame) -> None:
        self.writer_frame = ttk.Labelframe(frame, padding=LABELFRAME_PADDING, text="Writer")
        self.writer_combo = ComboBox(self.writer_frame, "Type:", self.writer_type, WRITERS)
        self.writer_combo.add_callback(self.writer_selected)
        self.writer_combo.draw(0, 0)

        self.csvWriter_flag = False
        self.csvWriter_options = False

        self.writer_frame.grid(column=0, row=1, sticky=(N, S, E, W))
        self.writer_frame.columnconfigure(0, weight=1)

    def writer_selected(self, event: Event) -> None:
        combo: ttk.Combobox = event.widget
        value = combo.get()
        self.csvWriter_flag = True if value == "csvWriter" else False

        # Add csvReader options
        if self.csvWriter_flag is True and self.csvWriter_options is False:
            self.csvWriter_options = True
            self.csvWriter_entry = EntryBox(self.writer_frame, "Filename:", self.csvwriter_filename)
            self.csvWriter_entry.draw(0, 1)
        
        # Destroy all
        if self.csvWriter_flag is False and self.csvWriter_options is True:
            self.csvWriter_entry.frame.destroy()
            self.csvWriter_options = False

    ######### LOGGER LOGIC ###########
    def init_logger(self, frame: ttk.Frame) -> None:
        self.logger_frame = ttk.Labelframe(frame, padding=LABELFRAME_PADDING, text="Logger")
        self.logger_type_combo = ComboBox(self.logger_frame, "Type:", self.logger_type, LOGGERS)
        self.logger_type_combo.draw(0, 0)
        self.logger_level_combo = ComboBox(self.logger_frame, "Log level:", self.logger_level, LOGGER_LEVELS)
        self.logger_level_combo.draw(0, 1)

        self.logger_frame.grid(column=0, row=2, sticky=(N, S, E, W))
        self.logger_frame.columnconfigure(0, weight=1)

    ######### OUTFOLDER LOGIC ###########
    def init_outfolder(self, frame: ttk.Frame) -> None:
        self.outfolder_frame = ttk.Labelframe(frame, padding=LABELFRAME_PADDING, text="Evidences directory")
        self.outfolder_browsebox = BrowseBox(self.outfolder_frame, "Folder:", self.outfolder, True)
        self.outfolder_browsebox.draw(0, 0)
        # info_message = "Default value \"null\" save the execution results in the \"./out\" directory (relative to the compiled binary)"
        # self.outfolder_info_label = ttk.Label(self.outfolder_frame, text=info_message, wraplength="10cm")

        self.outfolder_frame.grid(column=0, row=3, sticky=(N, S, E, W))
        self.outfolder_frame.columnconfigure(0, weight=1)

    ######### INTEGRATOR LOGIC ###########
    def init_integrator(self, frame: ttk.Frame) -> None:
        self.integrator_frame = ttk.Labelframe(frame, padding=LABELFRAME_PADDING, text="Integrator")

        self.integrator_combo = ComboBox(self.integrator_frame, "Type:", self.integrator_type, INTEGRATORS)
        self.integrator_combo.add_callback(self.integrator_selected)
        self.integrator_combo.draw(0, 0)

        self.integrator_mechanism_entry = ComboBox(self.integrator_frame, "Mechanism:", self.integrator_mechanism, get_mechanism_list())
        self.integrator_mechanism_entry.draw(0, 1)

        self.integrator_reltol_entry = EntryBox(self.integrator_frame, "Relative tolerance:", self.integrator_reltol)
        self.integrator_reltol_entry.draw(0, 2)

        self.integrator_abstol_entry = EntryBox(self.integrator_frame, "Absolute tolerance:", self.integrator_abstol)
        self.integrator_abstol_entry.draw(0, 3)

        self.integrator_pressure_entry = EntryBox(self.integrator_frame, "Pressure:", self.integrator_pressure)
        self.integrator_pressure_entry.draw(0, 4)

        self.integrator_time_entry = EntryBox(self.integrator_frame, "Integration time (dt):", self.integrator_time)
        self.integrator_time_entry.draw(0, 5)

        self.openmp_flag = False
        self.openmp_options = False

        # Display
        self.integrator_frame.grid(column=0, row=0, sticky=(N, S, E, W))
        self.integrator_frame.columnconfigure(0, weight=1)

    def integrator_selected(self, event: Event) -> None:
        combo: ttk.Combobox = event.widget
        value = combo.get()
        self.openmp_flag = value.endswith("OMP")

        # Draw OpenMP options
        if self.openmp_flag is True and self.openmp_options is False:
            self.openmp_options = True
            self.openmp_frame = ttk.Labelframe(self.integrator_frame, padding=LABELFRAME_PADDING, text="OpenMP")
            self.openmp_cpus_frame = ttk.Frame(self.openmp_frame, padding="2.5 0 0 0")
            self.openmp_schedule_frame = ttk.Frame(self.openmp_frame, padding="0 0 2.5 0")

            self.openmp_cpus_entry = EntryBox(self.openmp_cpus_frame, "Threads (CPUs):", self.openmp_cpus)
            self.openmp_cpus_entry.draw(0, 0)
            self.openmp_schedule_combobox = ComboBox(self.openmp_schedule_frame, "Schedule:", self.openmp_schedule, OMP_SCHEDULES)
            self.openmp_schedule_combobox.draw(0, 0)
            self.openmp_chunk_entry = EntryBox(self.openmp_schedule_frame, "Chunk:", self.openmp_chunk)
            self.openmp_chunk_entry.draw(0, 1)

            self.openmp_frame.grid(column=0, row=6, sticky=(N, S, E, W))
            self.openmp_frame.columnconfigure(0, weight=1)
            self.openmp_frame.columnconfigure(1, weight=1)
            self.openmp_schedule_frame.grid(column=0, row=0, sticky=(N, S, E, W))
            self.openmp_schedule_frame.columnconfigure(0, weight=1)
            self.openmp_cpus_frame.grid(column=1, row=0, sticky=(N, S, E, W))
            self.openmp_cpus_frame.columnconfigure(0, weight=1)

        # Delete OpenMP options
        if self.openmp_flag is False and self.openmp_options is True:
            self.openmp_frame.destroy()
            self.openmp_options = False

if __name__ == "__main__":
    window = Tk()
    configPage = ConfigYamlPage(window)
    window.mainloop()