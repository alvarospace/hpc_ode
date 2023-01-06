from tkinter import ttk
from tkinter import *

from constants import (
    READERS,
    WRITERS,
    LOGGERS,
    LOGGER_LEVELS,
    INTEGRATORS,
    OMP_SCHEDULES,
    INPUT_DATA
)
from common import get_absolute_path_input_data


class ComboBox:
    def __init__(self, parent, title, options):
        self.Frame = ttk.Frame(parent, padding=10)
        self.Label = ttk.Label(self.Frame, text=title+":", anchor=W, borderwidth=2).grid(
            column=0,
            row=0,
            sticky=(W, E)
        )

        # Combobox
        current_var = StringVar()
        self.Combo = ttk.Combobox(self.Frame, textvariable=current_var, values=options).grid(
            column=0,
            row=1
        )

class ConfigYamlPage:
    COMBO_TYPES = {
        "Reader": READERS,
        "Writer": WRITERS,
        "Integrator": INTEGRATORS,
        "Logger": LOGGERS
    }
    
    def __init__(self, root):
        self.input_data = get_absolute_path_input_data()

        self.Content = ttk.Frame(root, padding=20).grid(
            column=0,
            row=0,
            sticky=(N, W, E, S)
        )

        self.set_header()

        self.combos: list[ComboBox] = []
        self.set_options()
        # Select Integrators
        # integrators = ComboBox()

    def set_header(self):
        self.Header = ttk.Frame(self.Content).grid(
            column=0,
            columnspan=2,
            row=0,
            sticky=(N, W, E, S)
        )
        self.Title = ttk.Label(self.Header, text="YAML Configuration", font="TkHeadingFont").grid(
            column=0,
            row=0
        )

    def set_options(self):
        self.Options = ttk.Frame(self.Content).grid(
            column=1,
            row=1,
            sticky=(N, W, E, S)
        )

        for title, types in self.COMBO_TYPES.items():
            combo = ComboBox(self.Options, title, types)
            self.combos.append(combo)
            row = len(self.combos)
            combo.Frame.grid(column=0, row=row)



        


if __name__ == "__main__":
    window = Tk()
    window.title("ODEIntegrator Application")
    configPage = ConfigYamlPage(window)
    window.mainloop()