import os
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import filedialog, scrolledtext, messagebox

# TODO - Incorporate functionality from CLI e.g. save to csv, input validation exceptions etc.


def access_csv(p):
    path = filedialog.askopenfilename()
    p.set(path)


def threading_func():
    def cobreeder_for_ghost_wolves():
        try:
            arguments = [
                sys.executable, "ghost-cobreeder-v2.py", "run",
                individuals_file.get(),
                pairwise_relatedness_file.get(),
                num_pairs.get(),
                specify_pr.get(),
                obj_function.get(),
                unique_run_id.get(),
                weight_alleles.get(),
                weight_pr.get(),
                weight_prio.get(),
                global_pr_threshold.get(),
                exclude_disallow.get(),
                prio_calc_threshold.get()
            ]

            # Run ghost-cobreeder-v2.py as a subprocess
            ghost_cobreeder_process = subprocess.Popen(
                arguments,  # Inputs from GUI,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                cwd=os.path.dirname(os.path.abspath(__file__))
            )

            terminal.delete('1.0', tk.END)
            for line in iter(ghost_cobreeder_process.stdout.readline, ''):
                terminal.insert(tk.END, line)
                terminal.see(tk.END)
            ghost_cobreeder_process.stdout.close()
            ghost_cobreeder_process.wait()

        except Exception as e:
            messagebox.showerror("Error", str(e))

    t = threading.Thread(target=cobreeder_for_ghost_wolves)
    t.daemon = True
    t.start()


# Dictionary storing parameter labels
parameters = {
    "individuals_file": "Individuals Specification File",
    "pairwise_relatedness_file": "Scaled Pairwise Relatedness File",
    "num_pairs": "Number of Pairings",
    "obj_function": "Objective Function",
    "weight_pr": "Weight (PR)",
    "weight_alleles": "Weight (Alleles)",
    "weight_prio": "Weight (Priority)",
    "unique_run_id": "Unique Run ID",
    "global_pr_threshold": "Global PR Threshold",
    "specify_pr": "Custom PR Thresholds",
    "exclude_disallow": "Exclusions/Disallowed Pairings",
    "prio_calc_threshold": "Priority Calculations"
}

# GUI colour palette
SIDEBAR_COLOUR = "#44403E"
TERMINAL_COLOUR = "#0B0707"
BUTTON_COLOUR = "#410C11"
TEXT_COLOUR = "#FFFFFF"

# Initialise GUI
root = tk.Tk()
root.configure(bg=SIDEBAR_COLOUR)
root.title("CoBreeder for Ghost Wolves (v1)")
frame = tk.Frame(root, bg=SIDEBAR_COLOUR)
frame.pack(fill='both', expand=True)

# Add terminal display
terminal_frame = tk.Frame(frame, bg=SIDEBAR_COLOUR)
terminal_frame.pack(side='right', fill='both', expand=True)
terminal = scrolledtext.ScrolledText(terminal_frame, wrap=tk.WORD, height=40, padx=10, font=("Fixedsys", 15),
                                     bg=TERMINAL_COLOUR, fg=TEXT_COLOUR)
terminal.pack(fill='both', expand=True)

# Add parameter sidebar
sidebar = tk.Frame(frame, padx=10, pady=10, bg=SIDEBAR_COLOUR)
sidebar.pack(side='left', fill='y')
tk.Label(sidebar, text='PARAMETERS', fg=TEXT_COLOUR, bg=SIDEBAR_COLOUR).pack()

# Parameters and default values
individuals_file = tk.StringVar()
pairwise_relatedness_file = tk.StringVar()
num_pairs = tk.StringVar()
obj_function = tk.StringVar(value="Select Objective Function")
weight_pr = tk.StringVar(value="1")
weight_alleles = tk.StringVar(value="1")
weight_prio = tk.StringVar(value="1")
unique_run_id = tk.StringVar(value="experiment-name")
global_pr_threshold = tk.StringVar(value="0")
specify_pr = tk.StringVar(value="DEFAULT_PR")
exclude_disallow = tk.StringVar(value="ALL")
prio_calc_threshold = tk.StringVar(value="0")

for param in parameters:
    tk.Label(sidebar, text=parameters[param], fg=TEXT_COLOUR, bg=SIDEBAR_COLOUR).pack(anchor='w')

    # Browse local files to select CSV
    if param in ["individuals_file", "pairwise_relatedness_file"]:
        entry = tk.Entry(sidebar, textvariable=globals()[param], fg=TEXT_COLOUR, bg=TERMINAL_COLOUR,
                         insertbackground=TEXT_COLOUR, width=40)
        entry.pack()
        tk.Button(sidebar, text="Upload CSV", command=lambda p=globals()[param]: access_csv(p),
                  bg=BUTTON_COLOUR, fg=TEXT_COLOUR).pack()

    # Dropdown boxes for arguments with specific options
    elif param in ["obj_function", "specify_pr", "exclude_disallow"]:
        dropdown = tk.OptionMenu
        if param == "obj_function":
            dropdown = tk.OptionMenu(sidebar, obj_function, *["MIN_AV_PR", "MAX_TOTAL_ALLELES",
                                                              "MAX_TOTAL_PRIO", "MIN_PR_MAX_ALLELES",
                                                              "MIN_PR_MAX_ALLELES_MAX_PRIO"]
                                     )
        if param == "specify_pr":
            dropdown = tk.OptionMenu(sidebar, specify_pr, *["DEFAULT_PR", "CUSTOM_PR"])
        if param == "exclude_disallow":
            dropdown = tk.OptionMenu(sidebar, exclude_disallow, *["ALL", "EX"])
        dropdown.config(width=33, bg=TERMINAL_COLOUR, fg=TEXT_COLOUR, highlightbackground=SIDEBAR_COLOUR)
        dropdown.pack()

    # Text boxes for arguments requiring a number to be entered
    else:
        tk.Entry(sidebar, textvariable=globals()[param], width=40, bg=TERMINAL_COLOUR, fg=TEXT_COLOUR,
                 insertbackground=TEXT_COLOUR).pack()

# Button that runs the solver
tk.Button(sidebar, text="Allocate", fg=TEXT_COLOUR, bg=BUTTON_COLOUR, command=threading_func).pack(pady=30)

root.mainloop()
