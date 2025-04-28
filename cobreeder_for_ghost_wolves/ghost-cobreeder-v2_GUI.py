import markdown
import numpy as np
import os
import pandas as pd
from pandastable import Table
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import filedialog, scrolledtext, messagebox
from tkinterweb import HtmlFrame


def access_csv(p):
    """Browse files to identify CSV paths."""
    path = filedialog.askopenfilename()
    p.set(path)


def preview_individuals():
    """Preview the CSV or Excel file as a DataFrame in a new window."""
    try:
        p = individuals_file.get()
        if ('.csv' not in p) and ('.xlsx' not in p) and ('.xls' not in p):
            raise Exception("Files must be in CSV or Excel format.")
        else:
            df = pd.read_csv(p) if '.csv' in p else pd.read_excel(p)
            df.reset_index(inplace=True)
            df.rename(columns={'index': 'Index'}, inplace=True)
            df['Index'] = df['Index'].astype('Int64')
            preview_window = tk.Toplevel()
            preview_window.title("File Preview")
            preview_frame = tk.Frame(preview_window)
            preview_frame.pack(expand=True, fill="both")
            Table(preview_frame, dataframe=df).show()

    except Exception as e:
        messagebox.showerror("ERROR", str(e))


def open_readme():
    """Open README file in a new window."""
    file_path = "../README-V2_GUI.md"
    try:
        with open(file_path, 'r') as f:
            content = f.read()
            readme_window = tk.Toplevel()
            readme_window.title("CoBreeder for Ghost Wolves (v2) Usage Guidelines")
            html_frame_container = tk.Frame(readme_window)
            html_frame_container.pack(fill="both", expand=True)
            html_frame = HtmlFrame(html_frame_container)
            html_frame.pack(fill="both", expand=True)
            html_frame.load_html(markdown.markdown(content, extensions=['tables']))
    except Exception as e:
        messagebox.showerror("ERROR", str(e))


def input_validation():
    """Performs input validation for fields."""

    try:
        if int(num_pairs.get()) < 0:
            raise Exception("\nNumber of Pairings cannot be negative.")
        if int(global_pr_threshold.get()) < 0:
            raise Exception("\nGlobal PR Threshold cannot be negative.")
        if int(prio_calc_threshold.get()) < 0:
            raise Exception("\nPriority Calculations: Size of Priority Set cannot be negative.")

        # Files must be CSVs with no missing values
        pr = pd.read_csv(pairwise_relatedness_file.get(), delimiter=',', header=None, skiprows=1)
        individuals = pd.read_csv(individuals_file.get(), delimiter=',')
        if np.any(pr.isnull()) or np.any(individuals.isnull()):
            raise Exception("\nCSV files contain missing values.")

        # PR file must be symmetrical and only contain non-negative integers
        if not (pr >= 0).values.all():
            raise Exception("\nPR matrix contains negative values.")
        if not all(np.issubdtype(d, np.integer) for d in pr.dtypes):
            raise Exception("\nPR matrix contains non-integer values.")
        if not pr.equals(pr.T):
            raise Exception("\nPR matrix is not symmetrical.")

        # Individuals specification file must contain the required columns
        if list(individuals.columns) != ['Name', 'Male', 'Female', 'AssignToGroup', 'Alleles', 'Proven', 'Priority']:
            raise Exception("\nUnexpected column in individuals file.")
        # Alleles and AssignToGroup must be integers
        if not all(np.issubdtype(d, np.integer) for d in individuals[['AssignToGroup', 'Alleles']].dtypes):
            raise Exception("\nNon-integer value identified in individuals specification file.")
        # Alleles cannot be negative
        if not all(individuals['Alleles'] >= 0):
            raise Exception("\nIndividual specification file contains negative alleles.")
        # AssignToGroup can only be -1 or a group ID
        if not ((individuals['AssignToGroup'] >= -1) & (individuals['AssignToGroup'] < int(num_pairs.get()))).all():
            raise Exception("\nInvalid value in AssignToGroup column.")
        # Proven, Priority, Male, and Female can only take values 0 or 1
        for col in ['Proven', 'Priority', 'Male', 'Female']:
            if not set(individuals[col]) <= {0, 1}:
                raise Exception("\nProven, Priority, Male, and Female can only be 0 or 1.")
        # Individuals can only be male or female
        if not ((individuals['Male'] + individuals['Female']) == 1).all():
            raise Exception("\nIndividuals can only be male or female.")

        if int(prio_calc_threshold.get()) > int(num_pairs.get()):
            raise Exception("The size of the priority set cannot be larger than the number of pairs.")

        if ("PRIO" in str(obj_function.get())) and int(prio_calc_threshold.get()) == 0:
            raise Exception("MIN_PR_MAX_ALLELES_MAX_PRIO and MAX_PRIO require priority calculations to be enabled.")

    except Exception as e:
        messagebox.showerror("Error", str(e))


def threading_func():
    """Run CoBreeder for Ghost Wolves main logic."""
    def cobreeder_for_ghost_wolves():
        try:
            input_validation()
            weight_alleles, weight_pr, weight_prio = weights.get().split(",")

            arguments = [
                sys.executable, "ghost-cobreeder-v2_for_GUI.py", "run",
                individuals_file.get(),
                pairwise_relatedness_file.get(),
                num_pairs.get(),
                specify_pr.get(),
                obj_function.get(),
                unique_run_id.get(),
                weight_alleles,
                weight_pr,
                weight_prio,
                global_pr_threshold.get(),
                exclude.get(),
                disallow.get(),
                prio_calc_threshold.get(),
                prio_calc_ghost_weight.get(),
                save_csv.get()
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
    "weights": "Weight (Alleles,PR,Priority)",
    "unique_run_id": "Unique Run ID",
    "global_pr_threshold": "Global PR Threshold",
    "specify_pr": "Custom PR Thresholds",
    "exclude": "Exclusions",
    "disallow": "Disallowed Pairings",
    "prio_calc_threshold": "Priority Calculations: Size of Priority Set",
    "prio_calc_ghost_weight": "Priority Calculations: Weight for Alleles",
    "save_csv": "Save final solution to CSV?"
}

# GUI colour palette
SIDEBAR_COLOUR = "#44403E"
TERMINAL_COLOUR = "#0B0707"
BUTTON_COLOUR = "#410C11"
TEXT_COLOUR = "#FFFFFF"

# Initialise GUI
root = tk.Tk()
root.configure(bg=SIDEBAR_COLOUR)
root.title("CoBreeder for Ghost Wolves (v2)")
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
tk.Button(sidebar, text="USAGE GUIDELINES", fg=TEXT_COLOUR, bg=BUTTON_COLOUR, command=open_readme).pack(pady=5)

# Parameters and default values
individuals_file = tk.StringVar()
pairwise_relatedness_file = tk.StringVar()
num_pairs = tk.StringVar()
obj_function = tk.StringVar(value="MIN_PR_MAX_ALLELES")
weights = tk.StringVar(value="1,1,0")
unique_run_id = tk.StringVar(value="ghost-wolf-experiment-1")
global_pr_threshold = tk.StringVar(value="0")
specify_pr = tk.StringVar()
exclude = tk.StringVar()
disallow = tk.StringVar()
prio_calc_threshold = tk.StringVar(value="0")
prio_calc_ghost_weight = tk.StringVar(value="0.5")
save_csv = tk.StringVar(value="YES")

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
    elif param in ["obj_function", "save_csv"]:
        dropdown = tk.OptionMenu
        if param == "obj_function":
            dropdown = tk.OptionMenu(sidebar, obj_function, *["MIN_AV_PR", "MAX_TOTAL_ALLELES",
                                                              "MAX_TOTAL_PRIO", "MIN_PR_MAX_ALLELES",
                                                              "MIN_PR_MAX_ALLELES_MAX_PRIO"]
                                     )
        if param == "save_csv":
            dropdown = tk.OptionMenu(sidebar, save_csv, *["YES", "NO"])
        dropdown.config(width=33, bg=TERMINAL_COLOUR, fg=TEXT_COLOUR, highlightbackground=SIDEBAR_COLOUR)
        dropdown.pack()

    # Text boxes for arguments requiring a number to be entered
    else:
        if param == "exclude":
            tk.Button(sidebar, text="Preview Individuals", command=preview_individuals, bg=BUTTON_COLOUR,
                      fg=TEXT_COLOUR).pack()
        tk.Entry(sidebar, textvariable=globals()[param], width=40, bg=TERMINAL_COLOUR, fg=TEXT_COLOUR,
                 insertbackground=TEXT_COLOUR).pack()

# Button that runs the solver
tk.Button(sidebar, text="Allocate", fg=TEXT_COLOUR, bg=BUTTON_COLOUR, command=threading_func).pack(pady=30)

root.mainloop()
