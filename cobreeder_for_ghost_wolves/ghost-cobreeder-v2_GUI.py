import markdown
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
    file_path = "../README.md"
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


def threading_func():
    def cobreeder_for_ghost_wolves():
        try:
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
    "exclude": "Exclusions (leave blank to disable)",
    "disallow": "Disallowed Pairings (leave blank to disable)",
    "prio_calc_threshold": "Size of Priority Set",
    "prio_calc_ghost_weight": "Weight for Alleles (Priority Calculations)",
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
obj_function = tk.StringVar(value="Select Objective Function...")
weights = tk.StringVar(value="1,1,0")
unique_run_id = tk.StringVar(value="ghost-wolf-experiment-1")
global_pr_threshold = tk.StringVar(value="0")
specify_pr = tk.StringVar(value="e.g. 2-50,1-10")
exclude = tk.StringVar(value="e.g. 3,13,0")
disallow = tk.StringVar(value="e.g. 2-4,5-10")
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
