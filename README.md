# CoBreeder for Ghost Wolves (v2) Usage Guidelines

---

These usage guidelines pertain to the **Command Line Interface** (CLI) offered by _CoBreeder for Ghost Wolves_. 
Please use the script `ghost_cobreeder_v2_CLI.py` when using the CLI.

_CoBreeder for Ghost Wolves_ also offers a **Graphical User Interface** (GUI). This will be more accessible to users
that are not familiar with the terminal. To use the GUI, please run the script `ghost_cobreeder_v2_GUI.py`. This
script will use an adjusted version of the CLI for its logic, `ghost_cobreeder_v2_for_GUI.py`, but you should not need
to interact with this file.
GUI usage guidelines can be accessed from within the GUI itself but are also found in `README-GUI.md`.

## <u>Context</u>

This project builds upon the work completed by _Forshaw et al._, titled "Constraint Optimisation Approaches for 
Designing Group-Living Captive Breeding Programmes" [1]. The tool produced (_CoBreeder_) was designed with a focus on 
the Galápagos tortoise which required single-objective optimisation to minimise pairwise relatedness when allocating 
individuals to breeding groups. The aim of this project is to adapt and extend _CoBreeder_ to improve its suitability 
for breeding coyotes with red wolf ancestry. This requires multi-objective optimisation, minimising genetic relatedness 
and maximising the number of ghost alleles in each pairing. Additional functionality that could aid conservationists 
has also been implemented, such as the ability to exclude individuals and specify disallowed pairings without editing 
the original data. There is also the option to calculate a priority score for each individual based on the 
characteristics of the dataset. 

This initial adaptation focuses on assigning coyote pairings however this may not be entirely accurate to life as 
coyotes, though sometimes solitary animals, often live in packs with a single breeding pair. Future versions of this 
tool should take this into account, improving its flexibility. Secondly, the pairwise relatedness file currently 
assumes that values are scaled, with larger values indicating that two individuals are less related. This should be 
developed further to allow users to specify more intuitive thresholds (e.g. third cousins). 


## <u>Command Line Arguments</u>

This program is run via the command line. You will need to specify various arguments and provide relative file paths to 
the data that will be used (please see the "Expected File Formats" section for more details). The program expects 
the following arguments, in the order introduced:

- `individuals_file`
  - Relative path to the CSV file detailing the individuals.
- `pairwise_relatedness_file`
  - Relative path to the CSV file containing the scaled pairwise relatedness matrix.
- `num_pairs`
  - The number of pairings/groups to allocate. 
- `specify_pr`
  - Specify custom PR thresholds for certain groups with "CUSTOM_PR", or use the default threshold for all groups with 
    "DEFAULT_PR".
- `obj_function`
  - String specifying the objective function to use when solving. The options are:
    - `MIN_AV_PR` - Minimise the average pairwise relatedness across pairings.
    - `MAX_TOTAL_ALLELES` - Maximise the sum of ghost alleles in a solution.
    - `MAX_TOTAL_PRIO` - Maximise the sum of priority scores in a solution.
    - `MIN_PR_MAX_ALLELES` - Minimise pairwise relatedness and maximise the number of ghost alleles in each pairing 
    using the weights specified. 
    - `MIN_PR_MAX_ALLELES_MAX_PRIO` - Minimise pairwise relatedness, maximise ghost alleles, and 
    maximise priority scores in each pairing, using the weights specified.
  - Note that `MIN_PR_MAX_ALLELES_MAX_PRIO` and `MAX_TOTAL_PRIO` cannot be used if `prio_calc_threshold` is set to 0 as
    they use dynamic priority calculations. 
- `unique_run_id`
  - Unique string identifier for each run. It has no bearing on functionality, but is useful to identify results when 
    running the program multiple times.
- `weight_alleles`
  - Integer specifying the relative weight to place on alleles when performing multi-objective optimisation.
- `weight_pr`
  - Integer specifying the relative weight to place on pairwise relatedness when performing multi-objective optimisation.
- `weight_prio`
  - Integer specifying the relative weight to place on priority values when performing multi-objective optimisation.
- `pr_threshold`
  - Threshold for the scaled pairwise relatedness permitted in a pairing (0 by default).
  - Value not included.
- `exclude_disallow`
  - Exclude individuals or specify disallowed pairings with "EX", or skip this step with "ALL".
  - This is explained in more detail in the "Excluded Individuals & Disallowed Pairings" section.
- `prio_calc_threshold`
  - Threshold for priority calculations representing the number of individuals that can fall into the priority set.
  - Must be a positive integer.
  - Specify 0 to disable dynamic priority calculations and use manual priority assignments only.
  - E.g., 4 will enable priority calculations and select the top 4 scoring individuals for the priority set.
  - This is explained in more detail in the "Priority Calculations" section.

### Example:

`poetry run python ghost_wolf_cobreeder/ghost_cobreeder_v2_CLI.py run data/individuals.csv data/pr-scaled.csv 5
CUSTOM_PR MIN_PR_MAX_ALLELES ghost_experiment 1 1 0 50 ALL 4`
  - 5 pairings will be allocated. 
  - The user will be prompted to enter any custom PR thresholds for specific groups. 
  - The solver will maximise the number of ghost alleles and minimise pairwise relatedness between individuals in each 
  pairing. Equal weight is placed on each of these.
  - The name of the run is "ghost_experiment". This string will be displayed alongside results.
  - Any individuals with a PR greater than 50 can be paired.
  - All individuals in the individuals file will be considered and all pairings are allowed, provided PR is appropriate.
  - Priority calculations are enabled, and the top 4 individuals with the best priority values must be included in 
    solutions.


## <u>Expected File Formats</u>

The program expects input files to be in CSV format with the column names specified below. The program expects two 
files, detailing the individuals and the pairwise relatedness between them. Please ensure that these files are complete. 
If the pairwise relatedness between two individuals is unknown, marking these cells with 0 will ensure that the 
individuals are not paired together. 

Note that data can be randomly generated for testing purposes using the script in 
`supporting_scripts/data-generation-v2.py`. This script will generate two CSV files (a set of individuals and a PR
matrix) and two histograms showing the distribution of alleles and PR values for this dataset. 

### Individuals Specification File:

- `Name`
  - String describing the individual. May include various details such as sex and location.
- `Male` & `Female` 
  - Boolean columns taking values 0 or 1, where 1 indicates that the individual's sex matches the name of the column.
  - For any row, both columns cannot take the same value.
- `AssignToGroup` 
  - Specifies whether to assign this individual to a specific group.
  - This column should be -1 by default. If specifying a group, please ensure that the value used matches a group ID.
- `Alleles` 
  - Positive integers representing the number of ghost alleles that each individual has.
- `Proven` 
  - Boolean column taking values 0 or 1, where 1 indicates that the individual is proven.
- `Priority` 
  - Boolean column taking values 0 or 1, where 1 indicates that the individual is a 'priority individual'.
  - This column is used to manually specify priority individuals. Note that if priority calculations are enabled, the 
    values in this column will not be considered. 

**Example**:

| Name           | Male | Female | AssignToGroup | Alleles | Proven | Priority |
|----------------|------|--------|---------------|---------|--------|----------|
| Individual_0_M | 1    | 0      | -1            | 295     | 1      | 0        |
| Individual_1_M | 1    | 0      | -1            | 312     | 0      | 0        |
| Individual_2_F | 0    | 1      | 1             | 488     | 1      | 1        |
| ...            | ...  | ...    | ...           | ...     | ...    | ...      |


### Scaled PR Specification File:

This file represents the scaled pairwise relatedness matrix and must include values for each individual in the dataset.
Note that the values in this file are scaled; 0 represents the same individual and is also used to specify
banned pairings. A larger value indicates that two individuals are less related; this may not seem intuitive,
but when minimising pairwise relatedness the program is actually maximising these values. 

It is important that every individual specified in the individual specification file is included here, and that the 
values are mirrored. Remember that, if editing the PR between two different individuals, you will need to edit two 
cells. Note, however, that you do not need to manually edit the file to exclude individuals or create disallowed 
pairings; this can be done whilst running the program.

**Example**:

| Individual_0_M | Individual_1_F | Individual_2_F | Individual_3_M | Individual_4_M |
|----------------|----------------|----------------|----------------|----------------|
| 0              | 217            | 749            | 474            | 538            |
| 217            | 0              | 581            | 798            | 616            |
| 749            | 581            | 0              | 207            | 542            |
| 474            | 798            | 207            | 0              | 368            |
| 538            | 616            | 542            | 368            | 0              |


## <u>Excluded Individuals & Disallowed Pairings</u>

1. If "EX" was specified on the command line, you will be shown a preview of the individuals file in table format. Each 
individual will have an ID on the far left-hand side, starting at 0.
2. You will first be prompted to specify disallowed pairings. 
   - To disallow pairings between individual 4 and 8, and 0 and 1 you will type: 4-8, 0-1
   - Order does not matter here. For example, you could also type: 1-0, 8-4 
   - To skip this, leave it blank and tap the enter key.
3. You will then be prompted to specify individuals to exclude.
   - To exclude individuals 2, 5, and 6 you will type: 2,5,6
   - Again, order does not matter here, and you can skip this by leaving it blank and tapping the enter key.

Note: neither of these options edit the original files. Your original data will be left untouched and settings will 
only apply to the current run.


## <u>Priority Calculations</u>

_CoBreeder_ allowed for the manual specification of 'priority individuals', individuals with 
desirable characteristics that must be included in solutions. This program retains this feature, but also offers 
dynamic priority calculations that allow individuals to be ranked by priority. Individuals are compared to their peers 
and assigned a priority score between 0 and 100, where 100 is the maximum priority value. The top x individuals are 
selected to fall into the priority set, representing individuals that must be included in solutions, whilst the 
priority scores are maximised in some objective functions.

Priority calculations consider factors such as whether individuals are proven, and the number of mates and ghost 
alleles that they have compared to their peers. If priority calculations are enabled, you will be prompted to provide 
a weight for ghost alleles. This is a float between 0 and 1, where 1 represents completely ignoring number of mates 
in favour of ghost alleles. 0.5 will strike a balance between these. It is recommended to prioritise ghost alleles more 
heavily at the start and increase emphasis on the number of mates if there are unforeseen circumstances (e.g. 
individuals are proving hard to capture). This will result in more flexible solutions.

## <u>Outputs</u>

Each feasible solution found will be printed to the terminal as the code is running. Once the search has finished, 
you will be prompted to choose whether to save the best solution to a CSV file. If this option is chosen, the solution 
with the best objective value will be saved to the `results` folder, with a file name following the format
`best_solution_<unique_run_id>.csv`. This will also be confirmed in the terminal. 

This file may look something like this:

| Group | Ind_1_Name      | Ind_2_Name      | Ind_1_ID | Ind_2_ID | Ind_1_Alleles | Ind_2_Alleles | Pairwise_Relatedness | Priority_Sum |
|-------|-----------------|-----------------|----------|----------|---------------|---------------|----------------------|--------------|
| 0     | Individual_38_F | Individual_41_M | 38       | 41       | 357           | 155           | 967                  | 129          |
| 1     | Individual_22_F | Individual_37_M | 22       | 37       | 781           | 281           | 877                  | 161          |
| 2     | Individual_7_F  | Individual_17_M | 7        | 17       | 261           | 482           | 967                  | 143          |
| 3     | Individual_27_F | Individual_42_M | 27       | 42       | 401           | 440           | 927                  | 73           |
| ...   | ...             | ...             | ...      | ...      | ...           | ...           | ...                  | ...          |

Each row represents a pairing, and contains information about the individuals such as name, ID, and number of ghost 
alleles, as well as their pairwise relatedness. The `Priority_Sum` column is the sum of the priority scores of the 
individuals, and can reach a maximum of 200. Note that, if priority calculations are not enabled, the 
cells in this column will be "N/A".

---

## <u>References</u>

[1] M. Forshaw _et al._, "Constraint Optimisation Approaches for Designing Group-Living Captive Breeding Programmes," 
presented at _39th Ann. AAAI Conf. on Artificial Intelligence_, Philadelphia, Pennsylvania, USA, 2025.
