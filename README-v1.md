# CoBreeder for Ghost Wolves (v1) Usage Guidelines

---

This project builds upon the work completed by _Forshaw et al._, titled _Constraint Optimisation Approaches for 
Designing Group-Living Captive Breeding Programmes_ [1]. The tool produced (_CoBreeder_) was designed with a focus on 
the Galápagos tortoise which required single-objective optimisation to minimise pairwise relatedness when allocating 
individuals to breeding groups. The aim of this project is to adapt and extend _CoBreeder_ to improve its suitability 
for breeding coyotes with red wolf ancestry. This requires multi-objective optimisation, minimising genetic relatedness 
and maximising the number of ghost alleles in each pairing. Additional functionality that could aid conservationists 
has also been implemented, such as the ability to exclude individuals and specify disallowed pairings without editing 
the original data. There is also the option to calculate a priority score for each individual based on the 
characteristics of the dataset. 

This initial adaptation focuses on assigning coyote pairings however this may not be entirely accurate to life as 
coyotes, though sometimes solitary animals, often live in packs with a single breeding pair. Future versions of this 
tool should take this into account, improving its flexibility. Secondly, there is room for improvement where 
accessibility is concerned, as some users may find the current system for defining arguments confusing and cluttered.
Thirdly, the pairwise relatedness file currently assumes that values are scaled and a larger value means that two 
individuals are less related. This can be developed further to allow for users to specify more intuitive thresholds
(e.g. third cousins). Finally, input validation should be improved to counter human error.


## Command Line Arguments

This program is run via the command line. You will need to specify various arguments as well as provide relative paths 
to the data that will be used (please see the "Expected File Formats" section for more details). The program expects 
the following arguments, in the order introduced:

- `individuals_file`
  - Relative path to the CSV file detailing the individuals.
- `pairwise_relatedness_file`
  - Relative path to the CSV file containing the scaled pairwise relatedness matrix.
- `group_file`
  - Relative path to the CSV file detailing the group specifications.
- `obj_function`
  - String specifying the objective function to use when solving. The options are:
    - `MIN_AV_PR` - Minimise the average pairwise relatedness across pairings.
    - `MAX_TOTAL_ALLELES` - Maximise the total number of ghost alleles in a solution.
    - `MAX_TOTAL_PRIO` - Maximise the total priority in a solution.
    - `MIN_PR_MAX_ALLELES` - Minimise pairwise relatedness and maximise the number of ghost alleles in each pairing 
    using the weights specified. 
    - `MIN_PR_MAX_ALLELES_MAX_PRIO` - Minimise pairwise relatedness, maximise the number of ghost alleles, and 
    maximise priority in each pairing using the weights specified.
  - Note that `MIN_PR_MAX_ALLELES_MAX_PRIO` and `MAX_PRIO` cannot be used if `prio_calc_threshold` is set to 0 as they
    use dynamic priority calculations. 
- `unique_run_id`
  - Unique string identifier for each run. Has no bearing on functionality but is useful to identify results when 
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
  - E.g. 4 will enable priority calculations and select the top 4 scoring individuals for the priority set.
  - This is explained in more detail in the "Priority Calculations" section.

### EXAMPLE:

- `poetry run python ghost_wolf_cobreeder/ghost-cobreeder-v1.py run data/individuals.csv data/pr-scaled.csv 
data/groups.csv MIN_PR_MAX_ALLELES ghost_experiment 1 1 0 0 ALL 4`
  - The name of the run is "ghost_experiment". This string will be displayed alongside results.
  - The solver will maximise the number of ghost alleles and minimise pairwise relatedness between individuals in each 
  pairing. Equal weight is placed on each of these.
  - Any individuals with a PR greater than 0 can be paired (note that this would almost certainly need to be adjusted
  in a real-life scenario).
  - All individuals in the individuals file will be considered, and all opposite sex pairings (provided PR > 0) are 
  allowed.
  - Priority calculations are enabled and the top 4 individuals with the best priority values must be included in 
    solutions.


## Expected File Formats

All input files should be in CSV format, with the column names specified below. The program expects three files, 
detailing the individuals, groups, and the pairwise relatedness between the individuals. Please ensure that these
files are complete. If the pairwise relatedness between two individuals is unknown, marking these cells with 0 will
ensure that the individuals are not paired together. 

Note that data can be randomly generated for testing purposes using the script in `tests/data-generation.py`, specifying
number of individuals, groups, and lower/upper bounds on PR and ghost alleles. 

### INDIVIDUAL SPECIFICATION FILE:

- `Name`
  - String describing name of the individual. May include various details such as sex and location.
- `Male` & `Female` 
  - Boolean columns taking values 0 or 1, where 1 indicates that the individual's sex matches the name of the column.
  - For any row, both columns cannot take the same value.
- `AssignToGroup` 
  - Specifies whether to assign this individual to a specific group.
  - This column should be -1 by default. If specifying a group, ensure that the value used matches a group ID.
- `Alleles` 
  - A positive integer representing the number of ghost alleles that the individual has. 
- `Proven` 
  - Boolean column taking values 0 or 1, where 1 indicates that the individual is proven.
- `Priority` 
  - Boolean column taking values 0 or 1, where 1 indicates that the individual is a "priority individual".
  - This column is used to manually specify priority individuals. 
  - Note that if priority calculations are enabled, the values in this column will not be considered. 

Example:

| Name           | Male | Female | AssignToGroup | Alleles | Proven | Priority |
|----------------|------|--------|---------------|---------|--------|----------|
| Individual_0_M | 1    | 0      | -1            | 295     | 1      | 0        |
| Individual_1_M | 1    | 0      | -1            | 312     | 0      | 0        |
| Individual_2_F | 0    | 1      | 1             | 488     | 1      | 1        |
| ...            | ...  | ...    | ...           | ...     | ...    | ...      |


### GROUP SPECIFICATION FILE:

- `ID` 
  - Unique group ID, usually equal to the row index.
- `MinSize` & `MaxSize`
  - Minimum/maximum number of individuals to assign to this group (values included).
  - For pairings, both of these columns should be set to 2. 
- `NumMale` & `NumFemale` 
  - Number of males/females allowed in the group. The sum of these should not exceed the maximum capacity of the group.
  - For pairings, both of these columns should be set to 1. 
- `PRThreshold` 
  - The scaled PR threshold allowed in this group. 
  - Set to -1 to use the default value. 

All values in this file should be integers. 

Example:

| ID  | MinSize | MaxSize | NumMale | NumFemale | PRThreshold |
|-----|---------|---------|---------|-----------|-------------|
| 0   | 2       | 2       | 1       | 1         | -1          |
| 1   | 2       | 2       | 1       | 1         | -1          |
| 2   | 2       | 2       | 1       | 1         | -1          |
| ... | ...     | ...     | ...     | ...       | ...         |


### SCALED PR SPECIFICATION FILE:

This file represents the pairwise relatedness matrix and must include values for each individual in the individuals
file. Note that the values in this file are scaled; 0 represents the same individual and is also used to specify
banned pairings. A larger value indicates that the individuals are less related. This may not seem intuitive,
but when minimising pairwise relatedness the program is actually maximising these values. 

It is important that every individual specified in the individual specification file is included here, and that the 
values are mirrored. There are two cells for each pairing, therefore if editing the PR between two different 
individuals you will need to edit two cells.

Note that you do not need to manually edit the file to exclude individuals or create disallowed pairings, this can be
done whilst running the program. 

Example: 

| Individual_0_M | Individual_1_F | Individual_2_F | Individual_3_M | Individual_4_M |
|----------------|----------------|----------------|----------------|----------------|
| 0              | 217            | 749            | 474            | 538            |
| 217            | 0              | 581            | 798            | 616            |
| 749            | 581            | 0              | 207            | 542            |
| 474            | 798            | 207            | 0              | 368            |
| 538            | 616            | 542            | 368            | 0              |


## Excluded Individuals & Disallowed Pairings

1. If "EX" was specified on the command line, you will be shown a preview of the individuals file in table format. Each 
individual will have an ID on the far left-hand side, starting at 0.
2. You will first be prompted to specify disallowed pairings. 
   - To disallow pairings between individual 4 and 8, and 0 and 1 you will type: 4-8, 0-1
   - Order does not matter here. For example, you could also type: 1-0, 8-4 
   - To skip this, leave it blank and tap the enter key.
3. You will then be prompted to specify individuals to exclude.
   - To exclude individuals 2, 5, and 6 you will type: 2,5,6
   - Again, order does not matter here, and you can skip this by leaving it blank and tapping the enter key.

Note: neither of these options edit the original files. Your original data will be left untouched and this will only
apply to the current run.


## Priority Calculations

Previous tools allow breeding program managers to manually specify "priority individuals", individuals that have 
desirable characteristics and must be included in solutions. This program allows for this but also builds upon this
idea with dynamic priority calculations that allow individuals to be ranked by priority. Individuals are compared to
their peers and assigned a priority score between 0 and 100, where 100 is the maximum priority value.

The top x individuals are selected to fall into the priority set, representing individuals that must be included in 
solutions, whilst the more complex priority values can be used in some objective functions.
(NOTE: The PR threshold for priority calculations was removed as this didn't feel intuitive in practice.)

Priority calculations consider various factors including whether individuals are proven and the number of mates and 
ghost alleles that they have compared to their competitors. If priority calculations are enabled, you will be prompted 
to provide a weight for ghost alleles. This is a float between 0 and 1, where 1 represents completely ignoring 
number of mates in favour of ghost alleles, and 0 means only considering number of mates and ignoring the number
of ghost alleles. 0.5 will strike a balance between these. It is recommended to prioritise ghost alleles more heavily
at the start and increase emphasis on the number of mates if there are unforeseen circumstances (e.g. individuals are
proving hard to capture). This will result in more flexible solutions.

## Outputs

Each feasible solution found will be printed to the terminal as the code is running. Once the search has finished, 
you will be prompted to choose whether to save the best solution to a CSV file. If this option is chosen, the solution 
with the best objective value will be saved to the `results` folder, with a file name following the format
`best_solution_(unique_run_id).csv`. This will also be confirmed in the terminal. 

This file may look something like this:

| Group | Ind_1_Name      | Ind_2_Name      | Ind_1_ID | Ind_2_ID | Ind_1_Alleles | Ind_2_Alleles | Pairwise_Relatedness | Priority_Sum |
|-------|-----------------|-----------------|----------|----------|---------------|---------------|----------------------|--------------|
| 0     | Individual_38_F | Individual_41_M | 38       | 41       | 357           | 155           | 967                  | 129          |
| 1     | Individual_22_F | Individual_37_M | 22       | 37       | 781           | 281           | 877                  | 161          |
| 2     | Individual_7_F  | Individual_17_M | 7        | 17       | 261           | 482           | 967                  | 143          |
| 3     | Individual_27_F | Individual_42_M | 27       | 42       | 401           | 440           | 927                  | 73           |
| ...   | ...             | ...             | ...      | ...      | ...           | ...           | ...                  | ...          |

Each group is represented by a row, which contains information on the names, IDs, and number of ghost alleles of the 
individuals in the pairing as well as their pairwise relatedness. The `Priority_Sum` column is the sum of the priority
values of each individual and can reach a maximum of 200. Note that, if priority calculations are not enabled, the 
cells in this column will be "N/A".

---

## References

[1] M. Forshaw et al., "Constraint Optimisation Approaches for Designing Group-Living Captive Breeding Programmes," 
presented at 39th Ann. AAAI Conf. on Artificial Intelligence, Philadelphia, Pennsylvania, USA, 2025.
