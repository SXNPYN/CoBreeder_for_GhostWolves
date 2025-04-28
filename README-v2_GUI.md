# CoBreeder for Ghost Wolves (v2) Usage Guidelines (GUI)

---

These usage guidelines pertain to the **Graphical User Interface** (GUI) offered by _CoBreeder for Ghost Wolves_. To 
use the **Command Line Interface** instead, please see `README.md` and `ghost-cobreeder-v2_CLI.py` in the project 
directory.


## Context

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
tool should take this into account, improving its flexibility. Secondly, the pairwise relatedness file currently 
assumes that values are scaled and a larger value means that two individuals are less related. This can be developed 
further to allow for users to specify more intuitive thresholds (e.g. third cousins). 


## Fields (not the grassy kind)

Unfortunately, these fields aren't particularly biodiverse and must be filled with alphanumeric characters rather than 
flora and fauna. Some fields are mandatory but others may be left blank. Each field is explained below:

- `Individuals Specification File` 
  - Path to the CSV file detailing the individuals. This can be directly pasted into the field, but it is recommended to
    use the "Upload CSV" button to browse your files locally and select it from there. 
  - Please see the "Expected File Formats" section for more details. 

- `Scaled Pairwise Relatedness File`
  - Path to the CSV file containing the scaled pairwise relatedness matrix. This can be directly pasted into the field, 
  - but it is recommended to use the "Upload CSV" button to browse your files locally and select it from there. 
  - Please see the "Expected File Formats" section for more details. 

- `Number of Pairings`
  - The number of pairings to allocate. This should be an integer (e.g. `4`).

- `Objective Function`
  - Specifies the objective function to use when solving. The options offered are:
    - `MIN_AV_PR` - Minimise the average pairwise relatedness across pairings.
    - `MAX_TOTAL_ALLELES` - Maximise the total number of ghost alleles in a solution.
    - `MAX_TOTAL_PRIO` - Maximise the total priority in a solution.
    - `MIN_PR_MAX_ALLELES` - Minimise pairwise relatedness and maximise the number of ghost alleles in each pairing 
    using the weights specified. 
    - `MIN_PR_MAX_ALLELES_MAX_PRIO` - Minimise pairwise relatedness, maximise the number of ghost alleles, and 
    maximise priority in each pairing using the weights specified.
  - Note that `MIN_PR_MAX_ALLELES_MAX_PRIO` and `MAX_PRIO` cannot be used if `prio_calc_threshold` is set to 0 as they
    use dynamic priority calculations. 

- `Weight (Alleles,PR,Priority)`
  - List of integers specifying the relative weight to place on alleles, pairwise relatedness and priority values when 
    performing multi-objective optimisation.
  - You must specify these in this order, separated by commas. Do not use other separating characters.
    - Example: Allele weight = 3, PR weight = 2, Priority weight = 0 -> `3,2,0`

- `Unique Run ID`
  - Unique string identifier for each run. Has no bearing on functionality but is useful to identify results when 
    running the program multiple times.

- `Global PR Threshold`
  - Threshold for the scaled pairwise relatedness permitted in a pairing (0 by default).
  - Value not included.

- `Custom PR Thresholds`
  - This field is used to specify custom PR thresholds for certain groups.
  - It requires a list of group-PR pairs, separated by commas. Please remember that PR values are scaled integers. You 
    will need to use the group ID (index), which starts at zero. 
    - Example: If you would like the first group to have a custom PR threshold of 42 -> `0-42`
      - To extend this so that the 5th group has a PR threshold of 100 -> `0-42,4-100`
  - You do not need to specify a threshold for every group. Non-specified groups will take the global PR threshold, so
    this field can be left entirely blank.

- `Exclusions`
  - Exclude individuals by specifying a list of their IDs, separated by commas (e.g. `1,4,12`)
    - Use the `Preview Individuals` button to help when referencing indices.
  - This is explained in more detail in the "Excluded Individuals & Disallowed Pairings" section.
  - This field can be left blank.

- `Disallowed Pairings`
  - Define disallowed pairings by specifying a list of ID pairs, separated by commas (e.g. `2-5, 7-0`).
  - This is explained in more detail in the "Excluded Individuals & Disallowed Pairings" section.
  - This field can be left blank.

- `Size of Priority Set`
  - Threshold for priority calculations representing the number of individuals that can fall into the priority set.
    - E.g. 4 will enable priority calculations and select the top 4 scoring individuals for the priority set.
  - Must be a positive integer.
  - Specify 0 to disable dynamic priority calculations and use manual priority assignments only.
  - This is explained in more detail in the "Priority Calculations" section.

- `Weight for Alleles (Priority Calculations)`
  - Specifies the weight placed on ghost alleles when performing priority calculations. 
  - Must be a positive float between 0 and 1.
  - This is explained in more detail in the "Priority Calculations" section.

- `Save final solution to CSV?`
  - Selecting YES will save the best solution found to a CSV file in the `results` folder.


## Expected File Formats

Input files should be CSVs with the column names specified below. The program expects two files, detailing the 
individuals and the pairwise relatedness between them. Please ensure that these files are complete. If the pairwise 
relatedness between two individuals is unknown, marking these cells with 0 will ensure that the individuals are not 
paired together. 

Note that data can be randomly generated for testing purposes using the script in `tests/data-generation-v2.py`, 
specifying number of individuals and lower/upper bounds for PR and ghost alleles. 

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

- Every individual in the dataset will be assigned an index, starting from 0. The `Preview Individuals` button will 
display the individuals file in table format to help you when referencing indices. 
- Disallowed pairings:
  - To disallow pairings between individual 4 and 8, and 0 and 1 you will type: `4-8, 0-1`
  - Order does not matter here. For example, you could also type: `1-0, 8-4`
  - This field can be left blank.
- Exclusions
  - To exclude individuals 2, 5, and 6 you will type: `2,5,6`
  - Again, order does not matter here, and this field can be left blank.

Note: neither of these options edit the original files. Your original data will be left untouched and this will only
apply to the current run.

## Priority Calculations

Previous tools allow breeding program managers to manually specify "priority individuals", individuals that have 
desirable characteristics and must be included in solutions. This program allows for this but also builds upon this
idea with dynamic priority calculations that allow individuals to be ranked by priority. Individuals are compared to
their peers and assigned a priority score between 0 and 100, where 100 is the maximum priority value.

The top x individuals are selected to fall into the priority set, representing individuals that must be included in 
solutions, whilst the more complex priority values can be used in some objective functions.

Priority calculations consider various factors including whether individuals are proven and the number of mates and 
ghost alleles that they have compared to their competitors. If priority calculations are enabled, you will need to 
specify a weight for ghost alleles. This is a float between 0 and 1, where 1 represents completely ignoring 
number of mates in favour of ghost alleles, and 0 means only considering number of mates and ignoring the number
of ghost alleles. 0.5 will strike a balance between these. It is recommended to prioritise ghost alleles more heavily
at the start and increase emphasis on the number of mates if there are unforeseen circumstances (e.g. individuals are
proving hard to capture). This will result in more flexible solutions.

## Outputs

Each feasible solution found will be printed to the terminal as the code is running. Once the search has finished, 
you will be prompted to choose whether to save the best solution to a CSV file. If this option is chosen, the solution 
with the best objective value will be saved to the `results` folder, with a file name following the format
`best_solution_(unique_run_id).csv`. This will also be confirmed in the terminal display. 

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
