# CoBreeder for Ghost Wolves (V1) Usage Guidelines

---

Explain what it is and does.
Builds upon CoBreeder for Galapagos tortoises, reference Matt's work.
PR of 0 means not compatible 

## Inputs

The program is run via the terminal/command line. You will need to specify various arguments as well as provide paths 
to the files to use. 

### COMMAND LINE ARGUMENTS
```python
# Insert code here
```

## Expected File Formats

All input files should be in CSV format, with the column names specified below. The program expects three files, 
detailing the individuals, groups, and the pairwise relatedness between the individuals. Please ensure that these
files are complete. If the pairwise relatedness between two individuals is unknown, assigning these cells with 0 will
ensure that the individuals are not paired together. 

### INDIVIDUAL SPECIFICATION FILE:

- `Name`
  - String describing name of the individual. May include various details such as sex and location.
- `Male` & `Female` 
  - Boolean columns taking values 0 or 1, where 1 indicates that the individual's sex matches the name of the column.
  - For any row, both columns cannot take the same value.
- `AssignToGroup` 
  - Specifies whether to assign this individual to a specific group.
  - This column should be -1 by default.
  - 
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

This file represents the pairwise relatedness matrix, and must include values for each individual in the individuals
file. Note that the values in this file are scaled; 0 represents the same individual, and is also used to specify
banned pairings. A larger value indicates that the individuals are less related. This may not seem intuitive at first,
but when minimising pairwise relatedness, the program is actually maximising the values in this file. 

It is important that every individual specified in the individual specification file is included here, and that the 
values are mirrored. There are two cells for each pairing, therefore if editing the PR between two different 
individuals, you will need to edit two cells.

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


### Recommended weighting values & inputs whilst the program runs 

		- Recommend to not prioritise number of mates heavily at the start
		- Weightings etc
		- Recommend not to use prio optimisation if not done calculation

## Outputs

Each feasible solution found will be printed to the terminal as the code is running. Once the search has finished, 
you will be prompted to choose whether to save the best solution to a CSV file. If this option is chosen, the solution 
with the best objective value will be saved to the `results` folder, with a file name following the format
`best_solution_(experiment name).csv`. This will also be confirmed in the terminal. 

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
values of each individual, and can reach a maximum of 200. Note that, if priority calculations are not enabled, the 
cells in this column will be "N/A".
