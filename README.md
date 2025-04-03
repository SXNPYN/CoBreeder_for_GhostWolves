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

### Individual Specification File:

- `ID` - Unique group ID.
- `MinSize` - Minimum amount of individuals to assign to this group.
- `MaxSize` - Maximum capacity of the group.
- `NumMale` - Number of males allowed in the group.
- `NumFemale` - Number of females allowed in the group.
- `PRThreshold` - The scaled PR threshold allowed in this group. Set to -1 to use the default value. 

All values in these columns should be integers. 

Example:

| ID  | MinSize | MaxSize | NumMale | NumFemale | PRThreshold |
|-----|---------|---------|---------|-----------|-------------|
| 0   | 2       | 2       | 1       | 1         | -1          |
| 1   | 2       | 2       | 1       | 1         | -1          |
| 2   | 2       | 2       | 1       | 1         | -1          |
| ... | ...     | ...     | ...     | ...       | ...         |


### Group Specification File:

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


#### Pairwise Relatedness File:
- Col 1 Name
  - What it is
  - Data format
- Col 2

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
