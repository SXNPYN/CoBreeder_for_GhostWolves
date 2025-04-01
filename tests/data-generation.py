import pandas as pd
import random
import numpy as np
import os

NUM_INDIVIDUALS = 50
MIN_ALLELES = 150
MAX_ALLELES = 900
NUM_GROUPS = 10
MIN_PR = 1
MAX_PR = 1000


# Generate groups
groups = pd.DataFrame(columns=["ID", "MinSize", "MaxSize", "NumMale", "NumFemale", "PRThreshold"])
groups['ID'] = [i for i in range(NUM_GROUPS)]
groups['MinSize'] = [2]*NUM_GROUPS
groups['MaxSize'] = [2]*NUM_GROUPS
groups['NumMale'] = [1]*NUM_GROUPS
groups['NumFemale'] = [1]*NUM_GROUPS
groups['PRThreshold'] = [-1]*NUM_GROUPS
print(f'\n{groups.head()}')

# Generate individuals
individuals = pd.DataFrame(columns=["Name", "Male", "Female", "AssignToFirstGroup", "Alleles", "Proven", "Priority"])
individuals['Name'] = [f"Individual_{i+1}" for i in range(NUM_INDIVIDUALS)]
individuals['Male'] = [random.choice([0, 1]) for _ in range(NUM_INDIVIDUALS)]
individuals['Female'] = [1 - x for x in individuals['Male']]
individuals['AssignToFirstGroup'] = -1
individuals['Alleles'] = [random.randint(MIN_ALLELES, MAX_ALLELES) for _ in range(NUM_INDIVIDUALS)]
individuals['Proven'] = [random.choice([0, 1]) for _ in range(NUM_INDIVIDUALS)]
individuals['Priority'] = 0
print(f'\n{individuals.head()}')

# Generate PR matrix
random_matrix = np.random.randint(MIN_PR, MAX_PR, size=(NUM_INDIVIDUALS, NUM_INDIVIDUALS))
random_matrix = (random_matrix + random_matrix.T) / 2
np.fill_diagonal(random_matrix, 0)
pr = pd.DataFrame(random_matrix, columns=individuals['Name'])
pr = pr.astype(int)
print(f'\n{pr.head()}')

# Save data
data_dir_path = "data/wolves_randomly_generated"
os.makedirs(data_dir_path, exist_ok=True)
groups.to_csv(os.path.join("../", data_dir_path, "generated-individuals.csv"), index=False)
individuals.to_csv(os.path.join("../", data_dir_path, "generated-groups.csv"), index=False)
pr.to_csv(os.path.join("../", data_dir_path, "generated-pr-scaled.csv"), index=False)
print("Random data generated in data/wolves_randomly_generated.")
