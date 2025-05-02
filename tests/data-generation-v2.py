import pandas as pd
import random
import numpy as np
import os

pd.set_option("display.max_columns", None)

NUM_INDIVIDUALS = 25  # Number of individuals to generate
MIN_ALLELES = 100  # Lower bound for ghost alleles (value included)
MAX_ALLELES = 800  # Upper bound for ghost alleles (value included)
MIN_PR = 100  # Lower bound for PR permitted between two different individuals (value included)
MAX_PR = 1000  # Upper bound for PR permitted between two different individuals (value included)


# Generate individuals.
individuals = pd.DataFrame(columns=["Name", "Male", "Female", "AssignToGroup", "Alleles", "Proven", "Priority"])
individuals['Male'] = [random.choice([0, 1]) for _ in range(NUM_INDIVIDUALS)]
individuals['Female'] = [1 - x for x in individuals['Male']]
for i in range(NUM_INDIVIDUALS):
    if individuals.loc[i]['Male'] == 1:
        individuals.loc[i, 'Name'] = f"Ind_{i}_M"
    else:
        individuals.loc[i, 'Name'] = f"Ind_{i}_F"
if (random.choice([0, 1])) == 1:
    assign_col = [-1] * NUM_INDIVIDUALS
else:
    assign_col = [-1] * (NUM_INDIVIDUALS - 1) + [1]
    random.shuffle(assign_col)
individuals['AssignToGroup'] = assign_col
individuals['Alleles'] = [random.randint(MIN_ALLELES, MAX_ALLELES) for _ in range(NUM_INDIVIDUALS)]
individuals['Proven'] = [random.choice([0, 1]) for _ in range(NUM_INDIVIDUALS)]
individuals['Priority'] = 0
print("DATA PREVIEW:")
print(f'\n{individuals.head(10)}')

# Generate PR matrix.
random_matrix = np.random.randint(MIN_PR, MAX_PR, size=(NUM_INDIVIDUALS, NUM_INDIVIDUALS))
random_matrix = (random_matrix + random_matrix.T) / 2
np.fill_diagonal(random_matrix, 0)
pr = pd.DataFrame(random_matrix, columns=individuals['Name'])
pr = pr.astype(int)
print(f'\n{pr.head(10)}')

# Save data.
data_dir_path = "../data/randomly_generated_data"
os.makedirs(data_dir_path, exist_ok=True)
individuals.to_csv(os.path.join(data_dir_path, "generated-individuals.csv"), index=False)
pr.to_csv(os.path.join(data_dir_path, "generated-pr-scaled.csv"), index=False)
print(f"\nRandom data generated in {data_dir_path}.")
