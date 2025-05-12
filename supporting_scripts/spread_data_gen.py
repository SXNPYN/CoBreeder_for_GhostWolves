import pandas as pd
import random
import matplotlib.ticker as ticker
import numpy as np
import os
import seaborn as sns

pd.set_option("display.max_columns", None)
NUM_INDIVIDUALS = 25
# VALUES = [i for i in range(1, 1001) if i % 40 == 0] # Evenly distributed
VALUES = [1, 222, 278, 315, 345, 370, 392, 412, 431, 449, 467, 484, 500, 517, 534, 552, 570,
          589, 609, 631, 656, 686, 723, 779, 1000]  # Normal distribution
random.shuffle(VALUES)
random_matrix = np.random.randint(1, 1000, size=(NUM_INDIVIDUALS, NUM_INDIVIDUALS))

ax = sns.histplot(data=VALUES, bins=25, color='steelblue', kde=True)
ax.set(xlabel='Values')
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
figure = ax.get_figure()
figure.savefig("../results/normal_distribution.svg")


def generate_data(scale):
    individuals = pd.DataFrame(columns=["Name", "Male", "Female", "AssignToGroup", "Alleles", "Proven", "Priority"])
    individuals['Male'] = [random.choice([0, 1]) for _ in range(NUM_INDIVIDUALS)]
    individuals['Female'] = [1 - x for x in individuals['Male']]
    for i in range(NUM_INDIVIDUALS):
        if individuals.loc[i]['Male'] == 1:
            individuals.loc[i, 'Name'] = f"Ind_{i}_M"
        else:
            individuals.loc[i, 'Name'] = f"Ind_{i}_F"
    individuals['AssignToGroup'] = [-1] * NUM_INDIVIDUALS
    individuals['Alleles'] = VALUES
    individuals['Alleles'] = individuals['Alleles'].apply(lambda x: int(x * scale))
    individuals['Proven'] = [random.choice([0, 1]) for _ in range(NUM_INDIVIDUALS)]
    individuals['Priority'] = 0
    print("DATA PREVIEW:")
    print(f'\n{individuals.head(10)}')

    matrix = random_matrix * scale
    matrix = (matrix + matrix.T) / 2
    np.fill_diagonal(matrix, 0)
    pr = pd.DataFrame(matrix, columns=individuals['Name'])
    pr = pr.astype(int)
    print(f'\n{pr.head(10)}')

    # Save data.
    data_dir_path = "../data/dissertation_5.1/spread_investigations"
    os.makedirs(data_dir_path, exist_ok=True)
    individuals.to_csv(os.path.join(data_dir_path, f"norm_dist_individuals_{scale}.csv"), index=False)
    pr.to_csv(os.path.join(data_dir_path, f"norm_dist_pr-scaled_{scale}.csv"), index=False)
    print(f"\nRandom data generated in {data_dir_path}.")


generate_data(1)
generate_data(0.85)
generate_data(0.7)
generate_data(0.55)
generate_data(0.4)
generate_data(0.25)
generate_data(0.1)
