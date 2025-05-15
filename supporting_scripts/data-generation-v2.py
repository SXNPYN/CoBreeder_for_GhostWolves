import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import pandas as pd
import random
import seaborn as sns

pd.set_option("display.max_columns", None)

num_individuals = 45
mean_alleles = 500
standard_deviation_alleles = mean_alleles/3
mean_pr = 500
standard_deviation_pr = mean_pr/3


# Individuals DataFrame.
individuals = pd.DataFrame(columns=["Name", "Male", "Female", "AssignToGroup", "Alleles", "Proven", "Priority"])
individuals['Male'] = [random.choice([0, 1]) for _ in range(num_individuals)]
individuals['Female'] = [1 - x for x in individuals['Male']]
for i in range(num_individuals):
    if individuals.loc[i]['Male'] == 1:
        individuals.loc[i, 'Name'] = f"Ind_{i}_M"
    else:
        individuals.loc[i, 'Name'] = f"Ind_{i}_F"
# if (random.choice([0, 1])) == 1:
    # assign_col = [-1] * num_individuals
# else:
    # assign_col = [-1] * (num_individuals - 1) + [1]
    # random.shuffle(assign_col)
alleles = abs((np.random.normal(mean_alleles, standard_deviation_alleles, num_individuals)).astype(int))
print(f"Mean (alleles): {np.mean(alleles)}")
print(f"Standard deviation (alleles): {np.std(alleles)}")
individuals['Alleles'] = alleles
individuals['AssignToGroup'] = assign_col = [-1] * num_individuals
individuals['Proven'] = [random.choice([0, 1]) for _ in range(num_individuals)]
individuals['Priority'] = 0
print(f"DATA PREVIEW:\n{individuals.head(10)}")

# PR matrix
random_matrix = abs(np.random.normal(mean_pr, standard_deviation_pr,
                                     size=(num_individuals, num_individuals)).astype(int))
random_matrix = (random_matrix + random_matrix.T) / 2
print(f"Mean (PR): {np.mean(random_matrix)}")
print(f"Standard deviation (PR): {np.std(random_matrix)}")
np.fill_diagonal(random_matrix, 0)
pr = pd.DataFrame(random_matrix, columns=individuals['Name'])
pr = pr.astype(int)
print(f'\n{pr.head(10)}')

# Save data.
data_dir_path = "../data/randomly_generated_data"
os.makedirs(data_dir_path, exist_ok=True)
individuals.to_csv(os.path.join(data_dir_path, "individuals.csv"), index=False)
pr.to_csv(os.path.join(data_dir_path, "pr-scaled.csv"), index=False)

# Histograms for alleles and PR
ax = sns.histplot(data=alleles, bins=num_individuals, color='steelblue', kde=True)
ax.set(xlabel='Alleles')
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
figure = ax.get_figure()
figure.savefig("../data/randomly_generated_data/allele_distribution.svg")
plt.close()
ax = sns.histplot(data=random_matrix.flatten(), bins=num_individuals, color='steelblue', kde=True)
ax.set(xlabel='Pairwise Relatedness')
figure = ax.get_figure()
figure.savefig("../data/randomly_generated_data/pr_distribution.svg")
plt.close()

print(f"\nRandom data generated in {data_dir_path}.")
