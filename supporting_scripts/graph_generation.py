import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import seaborn as sns

sns.set_theme(context='paper', style='whitegrid', palette='Spectral')


# Graph visualising how the percentage of individuals allocated affects average time to first solution
df = pd.read_csv("../results/5-1_averages.csv", na_values='NaN')
ax = sns.lineplot(data=df, x="percent_allocated", y="av_time_to_first_solution_s",
                  hue='avg_across_x_runs', color='b', palette='Spectral', marker='D')
ax.set(xlabel='Percentage of individuals allocated to pairings', ylabel='Average time to identify first feasible '
                                                                        'solution /s')
ax.xaxis.set_major_formatter(ticker.PercentFormatter())
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
plt.xlim(0, 100)
ax.axhline(y=60, color='g', linestyle='--')  # 1 min mark
plt.legend(title="Number of runs averaged across")
figure = ax.get_figure()
figure.savefig("../results/percentage_allocated_graph.svg")
