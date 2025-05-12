import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import seaborn as sns

sns.set_theme(context='paper', style='whitegrid', palette='Spectral')


def generate_graph(graph):
    if graph == '5.1':
        df = pd.read_csv("../results/5-1_averages.csv", na_values='INFEASIBLE')
        df = df.dropna(subset='avg_across_x_runs')
        df['avg_across_x_runs'] = df['avg_across_x_runs'].astype(int)
        ax = sns.lineplot(data=df, x="percent_allocated", y="av_time_to_first_solution_s",
                          hue='avg_across_x_runs', palette='Spectral', marker='D')
        ax.set(xlabel='Percentage of individuals allocated to pairings',
               ylabel='Average time to identify first feasible solution /s')
        ax.xaxis.set_major_formatter(ticker.PercentFormatter())
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
        plt.xlim(0, 100)
        ax.axhline(y=60, color='gray', linestyle='--')  # 1 min mark
        ax.legend(title="Runs averaged")
        figure = ax.get_figure()
        figure.savefig("../results/percentage_allocated_graph.svg")

    if graph == '5.2_num_sols':
        df = pd.read_csv("../results/5-2_spread_averages.csv", na_values='INFEASIBLE')

        ax = sns.lineplot(data=df, x="spread", y="av_num_solutions", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)',
                ylabel='Average number of solutions across 5 runs')
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.legend(title="% Allocated")
        plt.xlim(0, 1050)
        figure = ax.get_figure()
        figure.savefig("../results/data_spread_num_sols.svg")

    if graph == '5.2_time':
        df = pd.read_csv("../results/5-2_spread_averages.csv", na_values='INFEASIBLE')
        ax = sns.lineplot(data=df, x="spread", y="av_time_to_first_sol_s", hue='percent_allocated', marker='D')
        #ax = sns.lineplot(data=df, x="spread", y="av_time_to_best_sol_s", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)', ylabel='Average time to first solution /s')
        plt.xlim(0, 1050)
        plt.ylim(0, 1.5)
        ax.legend(title="% Allocated")
        figure = ax.get_figure()
        figure.savefig("../results/data_spread_time_first.svg")


#generate_graph('5.1')
#generate_graph('5.2_num_sols')
#generate_graph('5.2_time')

