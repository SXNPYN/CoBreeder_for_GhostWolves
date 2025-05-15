import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import seaborn as sns
import numpy as np

sns.set_theme(context='paper', style='whitegrid', palette='Spectral')


def generate_graph(graph):

    if graph == '5.1_percent_allocated':

        df = pd.read_csv("../results/percent_allocated/5-1_averages.csv", na_values='INFEASIBLE')
        df = df.dropna(subset='runs_averaged')
        df['runs_averaged'] = df['runs_averaged'].astype(int)
        ax = sns.lineplot(data=df, x="percent_allocated", y="mean_time_to_first_solution_s",
                          hue='runs_averaged', palette='Spectral', marker='D')
        ax.set(xlabel='Percentage of individuals allocated to pairings',
               ylabel='Mean time to identify first feasible solution /s')
        ax.xaxis.set_major_formatter(ticker.PercentFormatter())
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
        plt.xlim(0, 100)
        ax.axhline(y=60, color='gray', linestyle='--')  # 1 min mark
        ax.legend(title="Runs averaged")
        figure = ax.get_figure()
        figure.savefig("../results/percent_allocated/percentage_allocated_mean.svg")
        plt.close()

        ax = sns.lineplot(data=df, x="percent_allocated", y="mean_time_to_first_solution_s", marker='D', label='Mean')
        sns.lineplot(data=df, x="percent_allocated", y="median_time_to_first_solution_s", marker='D',
                     ax=ax, label='Median')
        ax.set(xlabel='Percentage of individuals allocated to pairings',
               ylabel='Time to identify first feasible solution /s')
        ax.xaxis.set_major_formatter(ticker.PercentFormatter())
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
        plt.xlim(0, 100)
        ax.axhline(y=60, color='gray', linestyle='--')  # 1 min mark
        ax.legend()
        figure = ax.get_figure()
        figure.savefig("../results/percent_allocated/percentage_allocated_mean_med.svg")
        plt.close()

    if graph == '5.2_uniform':

        df = pd.read_csv("../results/spread_and_distribution/5-2_spread_averages.csv",
                         na_values='INFEASIBLE')
        df = df[df['distribution'] == 'uniform']

        ax = sns.lineplot(data=df, x="spread", y="mean_num_solutions", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)',
               ylabel='Mean number of solutions across 5 runs',
               title='(Uniform distribution)\n')
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.legend(title="% Allocated")
        plt.xlim(0, 1050)
        figure = ax.get_figure()
        figure.savefig("../results/spread_and_distribution/num_sols_uniform_dist.svg")
        plt.close()

        ax = sns.lineplot(data=df, x="spread", y="mean_time_to_first_sol_s", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)',
               ylabel='Mean time to first solution /s',
               title='(Uniform distribution)\n')
        plt.xlim(0, 1050)
        plt.ylim(0, 1.1)
        ax.legend(title="% Allocated")
        figure = ax.get_figure()
        figure.savefig("../results/spread_and_distribution/time_first_uniform_dist.svg")
        plt.close()

        ax = sns.lineplot(data=df, x="spread", y="mean_time_to_best_sol_s", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)',
               ylabel='Mean time to best solution /s',
               title='(Uniform distribution)\n')
        plt.xlim(0, 1050)
        plt.ylim(0, 20)
        ax.legend(title="% Allocated")
        figure = ax.get_figure()
        figure.savefig("../results/spread_and_distribution/time_best_uniform_dist.svg")
        plt.close()

    if graph == '5.2_normal':

        df = pd.read_csv("../results/spread_and_distribution/5-2_spread_averages.csv",
                         na_values='INFEASIBLE')
        df = df[df['distribution'] == 'gaussian']

        ax = sns.lineplot(data=df, x="spread", y="mean_num_solutions", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)',
               ylabel='Mean number of solutions across 5 runs',
               title='(Normal distribution)\n')
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.legend(title="% Allocated")
        plt.xlim(0, 1050)
        figure = ax.get_figure()
        figure.savefig("../results/spread_and_distribution/num_sols_normal_dist.svg")
        plt.close()

        ax = sns.lineplot(data=df, x="spread", y="mean_time_to_first_sol_s", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)',
               ylabel='Mean time to first solution /s',
               title='(Normal distribution)\n')
        plt.xlim(0, 1050)
        plt.ylim(0,)
        ax.legend(title="% Allocated")
        figure = ax.get_figure()
        figure.savefig("../results/spread_and_distribution/time_first_normal_dist.svg")
        plt.close()

        ax = sns.lineplot(data=df, x="spread", y="mean_time_to_best_sol_s", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)', ylabel='Mean time to best solution /s',
               title='(Normal distribution)\n')
        plt.xlim(0, 1050)
        ax.legend(title="% Allocated")
        figure = ax.get_figure()
        figure.savefig("../results/spread_and_distribution/time_best_normal_dist.svg")
        plt.close()

    if graph == '5.2_dist_comparison':

        df = pd.read_csv("../results/spread_and_distribution/5-2_spread_averages.csv",
                         na_values='INFEASIBLE')
        df = df[df['percent_allocated'] == 26.7]

        ax = sns.lineplot(data=df, x="spread", y="mean_num_solutions", hue='distribution', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)',
               ylabel='Mean number of solutions across 5 runs')
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.legend(title="Distribution")
        plt.xlim(0, 1050)
        figure = ax.get_figure()
        figure.savefig("../results/spread_and_distribution/num_sols_dist_comparison.svg")
        plt.close()

        ax = sns.lineplot(data=df, x="spread", y="mean_time_to_first_sol_s", hue='distribution', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)', ylabel='Mean time to first solution /s')
        plt.xlim(0, 1050)
        ax.legend(title="Distribution")
        figure = ax.get_figure()
        figure.savefig("../results/spread_and_distribution/time_first_dist_comparison.svg")
        plt.close()

        ax = sns.lineplot(data=df, x="spread", y="mean_time_to_best_sol_s", hue='distribution', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)', ylabel='Mean time to best solution /s')
        plt.xlim(0, 1050)
        ax.legend(title="Distribution")
        figure = ax.get_figure()
        figure.savefig("../results/spread_and_distribution/time_best_dist_comparison.svg")
        plt.close()

    if graph == '5.3_objective_values':

        df = pd.read_csv("../results/scalability_results/5-3_scalability_runs.csv")
        df = df.replace(np.inf, 999)
        df = df[df['priority_set_size'] == 0]
        datasets = df['dataset'].unique()

        for d in datasets:
            filtered_df = df[df['dataset'] == d]
            ax = sns.lineplot(data=filtered_df, x="time", y="objective_value", hue='run', marker='D',
                              palette="Set1")
            ax.set(xlabel='Elapsed real time /s', ylabel='Objective value',
                   title=f'Dataset: {d}, size of priority set: 0\n')
            ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
            plt.legend(title='Run ID')
            if d == 'synthetic_30_normal_dist':
                plt.ylim(45, 100)
            elif d == 'synthetic_45_normal_dist':
                plt.ylim(35, 100)
            else:
                plt.ylim(40, 100)
            figure = ax.get_figure()
            figure.savefig(f"../results/scalability_results/no_priority_inds/objective_value_{d}.svg")
            plt.close()

    if graph == '5.3_scatter':

        df = pd.read_csv("../results/scalability_results/5-3_scalability_runs.csv")
        df = df[df['best'] == 1]
        df = df.replace({"synthetic_30_normal_dist": "30 individuals (synthetic)",
                         "synthetic_45_normal_dist": "45 individuals (synthetic)",
                         "synthetic_60_normal_dist": "60 individuals (synthetic)"})

        ax = sns.stripplot(data=df, x="dataset", y="time", hue="priority_set_size", dodge=True, jitter=False,
                           palette=['darkcyan', 'maroon'], marker="D")
        ax.set(xlabel='Dataset', ylabel='Time to best solution /s')
        ax.yaxis.set_major_locator(ticker.MultipleLocator(100))
        plt.legend(title='Size of priority set')
        figure = ax.get_figure()
        figure.savefig("../results/scalability_results/no_priority_inds/strip_plot.svg")
        plt.close()

        ax = sns.boxplot(data=df, x="dataset", y="time", hue="priority_set_size", width=.5,
                         palette=['darkcyan', 'maroon'], medianprops={"color":"w"})
        ax.set(xlabel='Dataset', ylabel='Time to best solution /s')
        ax.yaxis.set_major_locator(ticker.MultipleLocator(100))
        ax.legend(title="Size of priority set")
        figure = ax.get_figure()
        figure.savefig("../results/scalability_results/no_priority_inds/box_plot.svg")
        plt.close()

    if graph == '5.3_av_optimality':

        pass

    if graph == '5.3_prio_set':

        df = pd.read_csv("../results/scalability_results/5-3_scalability_runs.csv")
        no_priority = df[df['priority_set_size'] == 0]
        priority = df[df['priority_set_size'] == 2]

        pass


# generate_graph('5.1_percent_allocated')
# generate_graph('5.2_uniform')
# generate_graph('5.2_normal')
# generate_graph('5.2_dist_comparison')
# generate_graph('5.3_objective_values')
# generate_graph('5.3_scatter')
generate_graph('5.3_av_optimality')
generate_graph('5.3_prio_set')
