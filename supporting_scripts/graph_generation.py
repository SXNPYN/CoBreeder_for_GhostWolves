import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import seaborn as sns

sns.set_theme(context='paper', style='whitegrid', palette='Spectral')


def generate_graph(graph):

    if graph == '5.1':

        df = pd.read_csv("../results/5-1_averages.csv", na_values='INFEASIBLE')
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
        figure.savefig("../results/percentage_allocated_mean.svg")
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
        figure.savefig("../results/percentage_allocated_mean_med.svg")
        plt.close()

    if graph == '5.2_uniform':

        df = pd.read_csv("../results/5-2_spread_averages.csv", na_values='INFEASIBLE')
        df = df[df['distribution'] == 'uniform']

        ax = sns.lineplot(data=df, x="spread", y="mean_num_solutions", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)',
               ylabel='Mean number of solutions across 5 runs',
               title='(Uniform distribution)')
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.legend(title="% Allocated")
        plt.xlim(0, 1050)
        figure = ax.get_figure()
        figure.savefig("../results/num_sols_uniform_dist.svg")
        plt.close()

        ax = sns.lineplot(data=df, x="spread", y="mean_time_to_first_sol_s", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)',
               ylabel='Mean time to first solution /s',
               title='(Uniform distribution)')
        plt.xlim(0, 1050)
        plt.ylim(0, 1.1)
        ax.legend(title="% Allocated")
        figure = ax.get_figure()
        figure.savefig("../results/time_first_uniform_dist.svg")
        plt.close()

        ax = sns.lineplot(data=df, x="spread", y="mean_time_to_best_sol_s", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)',
               ylabel='Mean time to best solution /s',
               title='(Uniform distribution)')
        plt.xlim(0, 1050)
        plt.ylim(0, 20)
        ax.legend(title="% Allocated")
        figure = ax.get_figure()
        figure.savefig("../results/time_best_uniform_dist.svg")
        plt.close()

    if graph == '5.2_normal':

        df = pd.read_csv("../results/5-2_spread_averages.csv", na_values='INFEASIBLE')
        df = df[df['distribution'] == 'normal']

        ax = sns.lineplot(data=df, x="spread", y="mean_num_solutions", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)',
               ylabel='Mean number of solutions across 5 runs',
               title='(Normal distribution)')
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.legend(title="% Allocated")
        plt.xlim(0, 1050)
        figure = ax.get_figure()
        figure.savefig("../results/num_sols_normal_dist.svg")
        plt.close()

        ax = sns.lineplot(data=df, x="spread", y="mean_time_to_first_sol_s", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)',
               ylabel='Mean time to first solution /s',
               title='(Normal distribution)')
        plt.xlim(0, 1050)
        plt.ylim(0,)
        ax.legend(title="% Allocated")
        figure = ax.get_figure()
        figure.savefig("../results/time_first_normal_dist.svg")
        plt.close()

        ax = sns.lineplot(data=df, x="spread", y="mean_time_to_best_sol_s", hue='percent_allocated', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)', ylabel='Mean time to best solution /s',
               title='(Normal distribution)')
        plt.xlim(0, 1050)
        # plt.ylim(0, 20)
        ax.legend(title="% Allocated")
        figure = ax.get_figure()
        figure.savefig("../results/time_best_normal_dist.svg")
        plt.close()

    if graph == '5.2_both':

        df = pd.read_csv("../results/5-2_spread_averages.csv", na_values='INFEASIBLE')
        df = df[df['percent_allocated'] == 26.7]

        ax = sns.lineplot(data=df, x="spread", y="mean_num_solutions", hue='distribution', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)',
               ylabel='Mean number of solutions across 5 runs')
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.legend(title="Distribution")
        plt.xlim(0, 1050)
        figure = ax.get_figure()
        figure.savefig("../results/num_sols_dist_comparison.svg")
        plt.close()

        ax = sns.lineplot(data=df, x="spread", y="mean_time_to_first_sol_s", hue='distribution', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)', ylabel='Mean time to first solution /s')
        plt.xlim(0, 1050)
        ax.legend(title="Distribution")
        figure = ax.get_figure()
        figure.savefig("../results/time_first_dist_comparison.svg")
        plt.close()

        ax = sns.lineplot(data=df, x="spread", y="mean_time_to_best_sol_s", hue='distribution', marker='D')
        ax.set(xlabel='Range (ghost alleles and pairwise relatedness)', ylabel='Mean time to best solution /s')
        plt.xlim(0, 1050)
        ax.legend(title="Distribution")
        figure = ax.get_figure()
        figure.savefig("../results/time_best_dist_comparison.svg")
        plt.close()


# generate_graph('5.1')
# generate_graph('5.2_uniform')
# generate_graph('5.2_normal')
generate_graph('5.2_both')
