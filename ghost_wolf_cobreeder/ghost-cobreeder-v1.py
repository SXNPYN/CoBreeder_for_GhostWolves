#!/usr/bin/env python3
from absl import app
import argparse
from enum import IntEnum
from ortools.sat.python import cp_model
import os
import pandas as pd
import sys
import time
from typing import Sequence


class GhostCobreederObjectiveFunction(IntEnum):
    MIN_PR = 1
    MAX_ALLELES = 2
    MAX_PRIO = 3
    MIN_PR_MAX_ALLELES = 4
    MIN_PR_MAX_ALLELES_MAX_PRIO = 5


class GhostCobreederPrinter(cp_model.CpSolverSolutionCallback):
    def __init__(self, seats, names, num_groups, num_individuals, paramstring, unique_id):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__solution_count = 0
        self.__best_solution = {}
        self.__start_time = time.time()
        self.__seats = seats
        self.__names = names
        self.__num_groups = num_groups
        self.__num_individuals = num_individuals
        self.__paramstring = paramstring
        self.__uniqueid = unique_id

    def on_solution_callback(self):
        current_time = time.time()
        objective = self.ObjectiveValue()
        self.__solution_count += 1

        print(
            "\nCOBREEDER-SOLUTION #%i, time elapsed: %f, objective value: %i, paramstring: %s, experiment: %s"
            % (self.__solution_count,
               current_time - self.__start_time,
               objective, self.__paramstring, self.__uniqueid)
        )

        for t in range(self.__num_groups):
            print("\tGroup %d: " % t)
            individuals = []
            for g in range(self.__num_individuals):
                if self.Value(self.__seats[(t, g)]):
                    print(f"\t\t{self.__names[g]}")
                    individuals.append((g, self.__names[g]))
            self.__best_solution[t] = individuals

    def num_solutions(self):
        return self.__solution_count

    def best_solution(self):
        return self.__best_solution


def save_solution_csv(args, connections, individual_allele_count, individual_priority_value, best_solution):
    """
    Saves the best solution to a CSV file detailing groups, individuals, ghost alleles, and pairwise relatedness.

    Args:
        :param args: Variable length argument list.
        :param list of lists connections: PR relatedness matrix.
        :param dict individual_allele_count: Dictionary mapping individual IDs to their number of ghost alleles.
        :param dict individual_priority_value: Dictionary mapping individual IDs to their priority value.
        :param dict best_solution: Dictionary mapping group number to a list of tuples in the form (individual ID, name)
    """

    solution_data = pd.DataFrame(columns=["Group", "Ind_1_Name", "Ind_2_Name", "Ind_1_ID", "Ind_2_ID",
                                          "Ind_1_Alleles", "Ind_2_Alleles", "Pairwise_Relatedness", "Priority_Sum"])
    i = 0  # Index of current row
    for group, individuals in best_solution.items():
        solution_data.loc[i, "Group"] = group
        solution_data.loc[i, "Pairwise_Relatedness"] = connections[individuals[0][0]][individuals[1][0]]
        solution_data.loc[i, "Ind_1_ID"] = individuals[0][0]
        solution_data.loc[i, "Ind_2_ID"] = individuals[1][0]
        solution_data.loc[i, "Ind_1_Name"] = individuals[0][1]
        solution_data.loc[i, "Ind_2_Name"] = individuals[1][1]
        solution_data.loc[i, "Ind_1_Alleles"] = individual_allele_count[individuals[0][0]]
        solution_data.loc[i, "Ind_2_Alleles"] = individual_allele_count[individuals[1][0]]
        if args.prio_calc_threshold == 0:
            solution_data.loc[i, "Priority_Sum"] = "N/A"
        else:
            prio_sum = individual_priority_value[individuals[0][0]] + individual_priority_value[individuals[1][0]]
            solution_data.loc[i, "Priority_Sum"] = prio_sum
        i += 1

    # Create results directory if it doesn't exist and save CSV
    results_dir_path = os.path.join(os.getcwd(), "results")
    os.makedirs(results_dir_path, exist_ok=True)
    out_file = f"best_solution_{args.unique_run_id}.csv"
    solution_data.to_csv(os.path.join(results_dir_path, out_file), index=False)
    print("Solution saved to results/%s." % out_file)


def calculate_priority(args, individuals, pr):
    """
    Dynamically calculates a priority value for each individual based on the dataset between 0 and 100, where 100 is
    the highest priority value possible.

    Args:
        :param args: Variable length argument list.
        :param individuals: DataFrame containing information about each individual.
        :param pr: DataFrame containing pairwise relatedness information for each individual.
        :return individuals: Edited DataFrame detailing individuals and their calculated priorities.
    """

    if args.prio_calc_threshold == 0:  # Use values from csv only
        priorities = individuals["Priority"].tolist()
        individuals['PriorityValue'] = [100 * p for p in priorities]  # All priority individuals are equal

    else:  # Dynamically calculate priority                           ... but some are more equal than others.
        while True:
            a = input("Weight to place on ghost alleles: ")
            try:
                a = float(a)
                if 0 <= a <= 1:
                    break
            except ValueError:
                print("Please enter a number between 0.0 and 1.0")
        b = 1.0 - a
        print("\nPlacing weight on ghost alleles and number of mates in a ratio of %i:%i" % (10 * a, 10 * b))

        female_individuals = individuals.query("Female == 1").copy()
        male_individuals = individuals.query("Male == 1").copy()
        female_g_max = max(female_individuals["Alleles"])  # Largest number of ghost alleles amongst females
        male_g_max = max(male_individuals["Alleles"])
        female_m_max = len(male_individuals)  # Max number of potential mates that a female can have
        male_m_max = len(female_individuals)

        # Calculate number of potential mates each individual has from the PR matrix
        pr_threshold = args.pr_threshold
        for i in female_individuals.index:
            num_mates = (pr.loc[male_individuals.index][i] > pr_threshold).sum()
            female_individuals.loc[i, "NumMates"] = num_mates
            print(f'Female {i} has {num_mates} potential mates.')
        for i in male_individuals.index:
            num_mates = (pr.loc[female_individuals.index][i] > pr_threshold).sum()
            male_individuals.loc[i, "NumMates"] = num_mates
            print(f'Male {i} has {num_mates} mates.')

        # Calculate priority value between 0 and 100.
        female_individuals['PriorityValue'] = (female_individuals['Proven'] *
                                               (((a * female_individuals['Alleles']) / female_g_max) +
                                               ((b * female_individuals['NumMates']) / female_m_max)) * 100).astype(int)
        male_individuals['PriorityValue'] = (male_individuals['Proven'] *
                                             (((a * male_individuals['Alleles']) / male_g_max) +
                                             ((b * male_individuals['NumMates']) / male_m_max)) * 100).astype(int)

        individuals = pd.concat([male_individuals, female_individuals]).sort_index()
        individuals['NumMates'] = individuals['NumMates'].astype(int)

        individuals['Priority'] = 0  # 0 by default
        sorted_individuals = individuals.sort_values(by='PriorityValue', ascending=False)
        top_priority = sorted_individuals.head(args.prio_calc_threshold)
        individuals.loc[top_priority.index, 'Priority'] = 1  # Select the top x individuals to be priority individuals

    return individuals


def build_data(args):
    """
    Takes the CSV files (detailing individuals, groups, and pairwise relatedness) provided by the user and converts the
    data into a format that can be used by solve_model(). Also allows for exclusion of specific individuals and
    disallowed pairings, if enabled from the command line.

    Args:
        :param args: Variable length argument list.
        :return: Various variables storing data from the CSV files, mostly columns in list format.
    """

    objective_function = GhostCobreederObjectiveFunction[args.obj_function]
    unique_id = args.unique_run_id
    print("\nRUN ID: %s \nOBJECTIVE FUNCTION: %i" % (unique_id, objective_function))

    pr = pd.read_csv(args.pairwise_relatedness_file, delimiter=',', header=None, skiprows=1)
    print("USING PAIRWISE RELATEDNESS FILE [%s]" % args.pairwise_relatedness_file)

    group_defs = pd.read_csv(args.group_file, delimiter=',')
    print(f"\nGROUP DEFINITIONS [{args.group_file}]: \n{group_defs}")

    individuals = pd.read_csv(args.individuals_file, delimiter=',')

    if args.exclude_disallow == "EX":
        # Show data to help practitioner know which indices to use
        print(f"\nSUMMARY OF INDIVIDUALS [{args.individuals_file}]: \n{individuals}")

        while True:
            disallowed_pairings = input("Specify disallowed parings? (List of ID pairs e.g. 3-5, 2-6, 1-7): ")
            if disallowed_pairings:
                try:
                    disallowed_pairings = [tuple(map(int, x.split('-'))) for x in
                                           disallowed_pairings.strip().split(",")]
                    # Set PR for disallowed combinations to 0
                    for i, j in disallowed_pairings:
                        pr.iloc[i, j] = 0
                        pr.iloc[j, i] = 0
                    break
                except (IndexError, ValueError):
                    print("Please ensure indices are valid (e.g. 0-1, 4-2).")
            else:
                break

        while True:
            exclusions = input("Exclude individuals? (List of IDs e.g. 0, 4, 6): ")
            if exclusions:
                try:
                    # Remove excluded individuals from individuals data
                    exclusions = list(set(int(x) for x in exclusions.strip().split(",")))
                    individuals.drop(exclusions, axis=0, inplace=True)
                    print(f"\n UPDATED SUMMARY OF INDIVIDUALS [{args.individuals_file}]: \n{individuals}")
                    # Update PR matrix to remove excluded individuals
                    pr.drop(exclusions, axis=0, inplace=True)
                    pr.drop(exclusions, axis=1, inplace=True)
                    break
                except (KeyError, ValueError):
                    print("Invalid input. Please enter the indices of the individuals to exclude (e.g. 0, 4, 2)")
            else:
                break

    individuals = calculate_priority(args, individuals, pr)
    print(f"\nSUMMARY OF INDIVIDUALS AFTER PROCESSING: \n{individuals}")
    print(f"\nPR MATRIX: \n{pr}")

    names = individuals["Name"].tolist()
    males = individuals["Male"].tolist()
    females = individuals["Female"].tolist()
    allocate_first_group = individuals["AssignToFirstGroup"].tolist()
    alleles = individuals["Alleles"].tolist()
    priorities = individuals["Priority"].tolist()
    priority_values = individuals["PriorityValue"].tolist()
    connections = pr.values.tolist()

    # Check that no individuals in individuals.csv are being silently ignored for not being in the PR file.
    if len(connections) != len(individuals):
        raise app.UsageError("There is a mismatch between the number of individuals and the size of the PR matrix.")

    return (connections, group_defs, names, males, females, allocate_first_group, alleles, priorities, priority_values,
            objective_function, unique_id)


def solve_model(args):
    """
    Uses the processed data to define decision variables, constraints and the objective functions, solve the model and
    output final statistics.

    Args:
        :param args: Variable length argument list.
    """

    (connections, group_defs, names, males, females, allocate_first_group, alleles, priorities, priority_values,
     objective_function, unique_id) = build_data(args)

    num_individuals = len(connections)
    num_groups = len(group_defs)
    all_groups = range(num_groups)
    all_individuals = range(num_individuals)

    # Create the CP model.
    model = cp_model.CpModel()

    # ------------------------------ DECISION VARIABLES ------------------------------ #

    seats = {}
    individual_must_be_allocated = {}
    individual_allele_count = {}
    individual_priority_values = {}
    colocated = {}
    same_group = {}
    opposing_sex = {}

    for g in all_individuals:
        individual_must_be_allocated[g] = 1 if priorities[g] == 1 else 0
        individual_allele_count[g] = alleles[g]
        individual_priority_values[g] = priority_values[g]

    for t in all_groups:
        for g in all_individuals:
            seats[(t, g)] = model.NewBoolVar("individual %i placed in group %i" % (g, t))

    for g1 in range(num_individuals - 1):
        for g2 in range(g1 + 1, num_individuals):
            colocated[(g1, g2)] = model.NewBoolVar("individual %i placed with individual %i" % (g1, g2))

    for g1 in range(num_individuals - 1):
        for g2 in range(g1 + 1, num_individuals):
            for t in all_groups:
                same_group[(g1, g2, t)] = model.NewBoolVar(
                    "Individual %i paired with individual %i in group %i" % (g1, g2, t)
                )

    for g1 in range(num_individuals - 1):
        for g2 in range(g1 + 1, num_individuals):
            opposing_sex[(g1, g2)] = 0 if males[g1] == males[g2] else 1

    # ------------------------------ OBJECTIVE FUNCTIONS ------------------------------ #

    # Sum of pairwise relatedness across all pairs in the solution.
    sum_pairs_pr = sum(
        connections[g1][g2] * colocated[g1, g2]
        for g1 in range(num_individuals - 1)
        for g2 in range(g1 + 1, num_individuals)
        if connections[g1][g2] > 0
    )
    # Sum of ghost alleles across the solution.
    total_alleles = sum(
        seats[(t, g)] * individual_allele_count[g]
        for g in range(num_individuals)
        for t in range(num_groups)
    )
    # Sum of priority values across the solution.
    total_priority = sum(
        seats[(t, g)] * priority_values[g]
        for g in range(num_individuals)
        for t in range(num_groups)
    )

    # Attempt at weighted Chebyshev method
    ideal_total_pr = max([pr for col in connections for pr in col]) * num_groups  # All pairs have the best PR
    ideal_total_alleles = max(individual_allele_count) * 2 * num_groups  # All individuals have max alleles
    ideal_total_priority = 100 * 2 * num_groups  # All individuals have a priority value of 100
    deviation = model.NewIntVar(0, 999999, "max_deviation")

    if objective_function == GhostCobreederObjectiveFunction.MIN_PR:
        model.Maximize(sum_pairs_pr)
    elif objective_function == GhostCobreederObjectiveFunction.MAX_ALLELES:
        model.Maximize(total_alleles)
    elif objective_function == GhostCobreederObjectiveFunction.MAX_PRIO:
        model.Maximize(total_priority)
    elif objective_function == GhostCobreederObjectiveFunction.MIN_PR_MAX_ALLELES:
        model.Add(deviation >= args.weight_pr * (ideal_total_pr - sum_pairs_pr))
        model.Add(deviation >= args.weight_alleles * (ideal_total_alleles - total_alleles))
        model.Minimize(deviation)
    elif objective_function == GhostCobreederObjectiveFunction.MIN_PR_MAX_ALLELES_MAX_PRIO:
        model.Add(deviation >= args.weight_pr * (ideal_total_pr - sum_pairs_pr))
        model.Add(deviation >= args.weight_alleles * (ideal_total_alleles - total_alleles))
        model.Add(deviation >= args.weight_prio * (ideal_total_priority - total_priority))
        model.Minimize(deviation)

    # ------------------------------ CONSTRAINTS ------------------------------ #

    # Allocate at least args.total_individuals individuals.
    total_allocated = sum(seats[(t, g)] for g in range(num_individuals) for t in range(num_groups))
    model.Add(total_allocated >= args.total_individuals)

    for g in all_individuals:
        if individual_must_be_allocated[g]:
            # Priority individuals must be allocated to exactly one group.
            print("Individual %s must be allocated." % g)
            model.Add(sum(seats[(t, g)] for t in all_groups) == 1)
        else:
            # Non-priority individuals may or may not be allocated.
            model.Add(sum(seats[(t, g)] for t in all_groups) <= 1)
        if allocate_first_group[g] != -1:
            # Allocate first individual to the group specified.
            model.Add(seats[(allocate_first_group[g], g)] == 1)

    for t in all_groups:

        # Each group is filled to the required capacity with a fixed number of males and females.
        if group_defs['MinSize'][t] == group_defs['MaxSize'][t] == 2:  # Expected to be the default for coyotes
            # print("Group %s is a M-F pairing." % t)
            model.Add(sum(seats[(t, g)] for g in all_individuals) == 2)
            model.Add(sum(males[g] * seats[(t, g)] for g in all_individuals) == 1)
            model.Add(sum(females[g] * seats[(t, g)] for g in all_individuals) == 1)
        else:
            # print("Group %s has capacity of %i to %i." % (t, group_defs['MinSize'][t], group_defs['MaxSize'][t]))
            model.Add(sum(seats[(t, g)] for g in all_individuals) >= group_defs['MaxSize'][t])
            model.Add(sum(seats[(t, g)] for g in all_individuals) <= group_defs['MinSize'][t])
            model.Add(sum(males[g] * seats[(t, g)] for g in all_individuals) >= group_defs['NumMale'][t])
            model.Add(sum(females[g] * seats[(t, g)] for g in all_individuals) >= group_defs['NumFemale'][t])

        # Add PR constraint for groups whose PRThreshold value != -1.
        if group_defs['PRThreshold'][t] != -1:
            model.Add(
                sum(
                    connections[g1][g2] * same_group[(g1, g2, t)]
                    for g1 in range(num_individuals - 1)
                    for g2 in range(g1 + 1, num_individuals)
                    if connections[g1][g2] < group_defs['PRThreshold'][t]
                ) < 1
            )
        else:  # Set to global value if not specified for the specific group.
            model.Add(
                sum(
                    connections[g1][g2] * same_group[(g1, g2, t)]
                    for g1 in range(num_individuals - 1)
                    for g2 in range(g1 + 1, num_individuals)
                    if connections[g1][g2] < args.pr_threshold
                ) < 1
            )

    for g1 in range(num_individuals - 1):
        for g2 in range(g1 + 1, num_individuals):
            for t in all_groups:
                model.AddBoolOr(
                    [
                        seats[(t, g1)].Not(),
                        seats[(t, g2)].Not(),
                        same_group[(g1, g2, t)],
                    ]
                )
                model.AddImplication(same_group[(g1, g2, t)], seats[(t, g1)])
                model.AddImplication(same_group[(g1, g2, t)], seats[(t, g2)])

            model.Add(sum(same_group[(g1, g2, t)] for t in all_groups) == colocated[(g1, g2)])

    # Breaking symmetry. First individual is placed in the first group.
    print("\nInitial group allocations:")
    for g1 in range(len(allocate_first_group)):
        if allocate_first_group[g1] != -1:
            print("\tIndividual %i is allocated to group %i" % (g1, allocate_first_group[g1]))
            model.Add(seats[(allocate_first_group[g1], g1)] == 1)

    # ------------------------------ SOLVE MODEL ------------------------------ #

    solver = cp_model.CpSolver()
    paramstring = "%i,%i,%i" % (num_groups, objective_function, num_individuals)
    solution_printer = GhostCobreederPrinter(seats, names, num_groups, num_individuals, paramstring, unique_id)

    # solver.parameters.max_time_in_seconds = 5
    # solver.parameters.log_search_progress = True
    # solver.parameters.num_workers = 5
    # solver.parameters.fix_variables_to_their_hinted_value = True

    status = solver.Solve(model, solution_printer)

    # Print solution statistics
    print("\nCOBREEDER-COMPLETION [%s]" % args.unique_run_id)
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        print("\tStatistics - Optimal") if status == cp_model.OPTIMAL else print("Statistics - Feasible")
        print("\t\t- conflicts    : %i" % solver.NumConflicts())
        print("\t\t- branches     : %i" % solver.NumBranches())
        print("\t\t- wall time    : %f s" % solver.WallTime())
        print("\t\t- num solutions: %i" % solution_printer.num_solutions())
        # Save best solution to CSV if desired
        save = input("\nSave final solution to CSV? (Y/N): ")
        if save.lower() == 'y':
            best_solution = solution_printer.best_solution()
            save_solution_csv(args, connections, individual_allele_count, individual_priority_values, best_solution)
    else:
        print("No solution found.")


def main(argv: Sequence[str]) -> None:
    parser = argparse.ArgumentParser(prog='CoBreeder_for_GhostWolves',
                                     description='Group-Living Captive Breeding Solver.', add_help=True)
    subparsers = parser.add_subparsers(help='sub-command help')
    run_parser = subparsers.add_parser('run')
    run_parser.add_argument('individuals_file', type=str, help='CSV file detailing individuals.')
    run_parser.add_argument('pairwise_relatedness_file', type=str,
                            help='Scaled pairwise relatedness matrix.')
    run_parser.add_argument('group_file', type=str, help='CSV file detailing groups.')
    run_parser.add_argument("obj_function", type=str,
                            choices=[e.name for e in GhostCobreederObjectiveFunction],
                            help='String specifying objective function.')
    run_parser.add_argument("unique_run_id", type=str, help='Unique string identifier for each run.')
    run_parser.add_argument("weight_alleles", type=int, help='Weight for alleles.')
    run_parser.add_argument("weight_pr", type=int, help='Weight for pairwise relatedness.')
    run_parser.add_argument("weight_prio", type=int, help='Weight for priority values.')
    run_parser.add_argument("total_individuals", type=int,
                            help='Minimum number of individual to allocate to a solution.')
    run_parser.add_argument("pr_threshold", type=int, default=0,
                            help='Threshold for scaled PR permitted in a pairing.')
    run_parser.add_argument("exclude_disallow", type=str, choices=["EX", "ALL"],
                            help='Exclude individuals or specify disallowed pairings.')
    run_parser.add_argument("prio_calc_threshold", type=int, choices=range(0, 101),
                            help='Threshold for priority calculation. 0 to disable and use manual priority assignments '
                                 'only.')

    args = parser.parse_args()
    if (args.obj_function == "MIN_PR_MAX_ALLELES_MAX_PRIO") and (args.prio_calc_threshold == 0):
        print("Error: MIN_PR_MAX_ALLELES_MAX_PRIO requires priority calculations to be enabled.")
        sys.exit(1)

    solve_model(args)


if __name__ == "__main__":
    sys.setrecursionlimit(0x100000)
    app.run(main)
