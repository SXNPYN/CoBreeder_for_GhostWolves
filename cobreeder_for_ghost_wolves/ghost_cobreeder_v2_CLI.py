#!/usr/bin/env python3
from absl import app
import argparse
from enum import IntEnum
import numpy as np
from ortools.sat.python import cp_model
import os
import pandas as pd
import sys
import time
from typing import Sequence


class GhostCobreederObjectiveFunction(IntEnum):
    MIN_AV_PR = 1
    MAX_TOTAL_ALLELES = 2
    MAX_TOTAL_PRIO = 3
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
        num_groups, obj_f, num_ind = self.__paramstring.split(",")

        print(
            f"\nCOBREEDER-SOLUTION #{self.__solution_count}: time elapsed: {(current_time - self.__start_time):.5f}, "
            f"objective value: {int(objective)} ({GhostCobreederObjectiveFunction(int(obj_f)).name}), "
            f"experiment: {self.__uniqueid} ({num_ind} individuals, {num_groups} groups)"
        )

        grouped_ids = []

        for t in range(self.__num_groups):
            print("\tGroup %d: " % t)
            individuals = []
            for g in range(self.__num_individuals):
                if self.Value(self.__seats[(t, g)]):
                    print(f"\t\t{self.__names[g]}")
                    individuals.append((g, self.__names[g]))
                    grouped_ids.append(g)
            self.__best_solution[t] = individuals

        # Print ungrouped individuals
        all_ids = list(range(self.__num_individuals))
        ungrouped_ids = [i for i in all_ids if i not in grouped_ids]
        print("\n\tUngrouped individuals: ")
        for g in ungrouped_ids:
            print(f"\t\t{self.__names[g]}")

    def best_solution(self):
        return self.__best_solution

    def num_solutions(self):
        return self.__solution_count


def save_solution_csv(args, connections, individual_allele_count, individual_priority_value, best_solution):
    """
    Saves the best solution to a CSV file detailing groups, individuals, ghost alleles, and pairwise relatedness.

    Args:
        :param args: Variable length argument list.
        :param list of lists connections: PR relatedness matrix.
        :param dict individual_allele_count: Dictionary mapping every individual's ID to number of ghost alleles.
        :param dict individual_priority_value: Dictionary mapping every individual's ID to priority value.
        :param dict best_solution: Dictionary mapping every group number to (individual ID, individual name)
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
        solution_data.loc[i, "Priority_Sum"] = "N/A" if args.prio_calc_threshold == 0 else (
                individual_priority_value[individuals[0][0]] + individual_priority_value[individuals[1][0]])
        i += 1

    # Create results directory if it doesn't exist and save CSV
    results_dir_path = os.path.join(os.getcwd(), "results")
    os.makedirs(results_dir_path, exist_ok=True)
    out_file = f"best_solution_{args.unique_run_id}.csv"
    solution_data.to_csv(os.path.join(results_dir_path, out_file), index=False)
    print("Solution saved to results/%s." % out_file)


def calculate_priority(args, individuals, pr):
    """
    Dynamically calculates a priority value between 0 and 100 for each individual based on the dataset, where 100 is
    the highest priority possible.

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
        pr_threshold = args.global_pr_threshold
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
                                                ((b * female_individuals['NumMates']) / female_m_max)) * 100).astype(
            int)
        male_individuals['PriorityValue'] = (male_individuals['Proven'] *
                                             (((a * male_individuals['Alleles']) / male_g_max) +
                                              ((b * male_individuals['NumMates']) / male_m_max)) * 100).astype(int)

        individuals = pd.concat([male_individuals, female_individuals]).sort_index()
        individuals['NumMates'] = individuals['NumMates'].astype(int)

        individuals['Priority'] = 0  # 0 by default
        sorted_individuals = individuals.sort_values(by='PriorityValue', ascending=False)
        top_priority = sorted_individuals.head(args.prio_calc_threshold)
        individuals.loc[top_priority.index, 'Priority'] = 1  # Select the top x individuals to be "priority individuals"

    return individuals


def check_file_format(args):
    """
    Performs input validation for CSV files.

    :param args: Variable length argument list.
    :return pr: DataFrame detailing scaled PR for each individual.
    :return individuals: DataFrame detailing individuals.
    """

    try:
        # Files must be CSVs with no missing values
        pr = pd.read_csv(args.pairwise_relatedness_file, delimiter=',', header=None, skiprows=1)
        individuals = pd.read_csv(args.individuals_file, delimiter=',')
        if np.any(pr.isnull()) or np.any(individuals.isnull()):
            raise Exception("\nERROR: CSV files contain missing values.")

        # PR file must be symmetrical and only contain non-negative integers
        if not (pr >= 0).values.all():
            raise Exception("\nERROR: PR matrix contains negative values.")
        if not all(np.issubdtype(d, np.integer) for d in pr.dtypes):
            raise Exception("\nERROR: PR matrix contains non-integer values.")
        if not pr.equals(pr.T):
            raise Exception("\nERROR: PR matrix is not symmetrical.")

        # Individuals specification file must contain the required columns
        if list(individuals.columns) != ['Name', 'Male', 'Female', 'AssignToGroup', 'Alleles', 'Proven', 'Priority']:
            raise Exception("\nERROR: Unexpected column in individuals file.")
        # Alleles and AssignToGroup must be integers
        if not all(np.issubdtype(d, np.integer) for d in individuals[['AssignToGroup', 'Alleles']].dtypes):
            raise Exception("\nERROR: Non-integer value identified in individuals specification file.")
        # Alleles cannot be negative
        if not all(individuals['Alleles'] >= 0):
            raise Exception("\nERROR: Individual specification file contains negative alleles.")
        # AssignToGroup can only be -1 or a group ID
        if not ((individuals['AssignToGroup'] >= -1) & (individuals['AssignToGroup'] < args.num_pairs)).all():
            raise Exception("\nERROR: Invalid value in AssignToGroup column.")
        # Proven, Priority, Male, and Female can only take values 0 or 1
        for col in ['Proven', 'Priority', 'Male', 'Female']:
            if not set(individuals[col]) <= {0, 1}:
                raise Exception("\nERROR: Proven, Priority, Male, and Female can only be 0 or 1.")
        # Individuals can only be male or female
        if not ((individuals['Male'] + individuals['Female']) == 1).all():
            raise Exception("\nERROR: Individuals can only be male or female.")

    except Exception as e:
        print(e)
        sys.exit()

    return pr, individuals


def build_data(args):
    """
    Takes the CSV files (detailing individuals and pairwise relatedness) provided by the user and converts the
    data into a format that can be used by solve_model(). Also allows for exclusion of specific individuals and
    disallowed pairings if enabled via the command line.

    Args:
        :param args: Variable length argument list.
        :return: Various variables storing data from the CSV files, mostly columns in list format.
    """

    objective_function = GhostCobreederObjectiveFunction[args.obj_function]
    unique_id = args.unique_run_id
    pr, individuals = check_file_format(args)

    print("\nRUN ID: %s" % unique_id)
    print("OBJECTIVE FUNCTION: %i" % objective_function)
    print("NUMBER OF PAIRINGS TO ALLOCATE: %i" % args.num_pairs)
    print("USING PAIRWISE RELATEDNESS FILE [%s]" % args.pairwise_relatedness_file)

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
                    print("Please ensure indices are valid (e.g. 0,6,17,3)")
            else:
                break

    individuals = calculate_priority(args, individuals, pr)
    print(f"\nSUMMARY OF INDIVIDUALS AFTER PROCESSING: \n{individuals}")
    print(f"\nPR MATRIX: \n{pr}")

    # Each group has default PR threshold unless specified
    group_prs = {key: -1 for key in range(args.num_pairs)}

    # Specify custom PR thresholds for groups
    if args.specify_pr == "CUSTOM_PR":

        print(f"\nGROUP IDs: {list(group_prs.keys())}")
        print("NOTE: You do not need to specify PR thresholds for all groups. Unspecified groups will take the default"
              f" value ({args.global_pr_threshold}).")
        print("Remember that PR values are scaled.")

        while True:
            custom_prs = input("Specify custom PR thresholds (List of group-PR pairs e.g. 3-50, 0-100): ")
            if custom_prs:
                try:
                    custom_prs = [tuple(map(int, x.split('-'))) for x in custom_prs.strip().split(",")]
                    # Set custom PRs in dictionary
                    for i, j in custom_prs:
                        if i not in list(range(args.num_pairs)):  # Ensure group ID is valid
                            raise Exception
                        group_prs[i] = j
                    break
                except (IndexError, ValueError, Exception):
                    print("\nInvalid input. Please enter valid group IDs and scaled PR values.")
            else:
                break

    names = individuals["Name"].tolist()
    males = individuals["Male"].tolist()
    females = individuals["Female"].tolist()
    group_prs = list(group_prs.values())
    allocate_first_group = individuals["AssignToGroup"].tolist()
    alleles = individuals["Alleles"].tolist()
    priorities = individuals["Priority"].tolist()
    priority_values = individuals["PriorityValue"].tolist()
    connections = pr.values.tolist()

    # Check that no individuals in individuals.csv are being silently ignored for not being in the PR file.
    if len(connections) != len(individuals):
        raise app.UsageError("There is a mismatch between the number of individuals and the size of the PR matrix.")

    return (connections, names, males, females, group_prs, allocate_first_group, alleles, priorities, priority_values,
            objective_function, unique_id)


def solve_model(args):
    """
    Uses the processed data to define decision variables, constraints and the objective functions, solve the model and
    output final statistics.

    Args:
        :param args: Variable length argument list.
    """

    (connections, names, males, females, group_prs, allocate_first_group, alleles, priorities, priority_values,
     objective_function, unique_id) = build_data(args)

    num_individuals = len(connections)
    num_groups = args.num_pairs
    all_groups = range(num_groups)
    all_individuals = range(num_individuals)

    # Create the CP model.
    model = cp_model.CpModel()

    # --------------------------------------------- DECISION VARIABLES --------------------------------------------- #

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

    # -------------------------------------------- OBJECTIVE FUNCTIONS -------------------------------------------- #

    # Sum of pairwise relatedness across all pairs in the solution.
    total_pr = sum(
        connections[g1][g2] * colocated[g1, g2]
        for g1 in range(num_individuals - 1)
        for g2 in range(g1 + 1, num_individuals)
        if connections[g1][g2] > 0
    )
    # Average PR across all pairs in the solution.
    av_pair_pr = model.NewIntVar(0, 100000000, "av_pair_pr")
    model.Add(total_pr == av_pair_pr * num_groups)  # av_pair_pr = total_pr / num_groups

    # Sum of ghost alleles across the solution.
    total_alleles = sum(
        seats[(t, g)] * individual_allele_count[g]
        for g in range(num_individuals)
        for t in range(num_groups)
    )
    # Average number of alleles across all pairs in the solution.
    av_pair_alleles = model.NewIntVar(0, 100000000, "av_pair_alleles")
    model.Add(total_alleles == av_pair_alleles * num_groups)

    # Sum of priority values across the solution.
    total_priority = sum(
        seats[(t, g)] * priority_values[g]
        for g in range(num_individuals)
        for t in range(num_groups)
    )
    # Average sum of priority values across all pairs in the solution.
    av_pair_priority = model.NewIntVar(0, 100000000, "av_pair_priority")
    model.Add(total_priority == av_pair_priority * num_groups)

    # Single-objective optimisation
    if objective_function == GhostCobreederObjectiveFunction.MIN_AV_PR:
        model.Maximize(av_pair_pr)
    elif objective_function == GhostCobreederObjectiveFunction.MAX_TOTAL_ALLELES:
        model.Maximize(total_alleles)
    elif objective_function == GhostCobreederObjectiveFunction.MAX_TOTAL_PRIO:
        model.Maximize(total_priority)

    # Multi-objective optimisation
    elif (objective_function == GhostCobreederObjectiveFunction.MIN_PR_MAX_ALLELES) or \
         (objective_function == GhostCobreederObjectiveFunction.MIN_PR_MAX_ALLELES_MAX_PRIO):

        # Calculate PR and sum of alleles/priority values for an ideal pair
        ideal_pair_pr = max([pr for col in connections for pr in col])  # Uses best PR in matrix as ideal value.
        ideal_pair_alleles = max(individual_allele_count.values()) * 2  # Uses best alleles as ideal value.
        ideal_pair_priority = 100 * 2  # Max priority for an individual is 100, hence max across a pair is 200.

        # PR
        scaled_pr_difference = model.NewIntVar(0, 100000000, "scaled_pr_difference")
        model.Add(scaled_pr_difference == (100 * ideal_pair_pr) - (100 * av_pair_pr))
        percent_deviation_pr = model.NewIntVar(0, 100000000, "percent_deviation_pr")
        model.Add(scaled_pr_difference == percent_deviation_pr * ideal_pair_pr)
        weighted_pr = model.NewIntVar(0, 100000000, "weighted_pr")
        model.AddMultiplicationEquality(weighted_pr, [percent_deviation_pr, args.weight_pr])

        # Alleles
        scaled_allele_difference = model.NewIntVar(0, 100000000, "scaled_allele_difference")
        model.Add(scaled_allele_difference == (100 * ideal_pair_alleles) - (100 * av_pair_alleles))
        percent_deviation_alleles = model.NewIntVar(0, 100000000, "percent_deviation_alleles")
        model.Add(scaled_allele_difference == percent_deviation_alleles * ideal_pair_alleles)
        weighted_alleles = model.NewIntVar(0, 100000000, "weighted_alleles")
        model.AddMultiplicationEquality(weighted_alleles, [percent_deviation_alleles, args.weight_alleles])

        if objective_function == GhostCobreederObjectiveFunction.MIN_PR_MAX_ALLELES:
            # deviation = model.NewIntVar(0, 1000000000, "max_deviation")  # Value to minimise
            # model.Add(deviation >= weighted_pr)
            # model.Add(deviation >= weighted_alleles)
            # model.Minimize(deviation)

            combined_pr_alleles = model.NewIntVar(0, 1000000000, "combined_pr_alleles")
            model.Add(combined_pr_alleles == weighted_pr + weighted_alleles)
            model.Minimize(combined_pr_alleles)

        elif objective_function == GhostCobreederObjectiveFunction.MIN_PR_MAX_ALLELES_MAX_PRIO:

            # Priority
            scaled_priority_difference = model.NewIntVar(0, 100000000, "scaled_priority_difference")
            model.Add(scaled_priority_difference == (100 * ideal_pair_priority) - (100 * av_pair_priority))
            percent_deviation_priority = model.NewIntVar(0, 100000000, "percent_deviation_priority")
            model.Add(scaled_priority_difference == percent_deviation_priority * ideal_pair_priority)
            weighted_priority = model.NewIntVar(0, 100000000, "weighted_priority")
            model.AddMultiplicationEquality(weighted_priority, [percent_deviation_priority,
                                                                args.weight_prio])

            combined_pr_alleles_prio = model.NewIntVar(0, 1000000000, "combined_pr_alleles_prio")
            model.Add(combined_pr_alleles_prio == weighted_pr + weighted_alleles + weighted_priority)
            model.Minimize(combined_pr_alleles_prio)

    # ----------------------------------------------- CONSTRAINTS ----------------------------------------------- #

    for g in all_individuals:
        if individual_must_be_allocated[g]:
            # Priority individuals must be allocated to exactly one group.
            print("%s (ID: %s) must be allocated." % (names[g], g))
            model.Add(sum(seats[(t, g)] for t in all_groups) == 1)
        else:
            # Non-priority individuals may or may not be allocated.
            model.Add(sum(seats[(t, g)] for t in all_groups) <= 1)

    for t in all_groups:
        # Each pair is filled to the required capacity with one male and one female.
        model.Add(sum(seats[(t, g)] for g in all_individuals) == 2)
        model.Add(sum(males[g] * seats[(t, g)] for g in all_individuals) == 1)
        model.Add(sum(females[g] * seats[(t, g)] for g in all_individuals) == 1)

    for group, pr in enumerate(group_prs):
        print(f"Group:{group}, PR threshold:{pr if pr != -1 else args.global_pr_threshold}")
        if pr != -1:  # Set custom PR constraint
            model.Add(
                sum(
                    connections[g1][g2] * same_group[(g1, g2, group)]
                    for g1 in range(num_individuals - 1)
                    for g2 in range(g1 + 1, num_individuals)
                    if connections[g1][g2] < pr
                ) < 1
            )
        else:  # Set to global value if no custom PR is specified
            model.Add(
                sum(
                    connections[g1][g2] * same_group[(g1, g2, group)]
                    for g1 in range(num_individuals - 1)
                    for g2 in range(g1 + 1, num_individuals)
                    if connections[g1][g2] < args.global_pr_threshold
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
    for g in range(len(allocate_first_group)):
        if allocate_first_group[g] != -1:
            print("\t%s (ID: %s) is allocated to group %i" % (names[g], g, allocate_first_group[g]))
            model.Add(seats[(allocate_first_group[g], g)] == 1)

    # --------------------------------------------- SOLVE MODEL --------------------------------------------- #

    solver = cp_model.CpSolver()
    paramstring = "%i,%i,%i" % (num_groups, objective_function, num_individuals)
    solution_printer = GhostCobreederPrinter(seats, names, num_groups, num_individuals, paramstring, unique_id)

    solver.parameters.max_time_in_seconds = 1800
    # solver.parameters.log_search_progress = True
    # solver.parameters.num_workers = 1
    # solver.parameters.fix_variables_to_their_hinted_value = True

    status = solver.Solve(model, solution_printer)

    # Print solution statistics
    print("\nCOBREEDER-COMPLETION [%s]" % args.unique_run_id)
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        print("\tStatistics")
        print("\t\t- Status       : %s" % solver.status_name(status))
        print("\t\t- Conflicts    : %i" % solver.NumConflicts())
        print("\t\t- Branches     : %i" % solver.NumBranches())
        print("\t\t- Wall time    : %f s" % solver.WallTime())
        print("\t\t- Num solutions: %i" % solution_printer.num_solutions())
        # Save best solution to CSV if desired
        save = input("\nSave final solution to CSV? (Y/N): ")
        if save.lower() == 'y':
            best_solution = solution_printer.best_solution()
            save_solution_csv(args, connections, individual_allele_count, individual_priority_values, best_solution)
    else:
        print("No solution found.")
        print("Status: %s" % solver.status_name(status))


def main(argv: Sequence[str]) -> None:
    parser = argparse.ArgumentParser(prog='CoBreeder_for_GhostWolves',
                                     description='Group-Living Captive Breeding Solver.', add_help=True)
    subparsers = parser.add_subparsers(help='sub-command help')
    run_parser = subparsers.add_parser('run')
    run_parser.add_argument('individuals_file', type=str, help='CSV file detailing individuals.')
    run_parser.add_argument('pairwise_relatedness_file', type=str,
                            help='Scaled pairwise relatedness matrix.')
    run_parser.add_argument('num_pairs', type=int, help='Number of pairings to allocate.')
    run_parser.add_argument("specify_pr", type=str, choices=["CUSTOM_PR", "DEFAULT_PR"],
                            help='Specify custom PR for specific groups.')
    run_parser.add_argument("obj_function", type=str,
                            choices=[e.name for e in GhostCobreederObjectiveFunction],
                            help='String specifying objective function.')
    run_parser.add_argument("unique_run_id", type=str, help='Unique string identifier for each run.')
    run_parser.add_argument("weight_alleles", type=int, help='Weight for alleles.')
    run_parser.add_argument("weight_pr", type=int, help='Weight for pairwise relatedness.')
    run_parser.add_argument("weight_prio", type=int, help='Weight for priority values.')
    run_parser.add_argument("global_pr_threshold", type=int, default=0,
                            help='Threshold for scaled PR permitted in a pairing.')
    run_parser.add_argument("exclude_disallow", type=str, choices=["EX", "ALL"],
                            help='Exclude individuals or specify disallowed pairings.')
    run_parser.add_argument("prio_calc_threshold", type=int, choices=range(0, 101),
                            help='Threshold for priority calculations representing the number of individuals that can'
                                 'fall into the priority set. 0 to disable and use manual priority assignments only.')

    args = parser.parse_args()

    if ((args.obj_function == "MIN_PR_MAX_ALLELES_MAX_PRIO" or args.obj_function == "MAX_PRIO") and
            (args.prio_calc_threshold == 0)):
        print("Error: MIN_PR_MAX_ALLELES_MAX_PRIO requires priority calculations to be enabled.")
        sys.exit(1)

    if args.prio_calc_threshold > args.num_pairs:
        print("Error: The size of the priority set cannot be larger than the number of pairs.")
        sys.exit(1)

    solve_model(args)


if __name__ == "__main__":
    sys.setrecursionlimit(0x100000)
    app.run(main)
