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

PR_THRESHOLD = 0  # TODO
MAX_TIME_SECONDS = 5  # TODO Threshold for max time allowed


class CobreederObjectiveFunction(IntEnum):
    MIN_PR = 1
    MAX_ALLELES = 2
    MAX_PRIO = 3
    MIN_PR_MAX_ALLELES = 4
    MIN_PR_MAX_ALLELES_MAX_PRIO = 5


class CobreederPrinter(cp_model.CpSolverSolutionCallback):
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


def save_solution_csv(args, connections, individual_allele_count, best_solution):
    """
    # TODO Add docstrings
    """
    solution_data = pd.DataFrame(columns=["Group", "Ind_1_Name", "Ind_2_Name", "Ind_1_ID", "Ind_2_ID",
                                          "Ind_1_Alleles", "Ind_2_Alleles", "Pairwise_Relatedness"])
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
        i += 1

    # Create results directory if it doesn't exist and save CSV
    results_dir_path = os.path.join(os.getcwd(), "results")
    os.makedirs(results_dir_path, exist_ok=True)
    out_file = f"best_solution_{args.unique_run_id}.csv"
    solution_data.to_csv(os.path.join(results_dir_path, out_file), index=False)
    print("Solution saved to results/%s." % out_file)


def calculate_priority(individuals, prio_threshold, pr):
    """
    # TODO Add docstrings
    """

    if prio_threshold == 0:  # Use values from csv only
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
        female_g_max = max(female_individuals["Alleles"])  # Highest/best number of ghost alleles amongst females
        male_g_max = max(male_individuals["Alleles"])  # Highest/best number of ghost alleles amongst males
        female_m_max = len(male_individuals)  # Max number of potential mates that a female can have
        male_m_max = len(female_individuals)  # Max number of potential mates that a male can have

        # Calculate number of potential mates each individual has from the PR matrix
        for i in female_individuals.index:
            num_mates = (pr.loc[male_individuals.index][i] > PR_THRESHOLD).sum()
            female_individuals.loc[i, "NumMates"] = num_mates
            print(f'Female {i} has {num_mates} potential mates.')
        for i in male_individuals.index:
            num_mates = (pr.loc[female_individuals.index][i] > PR_THRESHOLD).sum()
            male_individuals.loc[i, "NumMates"] = num_mates
            print(f'Male {i} has {num_mates} mates.')

        female_individuals['PriorityValue'] = (female_individuals['Proven'] *
                                               (((a * female_individuals['Alleles']) / female_g_max) +
                                               ((b * female_individuals['NumMates']) / female_m_max)) * 100).astype(int)

        male_individuals['PriorityValue'] = (male_individuals['Proven'] *
                                             (((a * male_individuals['Alleles']) / male_g_max) +
                                             ((b * male_individuals['NumMates']) / male_m_max)) * 100).astype(int)

        individuals = pd.concat([male_individuals, female_individuals]).sort_index()
        individuals['NumMates'] = individuals['NumMates'].astype(int)
        individuals['Priority'] = [1 if p > prio_threshold else 0 for p in individuals['PriorityValue']]

    return individuals


def build_data(args):
    """
    # TODO Add docstrings
    """

    objective_function = CobreederObjectiveFunction[args.obj_function]
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

    individuals = calculate_priority(individuals, args.prio_calc_threshold, pr)
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


def solve_with_discrete_model(args):
    """
    # TODO Add docstrings
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
    colocated = {}
    same_group = {}
    opposing_sex = {}

    for g in all_individuals:
        individual_must_be_allocated[g] = 1 if priorities[g] == 1 else 0
        individual_allele_count[g] = alleles[g]

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
    total_priorities = sum(
        seats[(t, g)] * priority_values[g]
        for g in range(num_individuals)
        for t in range(num_groups)
    )

    w_pr = args.weight_pr
    w_alleles = args.weight_alleles
    w_prio = 0.2*w_pr
    pr_alleles = w_pr * sum_pairs_pr + w_alleles * total_alleles + w_prio * total_priorities

    if objective_function == CobreederObjectiveFunction.MIN_PR:
        model.Maximize(sum_pairs_pr)
    elif objective_function == CobreederObjectiveFunction.MAX_ALLELES:
        model.Maximize(total_alleles)
    elif objective_function == CobreederObjectiveFunction.MAX_PRIO:
        model.Maximize(total_priorities)
    elif objective_function == CobreederObjectiveFunction.MIN_PR_MAX_ALLELES:
        model.Maximize(total_priorities)  # TODO Maximise ghost while minimising pr
    elif objective_function == CobreederObjectiveFunction.MIN_PR_MAX_ALLELES_MAX_PRIO:
        model.Maximize(pr_alleles)  # TODO Maximise ghost while minimising pr whilst maximising priority values?

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
                    if connections[g1][g2] < PR_THRESHOLD
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
    solution_printer = CobreederPrinter(seats, names, num_groups, num_individuals, paramstring, unique_id)

    solver.parameters.max_time_in_seconds = MAX_TIME_SECONDS
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
            save_solution_csv(args, connections, individual_allele_count, best_solution)
    else:
        print("No solution found.")


def main(argv: Sequence[str]) -> None:
    parser = argparse.ArgumentParser(prog='CoBreeder_for_GhostWolves',
                                     description='Group-Living Captive Breeding Solver.', add_help=True)
    subparsers = parser.add_subparsers(help='sub-command help')
    run_parser = subparsers.add_parser('run')
    run_parser.add_argument('individuals_file', type=str,
                            help='CSV file detailing individuals.')
    run_parser.add_argument('pairwise_relatedness_file', type=str,
                            help='Scaled pairwise relatedness matrix.')
    run_parser.add_argument('group_file', type=str,
                            help='CSV file detailing groups.')
    run_parser.add_argument("obj_function", type=str,
                            choices=[e.name for e in CobreederObjectiveFunction],
                            help='String specifying objective function.')
    run_parser.add_argument("unique_run_id", type=str,
                            help='Unique string identifier for each run.')
    run_parser.add_argument("weight_alleles", type=int,
                            help='Weight for alleles.')
    run_parser.add_argument("weight_pr", type=int,
                            help='Weight for PR.')
    run_parser.add_argument("total_individuals", type=int,
                            help='Minimum number of individuals allocated to a solution.')
    run_parser.add_argument("exclude_disallow", type=str, choices=["EX", "ALL"],
                            help='Exclude individuals or specify disallowed pairings with "EX", or use all with "ALL".')
    run_parser.add_argument("prio_calc_threshold", type=int, choices=range(0, 101),
                            help='Threshold for priority calculation. 0 to disable and use manual priority assignments '
                                 'only.')
    run_parser.add_argument('subst', nargs='?', default=0, type=int)

    args = parser.parse_args()
    solve_with_discrete_model(args)


if __name__ == "__main__":
    sys.setrecursionlimit(0x100000)
    app.run(main)
