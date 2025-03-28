#!/usr/bin/env python3
from absl import app
import argparse
from enum import IntEnum
from ortools.sat.python import cp_model
import pandas as pd
import sys
import time
from typing import Sequence

PR_THRESHOLD = 0 # TODO Threshold for min pairwise relatedness permitted (value not included)
MAX_TIME_SECONDS = -1

#TODO
class CobreederObjectiveFunction(IntEnum):
    ALL_PAIRS = 11
    MALE_FEMALE = 2
    ALLELES_MIN = 3
    ALLELES_MAX = 4
    ALL_PAIRS_PR_MIN = 5
    ALL_PAIRS_PR_MAX = 6
    WEIGHTED_ALLELES_PR_50_50 = 7
    ALL_PAIRS_PR_MIN_SQUARED = 8
    ALL_PAIRS_PR_MAX_SQUARED = 9
    MALE_FEMALE_SQUARED = 10
    SWINGER = 11


class CobreederPrinter(cp_model.CpSolverSolutionCallback):
    """
    # TODO Add docstrings
    """

    def __init__(self, seats, names, num_groups, num_individuals, paramstring, unique_id):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__solution_count = 0
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
            for g in range(self.__num_individuals):
                if self.Value(self.__seats[(t, g)]):
                    print(f"\t\t{self.__names[g]}" )

    def num_solutions(self):
        return self.__solution_count


def save_solution_csv(args):
    """
    # TODO Add docstrings
    """
    out_file = f"best_allocation_{args.unique_run_id}.csv"

    df = pd.DataFrame(columns=["Group", "Ind_M", "Ind_F","Alleles_M", "Alleles_F", "PairwiseRelatedness"])
    # TODO Write data

    df.to_csv(out_file, index=False)
    print("Solution saved to %s" % out_file)


def calculate_priority(individuals, prio_threshold, pr):
    """
    # TODO Add docstrings
    """

    if prio_threshold == 0: # Use values from csv only
        priorities = individuals["Priority"].tolist()
        individuals['PriorityValue'] = [100 * p for p in priorities] # All priority individuals are equal

    else: # Dynamically calculate priority                            ... but some are more equal than others.
        while True:
            a = input("Weight to place on ghost alleles: ")
            try:
                a = float(a)
                if 0 <= a <= 1: break
            except ValueError:
                print("Please enter a number between 0.0 and 1.0")
        b = 1.0 - a
        #print(f"\nPlacing weight on ghost alleles and number of mates in a ratio of {int(10*a)}:{int(10*b)}.")
        print("\nPlacing weight on ghost alleles and number of mates in a ratio of %i:%i" % (10*a, 10*b))

        female_individuals = individuals.query("Female == 1").copy()
        male_individuals = individuals.query("Male == 1").copy()
        female_g_max = max(female_individuals["Alleles"]) # Highest number of ghost alleles amongst females
        male_g_max = max(male_individuals["Alleles"]) # Highest number of ghost alleles amongst males
        female_m_max = len(male_individuals) # Max number of potential mates that a female can have
        male_m_max = len(female_individuals) # Max number of potential mates that a male can have

        # Calculate number of potential mates each individual has from the PR matrix
        pr_f = pr.loc[female_individuals.index]
        pr_m = pr.loc[male_individuals.index]

        for i in female_individuals.index:
            num_mates = (pr_m[i] > PR_THRESHOLD).sum()
            female_individuals.loc[i, "NumMates"] = num_mates
            #print(f'Female {i} has {num_mates} potential mates.')
        for i in male_individuals.index:
            num_mates = (pr_f[i] > PR_THRESHOLD).sum()
            male_individuals.loc[i, "NumMates"] = num_mates
            #print(f'Male {i} has {num_mates} mates.')

        female_individuals['PriorityValue'] = (female_individuals['Proven'] *
                                          (((a * female_individuals['Alleles']) / female_g_max)
                                           + ((b * female_individuals['NumMates']) / male_m_max)) * 100).astype(int)
        female_individuals['Priority'] = (female_individuals['PriorityValue'] > prio_threshold).astype(int)

        male_individuals['PriorityValue'] = (male_individuals['Proven'] *
                                               (((a * male_individuals['Alleles']) / male_g_max)
                                                + ((b * male_individuals['NumMates']) / female_m_max)) * 100).astype(int)
        male_individuals['Priority'] = (male_individuals['PriorityValue'] > prio_threshold).astype(int)

        priorities = pd.concat([male_individuals['Priority'],
                                female_individuals['Priority']]).sort_index()
        priority_values = pd.concat([male_individuals['PriorityValue'],
                                     female_individuals['PriorityValue']]).sort_index()

        individuals = pd.concat([male_individuals, female_individuals]).sort_index()
        individuals['NumMates'] = individuals['NumMates'].apply(lambda x: int(x))

    return individuals


def build_data(args):
    """
    # TODO Add docstrings
    """

    objective_function = CobreederObjectiveFunction[args.obj_function]
    unique_id = args.unique_run_id
    print("\nRUN ID: %s \nOBJECTIVE FUNCTION: %i" % (args.unique_run_id, objective_function))

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
                    disallowed_pairings = [tuple(map(int, x.split('-'))) for x in disallowed_pairings.strip().split(",")]
                    # Set PR for disallowed combinations to 0
                    for i, j in disallowed_pairings:
                        pr.iloc[i, j] = 0
                        pr.iloc[j, i] = 0
                    break
                except (IndexError, ValueError):
                    print("Invalid input. Please only include indices matching to the individuals above (e.g. 0-1, 4-2)")
            else: break


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
            else: break

    individuals = calculate_priority(individuals, args.prio_calc_threshold, pr)
    print(f"\nSUMMARY OF INDIVIDUALS AFTER PROCESSING: \n{individuals}")
    print(f"\nPR MATRIX: \n{pr}")

    names = individuals["Name"].tolist()
    males = individuals["Male"].tolist()
    females = individuals["Female"].tolist()
    allocate_first_group = individuals["AssignToFirstCorral"].tolist()
    species = individuals["Species"].tolist()
    alleles = individuals["Alleles"].tolist()
    priorities = individuals["Priority"].tolist()
    priority_values = individuals["PriorityValue"].tolist()

    # Check that no individuals in individuals.csv are being silently ignored for not being in the PR file.
    connections = pr.values.tolist()
    if len(connections) != len(individuals):
        raise app.UsageError("There is a mismatch between the number of individuals and the size of the PR matrix.")

    return (connections, group_defs, names, males, females, allocate_first_group, species, alleles, priorities,
            priority_values, objective_function, unique_id)

# TODO
def solve_with_discrete_model(args):
    """
    # TODO Add docstrings
    """

    (connections, group_defs, names, males, females, allocate_first_group, species, alleles, priorities,
     priority_values, objective_function, unique_id) = build_data(args)

    num_individuals = len(connections)
    num_groups = len(group_defs)
    all_groups = range(num_groups)
    all_individuals = range(num_individuals)

    # Create the CP model.
    model = cp_model.CpModel()



    x = input("\n\n=^..^=   =^..^=   =^..^=    =^..^=    =^..^=    =^..^=    =^..^=    =^..^=    =^..^=    =^..^=")

    # TODO priority ranking

    #
    # Decision variables
    #
    seats = {}
    individual_corral_compatibility = {}
    optional_group_allocation = {}
    compulsory_match_individual_corral_violated = {}
    individual_must_be_allocated = {}
    individual_allele_count = {}
    print("Corral, Individual, CompGroup, OptionalGroup")

    for g in all_individuals:
        individual_must_be_allocated[g] = 1 if species[g] != "R" else 0
        individual_allele_count[g] = alleles[g]

    for t in all_groups:
        for g in all_individuals:
            seats[(t, g)] = model.NewBoolVar("individual %i placed in corral %i" % (g, t))
            print("%s %s %s %s" % (t, g, species[g] == group_defs['CompGroup'][t],
                                   species[g] == group_defs['OptionalGroup'][t]))
            individual_corral_compatibility[(t, g)] = 1 if species[g] == group_defs['CompGroup'][t] or \
                                                           species[g] == group_defs['OptionalGroup'][t] else 0
            optional_group_allocation[(t, g)] = 1 if species[g] == group_defs['OptionalGroup'][t] else 0
            compulsory_match_individual_corral_violated[(t, g)] = 1 if species[g] != group_defs['CompGroup'][t] and \
                                                                       species[g] != "R" else 0

    colocated = {}
    for g1 in range(num_individuals - 1):
        for g2 in range(g1 + 1, num_individuals):
            colocated[(g1, g2)] = model.NewBoolVar(
                "individual %i colocated with individual %i" % (g1, g2)
            )

    same_corral = {}
    for g1 in range(num_individuals - 1):
        for g2 in range(g1 + 1, num_individuals):
            for t in all_groups:
                same_corral[(g1, g2, t)] = model.NewBoolVar(
                    "guest %i seats with guest %i on corral %i" % (g1, g2, t)
                )

    opposing_sex = {}
    for g1 in range(num_individuals - 1):
        for g2 in range(g1 + 1, num_individuals):
            opposing_sex[(g1, g2)] = 0 if males[g1] == males[g2] else 1

    # TODO implement strategy pattern for CobreederObjectiveFunction

    alleles = sum(
        seats[(t, g)] * individual_allele_count[g]
        for g in range(num_individuals)
        for t in range(num_groups)
    )
    all_pairs_pr = sum(
        connections[g1][g2] * colocated[g1, g2]
        for g1 in range(num_individuals - 1)
        for g2 in range(g1 + 1, num_individuals)
        if connections[g1][g2] > 0
    )
    opposing_sex_pr = sum(
        connections[g1][g2] * colocated[g1, g2] * opposing_sex[g1, g2]
        for g1 in range(num_individuals - 1)
        for g2 in range(g1 + 1, num_individuals)
        if connections[g1][g2] > 0
    )

    all_pairs_pr_squared = sum(
        connections[g1][g2] * connections[g1][g2] * colocated[g1, g2]  # * colocated[g1, g2]
        for g1 in range(num_individuals - 1)
        for g2 in range(g1 + 1, num_individuals)
        if connections[g1][g2] > 0
    )

    swingerA = sum(
        connections[g1][g2] * connections[g1][g2] * colocated[g1, g2]  # * colocated[g1, g2]
        for g1 in range(num_individuals - 1)
        for g2 in range(g1 + 1, num_individuals)
        if connections[g1][g2] > 0
    )
    swingerB = sum(
        -1 * colocated[g1, g2]  # * colocated[g1, g2]
        for g1 in range(num_individuals - 1)
        for g2 in range(g1 + 1, num_individuals)
        if connections[g1][g2] > 0
    )

    opposing_sex_pr_squared = sum(
        connections[g1][g2] * connections[g1][g2] * colocated[g1, g2] * opposing_sex[g1, g2]
        for g1 in range(num_individuals - 1)
        for g2 in range(g1 + 1, num_individuals)
        if connections[g1][g2] > 0
    )

    if objective_function == CobreederObjectiveFunction.ALL_PAIRS:
        model.Minimize(all_pairs_pr)
    elif objective_function == CobreederObjectiveFunction.MALE_FEMALE:
        model.Minimize(opposing_sex_pr)
    elif objective_function == CobreederObjectiveFunction.MALE_FEMALE_SQUARED:
        model.Minimize(opposing_sex_pr_squared)
    elif objective_function == CobreederObjectiveFunction.ALLELES_MIN:
        model.Minimize(alleles)
    elif objective_function == CobreederObjectiveFunction.ALLELES_MAX:
        model.Maximize(alleles)
    elif objective_function == CobreederObjectiveFunction.ALL_PAIRS_PR_MIN:
        model.Minimize(all_pairs_pr)
    elif objective_function == CobreederObjectiveFunction.ALL_PAIRS_PR_MIN_SQUARED:
        model.Minimize(all_pairs_pr_squared)
    elif objective_function == CobreederObjectiveFunction.ALL_PAIRS_PR_MAX:
        model.Maximize(all_pairs_pr)
    elif objective_function == CobreederObjectiveFunction.ALL_PAIRS_PR_MAX_SQUARED:
        model.Maximize(all_pairs_pr_squared)
    elif objective_function == CobreederObjectiveFunction.SWINGER:
        # model.Minimize(swingerA * (-1 * swingerB))
        model.Minimize(swingerA + swingerB)
    elif objective_function == CobreederObjectiveFunction.WEIGHTED_ALLELES_PR_50_50:
        # Hard coded values to be abstracted for camera ready artifact for AAAI25
        weighted_alleles = args.weight_alleles * ((alleles - 2721) * (2721 - 4069))
        weighted_pr = args.weight_pr * ((all_pairs_pr_squared - 7171) * (7171 - 640474))
        print("weighted_alleles = %s, weighted_pr = %s" % (weighted_alleles, weighted_pr))
        model.Maximize(weighted_pr + (-1 * weighted_alleles))

    # Constraints

    total_allocated = sum(
        seats[(t, g)]
        for g in range(num_individuals)
        for t in range(num_groups)
    )
    model.Add(total_allocated >= args.total_individuals)

    # Setting individual-centric constraints.
    for g in all_individuals:

        # All individuals which must be placed are allocated to a corral.
        # Else, each individual is placed in most one corral.
        if individual_must_be_allocated[g]:
            print("%s must be alloc" % g)
            model.Add(sum(seats[(t, g)] for t in all_groups) == 1)
        else:
            model.Add(sum(seats[(t, g)] for t in all_groups) <= 1)

        # model.Add(sum(seats[(t, g)] for t in all_groups) == 1)

        # All individuals which are pre-allocated to a corral, are placed accordingly.
        if allocate_first_group[g] != -1:
            model.Add(seats[(allocate_first_group[g], g)] == 1)

    # Setting corral-centric constraints
    for t in all_groups:
        print("Corral %s is of size %i to %i; compulsory %s and optional %s" % (
            t, group_defs['MinSize'][t], group_defs['MaxSize'][t], group_defs['CompGroup'][t],
            group_defs['OptionalGroup'][t],))

        # Each corral is filled to the required capacity.
        model.Add(sum(seats[(t, g)] for g in all_individuals) >= group_defs['MaxSize'][t])
        model.Add(sum(seats[(t, g)] for g in all_individuals) <= group_defs['MinSize'][t])

        # Each corral is allocated a fixed number of male individuals.
        model.Add(sum(males[g] * seats[(t, g)] for g in all_individuals) >= group_defs['NumMale'][t])

        # Each corral is allocated a fixed number of female individuals.
        model.Add(sum(females[g] * seats[(t, g)] for g in all_individuals) >= group_defs['NumFemale'][t])

        print(sum(individual_corral_compatibility[(t, g)] * seats[(t, g)] for g in all_individuals))

        # Cap the number of 'optional' group allocations
        if group_defs['MaxNumNonComp'][t] != -1:
            print("circulate")
            print(sum(optional_group_allocation[(t, g)] * seats[(t, g)] for g in all_individuals))
            # model.Add(sum(optional_group_allocation[(t, g)] * seats[(t, g)] for g in all_individuals) <=
            #           group_defs['MaxNumNonComp'][t])

        # Add 'maximum pairwise relatedness' constraint for corrals whose MaxPR value != -1.
        if group_defs['MaxPR'][t] != -1:
            model.Add(
                sum(
                    connections[g1][g2] * same_corral[(g1, g2, t)]
                    for g1 in range(num_individuals - 1)
                    for g2 in range(g1 + 1, num_individuals)
                    if connections[g1][g2] > group_defs['MaxPR'][t]
                ) < 1
            )

    # Link colocated with seats
    for g1 in range(num_individuals - 1):
        for g2 in range(g1 + 1, num_individuals):
            for t in all_groups:
                # Link same_corral and seats.
                model.AddBoolOr(
                    [
                        seats[(t, g1)].Not(),
                        seats[(t, g2)].Not(),
                        same_corral[(g1, g2, t)],
                    ]
                )
                model.AddImplication(same_corral[(g1, g2, t)], seats[(t, g1)])
                model.AddImplication(same_corral[(g1, g2, t)], seats[(t, g2)])

            # Link colocated and same_table.
            model.Add(
                sum(same_corral[(g1, g2, t)] for t in all_groups) == colocated[(g1, g2)]
            )

    # Symmetry breaking. First individual is placed in the first group.
    print("Start of initial corral allocations.")
    for g1 in range(len(allocate_first_group)):
        if allocate_first_group[g1] != -1:
            print("\t\tIndividual %i is allocated to corral %i" % (g1, allocate_first_group[g1]))
            model.Add(seats[(allocate_first_group[g1], g1)] == 1)
    print("End of initial corral allocations.")

    paramstring = "%i,%i,%i" % (num_groups, objective_function, num_individuals)

    # Solve model.
    solver = cp_model.CpSolver()
    solution_printer = CobreederPrinter(seats, names, num_groups, num_individuals, paramstring, unique_id)

    # solver.parameters.max_time_in_seconds = 1.0
    # solver.parameters.log_search_progress = True
    # solver.parameters.num_workers = 1
    # solver.parameters.fix_variables_to_their_hinted_value = True

    # TODO check what solver.SearchForAllSolutions() is.
    status = solver.Solve(model, solution_printer)

    print("\nCOBREEDER-COMPLETION [%s]" % args.unique_run_id)

    # Print solution statistics
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        print("\tStatistics - Optimal") if status == cp_model.OPTIMAL else print("Statistics - Feasible")
        print("\t- conflicts    : %i" % solver.NumConflicts())
        print("\t- branches     : %i" % solver.NumBranches())
        print("\t- wall time    : %f s" % solver.WallTime())
        print("\t- num solutions: %i" % solution_printer.num_solutions())
        # Save best solution to CSV if desired
        save = input("\nSave final solution to CSV? (Y/N): ")
        if save.lower() == 'y':
            save_solution_csv(args)
    else:
        print("No solution found.")


def main(argv: Sequence[str]) -> None:
    parser = argparse.ArgumentParser(prog='CoBreeder_for_GhostWolves',
                                     description='Group-Living Captive Breeding Solver.', add_help=True)
    subparsers = parser.add_subparsers(help='sub-command help')
    run_parser = subparsers.add_parser('run')
    run_parser.add_argument('individuals_file', type=str,
                            help='List of individuals')
    run_parser.add_argument('pairwise_relatedness_file', type=str,
                            help='Scaled pairwise relatedness matrix')
    run_parser.add_argument('group_file', type=str,
                            help='Group specification')
    run_parser.add_argument("obj_function", type=str,
                            choices=[e.name for e in CobreederObjectiveFunction],
                            help='String for objective function')
    run_parser.add_argument("unique_run_id", type=str,
                            help='Unique identifier for each solver run.')
    run_parser.add_argument("weight_alleles", type=int,
                            help='Weight for alleles.')
    run_parser.add_argument("weight_pr", type=int,
                            help='Weight for PR.')
    run_parser.add_argument("total_individuals", type=int,
                            help='The minimum number of individuals allocated to a solution.')
    run_parser.add_argument("exclude_disallow", type=str, choices=["EX", "ALL"],
                            help='Exclude individuals or specify disallowed pairings with "EX".')
    run_parser.add_argument("prio_calc_threshold", type=int, choices=range(0, 101),
                            help='Threshold for priority calculation. 0 to disable priority calculation and use'
                                 'manual priority assignments only.')
    run_parser.add_argument('subst', nargs='?', default=0, type=int)

    args = parser.parse_args()
    solve_with_discrete_model(args)


if __name__ == "__main__":
    sys.setrecursionlimit(0x100000)
    app.run(main)
