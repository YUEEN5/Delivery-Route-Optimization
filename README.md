Project Overview
This project addresses the logistical challenges of modern online shopping by optimizing delivery routes
. Modelled as a Constraint Optimization Problem (COP), it aims to minimize the number of trucks and total distance travelled while maximizing daily delivery locations
.
Evolutionary Algorithm Design
Genetic Algorithm (GA): Employs a nature-inspired GA to evolve efficient routing solutions
.
Core Mechanisms:
Chromosome Representation: Uses integer (locations) and alphabet (truck IDs) encodings
.
Operations: Implements Tournament selection, N-Point crossover, and Swap mutation to maintain genetic diversity
.
Fitness Function: A multi-objective function that weights distance (30%), location coverage (40%), and resource usage (30%)
.
Constraints: Includes a strict 250kg weight limit per truck, with penalties applied to infeasible solutions
.
Results
Testing over 50 generations demonstrates consistent convergence toward higher fitness values, providing stable and cost-effective delivery paths
.
