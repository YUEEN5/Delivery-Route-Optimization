# Parcel Delivery Route Optimization using Genetic Algorithm

## Overview
This project focuses on optimizing parcel delivery routes to reduce transportation cost and improve logistics efficiency. The problem is modeled as a constraint optimization problem and solved using a Genetic Algorithm (GA).

The system searches for efficient delivery routes while considering real-world operational constraints.

---

## Problem Statement
Parcel delivery companies must determine efficient delivery routes while managing constraints such as truck capacity and delivery coverage.

Traditional route planning methods may not scale well with complex delivery networks.

This project aims to develop an optimization algorithm capable of generating cost-efficient delivery routes.

---

## Solution
A **Genetic Algorithm (GA)** is used to search for optimal delivery routes by iteratively evolving candidate solutions. The algorithm evaluates route efficiency using a fitness function that considers travel distance, number of trucks used, and delivery coverage.

Constraint handling ensures that truck capacity limits are respected.

---

## Key Features
- Delivery route optimization using evolutionary algorithms
- Multi-objective fitness function
- Constraint handling for truck capacity
- Efficient exploration of large search spaces

---

## Methodology

### 1. Problem Representation
Delivery routes are represented as chromosomes in the genetic algorithm.

### 2. Population Initialization
An initial population of possible route configurations is generated.

### 3. Selection
**Tournament selection** is used to select high-quality solutions.

### 4. Crossover
**N-point crossover** is applied to combine parent solutions.

### 5. Mutation
**Swap mutation** introduces variation in routes to explore new solutions.

### 6. Fitness Evaluation
Solutions are evaluated based on:

- Total travel distance
- Number of trucks used
- Number of delivery locations served

Penalty functions ensure that truck weight limits (250kg) are not exceeded.

---

## Technologies Used

- Python
- Genetic Algorithm
- Optimization Algorithms
- NumPy
- Data Structures

---

## Results
Experimental testing demonstrated that the genetic algorithm consistently improves delivery route efficiency across multiple generations, producing stable and optimized solutions.

---

## Project Structure
