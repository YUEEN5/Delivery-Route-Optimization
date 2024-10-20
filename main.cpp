#include <iostream>
#include <vector>
#include <set>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <algorithm>
using namespace std;

const int V = 31;
const int POPSIZE = 30;
const int GENE = 20;
const int WEIGHT[V] = { 0, 31, 29, 53, 38, 56, 45, 30, 15, 20, 27, 32, 25, 18, 9, 17, 30, 18, 25, 39, 22, 10, 35, 20, 11, 20, 17, 28, 40, 55, 31};//inside [] put value or don't pur also can
const int CAPACITY = 250;
const int NUMOFDRIVERS = 5;
const char DRIVERS[NUMOFDRIVERS] = { 'A', 'B', 'C', 'D', 'E'};
const float PENALTY = 0.05;
const double XOVER_PROB = 0.9; //probability of parents have children, coz not every parents will have children
const double MUTATE_PROB = 0.2; //use 0.9 to test, if the mutate works, then only change back to 0.1
const int MAX_GENERATION = 30;
//variables
//int chromosome1[POPSIZE][GENE];
//char chromosome2[POPSIZE][GENE];
//int chromosome[POPSIZE][GENE];
float fitness[POPSIZE]; //index for every chromosome
int indexParents[2];
//1. A temporary variable to store the new chromosome data structure
int newChromosome[POPSIZE][GENE];
//2. A counter for counting how many new chromosomes there are so far (array counter)
//int countNewChromo;
//Declare three output files for average fitness, best fitness and best chromosome
ofstream avgFitnessFile, bestFitnessFile, bestChromoFile;
//Declare avgFitness and bestFitness variables and give an initial value
float avgFitness, bestFitness = -1;
//int bestChromo1[GENE]; //declare here
//char bestChromo2[GENE];
//chromo1=location; chromo2=driver;
struct Gene {
    int chromosome1[GENE];
    char chromosome2[GENE];
};

Gene population[POPSIZE];
Gene parents[2], offsprings[2];
Gene bestChromo;

int map[][V] = {
    {0,	21,	35,	13,	16,	2,	7,	15,	30,	27,	31,	25,	30,	21,	16, 20,	14,	8,	25,	11,	21,	10,	19,	20,	7,	29,	18,	30,	22,	26,	18},
    {21, 0,	22,	9,	14,	27,	34, 18,	23,	40,	17,	27,	21,	19,	32, 38,	40,	12,	3,	39,	26,	32,	18,	9,	14,	7,	8,	29,	14,	41,	21},
    {35, 22, 0,	23,	36,	21,	25,	31,	30,	27,	6,	18,	18,	5,	29, 23,	19,	19,	35,	39,	13,	38,	20,	35,	38,	18,	33,	31,	12,	36,	30},
    {13, 9,	23,	0,	25,	18,	20,	21,	23,	24,	17,	11,	19,	15,	17, 16,	32,	27,	31,	38,	14,	41,	20,	30,	39,	37,	4,	40,	3,	32,	13},
    {16, 14, 36, 25,0,	15,	39,	41,	24,	15,	9,	14,	31,	32,	24, 20,	6,	10,	6,	25,	15,	20,	29,	2,	21,	6,	11,	14,	40,	40,	7},
    {2,	27,	21,	18,	15,	0,	14,	10,	11,	22,	15,	21,	22,	20,	17, 35,	24,	11,	1,  33,	38,	3,	21,	15,	36,	41,	9,	16,	12,	11,	26},
    {7,	34,	25,	20,	39,	14,	0,	32,	20,	33,	23,	32,	33,	10,	11, 10,	26,	14,	14,	15,	22,	27,	4,	32,	13,	2,	24,	11,	4,	35,	12},
    {15, 18, 31, 21,41,	10,	32,	0,  24,	28,	24,	8,	12,	15,	28, 2,	22,	34,	9,	11,	36,	11,	9,	22,	10,	18,	36,	10,	40,	33,	17},
    {30, 23,30,	23,	24,	11,	20,	24,	0,	19,	12,	40,	3,	35,	21, 12,	24,	10,	28,	3,	8,	3,	21,	28,	13,	6,	34,	5,	10,	32,	18},
    {27, 40,27,	24,	15,	22,	33,	28,	19,	0,	21,	11,	18,	26,	18, 5,	17,	12,	33,	9,	15,	26,	18,	33,	8,	40,	28,	14,	18,	16,	20},
    {31, 17, 6, 17,	9,	15,	23,	24,	12,	21,	0,	23,	24,	17,	29, 30,	12,	24,	11,	39,	11,	22,	32,	10,	36,	30,	24,	11,	21,	15,	5},
    {25, 27, 18, 11, 14, 21, 32, 8, 40,	11,	23,	0,	29,	12,	9, 37,  11,	17,	12,	1,	25,	12,	21,	7,	24,	26,	5,	2,	8,	32,	30},
    {30, 21, 18, 19, 31, 22, 33, 12, 3,	18,	24,	29,	0,	34,	23, 33,	20,	39,	18,	27,	32,	37,	11,	12,	20,	23,	32,	11,	17,	20,	21},
    {21, 19, 5,	15,	32,	20,	10,	15,	35,	26,	17,	12,	34,	0,	17, 9,	23,	20,	23,	27,	34,	22,	24,	23,	37,	4,	27,	36,	19,	4,	19},
    {16, 32, 29, 17, 24, 17, 11, 28, 21, 18, 29, 9, 23,	17,	0, 17,	13,	8,	9,	6,	41,	20,	8,	32,	27,	35,	13,	14,	12,	23,	10},
    {20, 38, 23, 16, 20, 35, 10, 2,	12,	5,	30,	37,	33,	9,	17,	0,	30,	21,	13,	10,	6,	18,	30,	25,	29,	8,	32,	26,	42,	37,	14},
    {14, 40, 19, 32, 6,	24,	26,	22,	24,	17,	12,	11,	20,	23,	13,	30,	0,  42,	7,	4,	13,	13,	3,	23,	29,	32,	17,	21,	39,	9,	20},
    {8,	12,	19,	27,	10,	11,	14,	34,	10,	12,	24,	17,	39,	20,	8,	21,	42,	0,	5,	18,	8,	18,	16,	31,	3,	14,	13,	12,	20,	6,	17},
    {25, 3,	35,	31,	6,	1,	14,	9,	28,	33,	11,	12,	18,	23,	9,	13,	7,	5,	0,	14,	35,	31,	4,	41,	32,	8,	5,	19,	9,	26,	19},
    {11, 39, 39, 38, 25, 33, 15, 11, 3,	9,	39, 1,	27,	27,	6,	10,	4,	18,	14,	0,	17,	32,	32,	37,	18,	13,	27,	33,	30,	22,	30},
    {21, 26, 13, 14, 15, 38, 22, 36, 8,	15,	11,	25,	32,	34,	41,	6,	13,	8,	35,	17,	0,	35,	15,	27,	41,	38,	9,	10,	25,	3,	26},
    {10, 32, 38, 41, 20, 3,	27,	11,	3,	26,	22,	12,	37,	22,	20,	18,	13,	18,	31,	32,	35,	0,	24,	16,	9,	10,	29,	23,	9,	16,	21},
    {19, 18, 20, 20, 29, 21, 4,	9,	21,	18,	32,	21,	11,	24,	8,	30,	3,	16,	4,	32,	15,	24,	0,	34,	10,	4,	17,	26,	40,	31,	11},
    {20, 9,	35,	30,	2,	15,	32,	22,	28,	33,	10,	7,	12,	23,	32,	25,	23,	31,	41,	37,	27,	16,	34,	0,	32,	20,	9,	13,	14,	11,	33},
    {7,	14,	38,	39,	21,	36,	13,	10,	13,	8,	36,	24,	20,	37,	27,	29,	29,	3,	32,	18,	41,	9,	10,	32,	0,	25,	28,	32,	16,	25,	36},
    {29, 7,	18,	37,	6,	41,	2,	18,	6,	40,	30,	26,	23,	4,	35,	8,	32,	14,	8,	13,	38,	10,	4,	20,	25,	0,	27,	13,	36,	21,	31},
    {18, 8,	33,	4,	11,	9,	24,	36,	34,	28,	24,	5,	32,	27,	13,	32,	17,	13,	5,	27,	9,	29,	17,	9,	28,	27,	0,	30,	22,	7,	40},
    {30, 29, 31, 40, 14, 16, 11, 10, 5,	14,	11,	2,	11,	36,	14,	26,	21,	12,	19,	33,	10,	23,	26,	13,	32,	13,	30,	0,	12,	10,	9},
    {22, 14, 12, 3,	40,	12,	4,	40,	10,	18,	21,	8,	17,	19,	12,	42,	39,	20,	9,	30,	25,	9,	40,	14,	16,	36,	22,	12,	0,	29,	18},
    {26, 41, 36, 32, 40, 11, 35, 33, 32, 16, 15, 32, 20, 4,	23,	37,	9,	6,	26,	22,	3,	16,	31,	11,	25,	21,	7,	10,	29,	0,	20},
    {18, 21, 30, 13, 7,	26,	12,	17,	18,	20,	5,	30,	21,	19,	10,	14,	20,	17,	19,	30,	26,	21,	11,	33,	36,	31,	40,	9,	18,	20,	0 }

};

void initializePopulation() {
    //Chromosome chromo; //cannot put here
    srand(time(NULL)); // Seed the random number generator
    for (int i = 0; i < POPSIZE; i++) {
        Gene& chromo = population[i]; // Access the global array
        //chromo.chromosome1[0] = 1; // Set the starting point to 1
        set<int> uniqueGenes; // Set to track unique genes in a chromosome
        //uniqueGenes.insert(1); // Insert the starting point into the set
        int geneCount = 0; // Start from 1 because the first gene is already set
        int gene;
        while (geneCount < GENE) {
            do {
                gene = rand() % (V + 1); // Generate a gene
            } while (gene == 1);
            if (gene == 0 || uniqueGenes.find(gene) == uniqueGenes.end()) { // Ensure the gene is unique if it is not 0
                uniqueGenes.insert(gene);
                chromo.chromosome1[geneCount] = gene;
                geneCount++;
            }
        }
        for (int j = 0; j < GENE; j++)
        {
            chromo.chromosome2[j] = DRIVERS[rand() % NUMOFDRIVERS]; // Randomly assign a driver
        }
    }
}



void printChromosome() {
    //Gene chromo;
    for (int i = 0; i < POPSIZE; i++)
    {
        Gene& chromo = population[i]; // Access the global array
        cout << "\tChromosome " << i << ": START -> ";
        for (int j = 0; j < GENE; j++)
        {
            cout << chromo.chromosome1[j] << " -> ";
        }
        cout << "END";
        cout << endl;
        cout << "\t\t      START -> ";
        for (int j = 0; j < GENE; j++)
        {
            cout << chromo.chromosome2[j] << " -> ";
        }
        cout << "END";
        cout << endl;
    }
}

void evaluateChromosome() {
    int totalW, totalD, driverUsed, numOfLocation;
    float totalPenalty;
    for (int i = 0; i < POPSIZE; i++) {

        Gene& chromo = population[i];
        vector<int> paths[NUMOFDRIVERS];
        int weights[NUMOFDRIVERS] = { 0, 0, 0, 0, 0};
        int distances[NUMOFDRIVERS] = { 0, 0, 0, 0, 0};
        driverUsed = 0;
        totalD = 0;
        totalW = 0;
        numOfLocation = 0;
        totalPenalty = 0;

        for (int j = 0; j < GENE; j++) {
            if (chromo.chromosome1[j] != 0) {
                int driverIndex = chromo.chromosome2[j] - 'A'; //to convert A, B, C to 0, 1, 2 using ASCII
                paths[driverIndex].push_back(chromo.chromosome1[j]);
                numOfLocation += 1;
            }
        }

        cout << "\n\tChromosome " << i << ": " << endl;

        //const char* pathNames[3] = { "A", "B", "C" };
        for (int d = 0; d < NUMOFDRIVERS; d++) {
            if (!paths[d].empty()) {
                cout << "\tPath " << DRIVERS[d] << ": ";
                for (size_t j = 0; j < paths[d].size(); j++) {
                    weights[d] += WEIGHT[paths[d][j] - 1];
                    cout << paths[d][j] << " -> ";
                }
                cout << "END";
                cout << "\tParcel Weights: " << weights[d];
                if (weights[d] > CAPACITY)
                {
                    totalPenalty += PENALTY;
                }
                // Calculate distance starting from the starting point (index 0)
                distances[d] = map[0][paths[d][0] - 1];

                // Calculate distances between subsequent locations in the path
                for (size_t j = 0; j < paths[d].size() - 1; j++) {
                    int loc1 = paths[d][j];
                    int loc2 = paths[d][j + 1];
                    distances[d] += map[loc1 - 1][loc2 - 1];
                }
                // Calculate distance from the last point go back to the starting point (index 0)
                distances[d] += map[paths[d].back() - 1][0];

                cout << "\tDistance Path " << DRIVERS[d] << ": " << distances[d] << endl;
                driverUsed += 1;
                totalD += distances[d];
            }
        }
        cout << "\n\tNumber of drivers: " << driverUsed << "\tTotal Distance: " << totalD << "\tTotal location visited: " << numOfLocation << "\tPenalty applied: " << totalPenalty;
        fitness[i] = ((numOfLocation * 1.0 / (V-1)) * 0.4) + ((1.0 / driverUsed) * 0.3) + ((1.0 / totalD) * 0.3) - totalPenalty;
        cout << "\tFitness: " << fitness[i] << endl << endl;
    }
}

void tournamentSelection() {
    int players[3];
    float max;
    //avoid get the same parents
    do {
        for (int i = 0; i < 2; i++)
        {
            //avoid same players being chosen
            do {
                max = -1;
                for (int j = 0; j < 3; j++) {
                    players[j] = rand() % POPSIZE;
                    if (fitness[players[j]] > max) {
                        max = fitness[players[j]];
                        indexParents[i] = players[j];
                    }
                }
                if (players[2] != players[1] && players[2] != players[0] && players[0] != players[1])
                {
                    for (int j = 0; j < 3; j++) {
                        cout << "\tPlayer " << j << ": chromosome " << players[j] << "\t\tFitness: " << fitness[players[j]] << endl;
                    }
                }
            } while (players[2] == players[1] || players[2] == players[0] || players[0] == players[1]);

            cout << "\n\tWinner: chromosome " << indexParents[i] << " with fitness value: " << max << endl << endl;

            //return best; // Return index of the selected parent
            if (indexParents[0] == indexParents[1])
            {
                cout << "\n\tSame parents are selected. The above set is no longer be used.\n";
            }
        }
    } while (indexParents[0] == indexParents[1]);

    cout << "\n\tSelected Parents:";
    for (int i = 0; i < 2; i++)
    {
        cout << "\n\tParent " << i + 1 << ": Chromosome " << indexParents[i] << endl;
        parents[i] = population[indexParents[i]];
        //Gene& parentsChromo[i] = parents[i];
        cout << "\t";
        for (int j = 0; j < GENE; j++)
        {
            cout << parents[i].chromosome1[j] << "  ";
        }
        cout << endl << "\t";
        for (int j = 0; j < GENE; j++)
        {
            cout << parents[i].chromosome2[j] << "  ";
        }
    }
}

void fixDuplicates(int chromosome[], int length) {
    set<int> seen;
    vector<int> duplicates;

    // Identify duplicates, but allow multiple 0s
    for (int i = 0; i < length; i++) {
        if (chromosome[i] != 0 && seen.find(chromosome[i]) != seen.end()) {
            duplicates.push_back(i);
        }
        else if (chromosome[i] != 0) {
            seen.insert(chromosome[i]);
        }
    }

    // Replace duplicates with unique random values
    for (int index : duplicates) {
        int newGene;
        do {
            newGene = rand() % (V + 1); // Generate a gene from 0 to V
        } while (newGene == 1 || (newGene != 0 && seen.find(newGene) != seen.end()));
        chromosome[index] = newGene;
        if (newGene != 0) {
            seen.insert(newGene);
        }
    }
}

void nPointCrossover() {
    int crossoverPoints[3];
    int n = 3;
    // Generate n random crossover points
    do {
        for (int i = 0; i < n; i++) {
            crossoverPoints[i] = rand() % GENE;
        }
    } while (crossoverPoints[0] == crossoverPoints[1] || crossoverPoints[0] == crossoverPoints[2] || crossoverPoints[1] == crossoverPoints[2]);

    // Sort the crossover points in ascending order
    sort(crossoverPoints, crossoverPoints + n);
    cout << "\n\tThe crossover points are: ";
    for (int i = 0; i < n; i++)
    {
        cout << crossoverPoints[i] << ", ";
    }
    cout << endl;
    // Initialize offspring with parents
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < GENE; j++)
        {
            offsprings[i].chromosome1[j] = parents[i].chromosome1[j];
            offsprings[i].chromosome2[j] = parents[i].chromosome2[j];
        }
    }
    cout << "\n\tBefore Crossover:" << endl;
    for (int i = 0; i < 2; i++)
    {
        cout << "\tChildren " << i + 1 << ": " << endl;
        cout << "\t";
        for (int j = 0; j < GENE; j++)
        {
            cout << offsprings[i].chromosome1[j] << offsprings[i].chromosome2[j] << " ";
        }
        cout << endl;
    }
    // Apply N-Point Crossover
    bool swap = false;
    int previousPoint = 0;
    for (int i = 0; i < n; i++) {
        int currentPoint = crossoverPoints[i];
        if (swap) {
            for (int j = previousPoint; j < currentPoint; j++) {
                std::swap(offsprings[0].chromosome1[j], offsprings[1].chromosome1[j]);
                std::swap(offsprings[0].chromosome2[j], offsprings[1].chromosome2[j]);
            }
        }
        swap = !swap;
        previousPoint = currentPoint;
    }

    // Handle the section after the last crossover point
    if (swap) {
        for (int j = previousPoint; j < GENE; j++) {
            std::swap(offsprings[0].chromosome1[j], offsprings[1].chromosome1[j]);
            std::swap(offsprings[0].chromosome2[j], offsprings[1].chromosome2[j]);
        }
    }
    cout << "\n\tAfter Crossover: " << endl;
    for (int i = 0; i < 2; i++)
    {
        cout << "\tChildren " << i + 1 << ": " << endl << "\t";
        for (int j = 0; j < GENE; j++)
        {
            cout << offsprings[i].chromosome1[j] << offsprings[i].chromosome2[j] << " ";
        }
        cout << endl;
    }

    // Fix duplicates in the offspring
    fixDuplicates(offsprings[0].chromosome1, GENE);
    fixDuplicates(offsprings[1].chromosome1, GENE);
    cout << "\n\tAfter Fixing repeated location: " << endl;
    for (int i = 0; i < 2; i++)
    {
        cout << "\tChildren " << i + 1 << ": " << endl << "\t";
        for (int j = 0; j < GENE; j++)
        {
            cout << offsprings[i].chromosome1[j] << offsprings[i].chromosome2[j] << " ";
        }
        cout << endl;
    }
}

void swapMutation() {

    // Perform mutation by swapping two genes if mutation condition is met
    vector<int> mutationPoints;
    for (int i = 0; i < 2; i++)
    {
        cout << "\n\tChildren " << i + 1 << ": ";
        // Display children chromosome before mutation

        int swapIndex1 = rand() % GENE;
        int swapIndex2 = rand() % GENE;
        while (swapIndex1 == swapIndex2) {
            swapIndex2 = rand() % GENE;
        }
        cout << "\n\tThe swap index are " << swapIndex1 << " and " << swapIndex2 << endl;

        cout << "\tBefore Mutation:\n\t";

        for (int j = 0; j < GENE; j++) {
            cout << offsprings[i].chromosome1[j] << offsprings[i].chromosome2[j] << " ";
        }
        cout << endl;

        // Swap the genes
        int temp = offsprings[i].chromosome1[swapIndex1];
        offsprings[i].chromosome1[swapIndex1] = offsprings[i].chromosome1[swapIndex2];
        offsprings[i].chromosome1[swapIndex2] = temp;

        // Swap the corresponding drivers
        char tempDriver = offsprings[i].chromosome2[swapIndex1];
        offsprings[i].chromosome2[swapIndex1] = offsprings[i].chromosome2[swapIndex2];
        offsprings[i].chromosome2[swapIndex2] = tempDriver;

        // Store mutation points for display
        mutationPoints.push_back(swapIndex1);
        mutationPoints.push_back(swapIndex2);

        // Display children chromosome after mutation
        cout << "\tAfter Mutation:\n\t";
        for (int j = 0; j < GENE; j++) {
            cout << offsprings[i].chromosome1[j] << offsprings[i].chromosome2[j] << " ";
        }
        cout << endl;

    }
}

void evaluateMutatedChromosome() {

    cout << "\nChildren Evaluation: " << endl;
    for (int i = 0; i < 2; i++)
    {
        vector<int> paths[NUMOFDRIVERS];
        int weights[NUMOFDRIVERS] = {0};
        int distances[NUMOFDRIVERS] = {0};
        int driverUsed = 0;
        int numOfLocation = 0;
        float totalPenalty = 0.0;
        for (int j = 0; j < GENE; j++) {
            if (offsprings[i].chromosome1[j] != 0) {
                int driverIndex = offsprings[i].chromosome2[j] - 'A'; // Convert A, B, C to 0, 1, 2 using ASCII
                paths[driverIndex].push_back(offsprings[i].chromosome1[j]);
                numOfLocation++;
            }
        }


        cout << "\tChildren " << i+1 << " chromosome " << indexParents[i] << " : " << endl;

        for (int d = 0; d < NUMOFDRIVERS; d++) {
            if (!paths[d].empty()) {
                cout << "\tPath " << DRIVERS[d] << ": ";
                for (size_t j = 0; j < paths[d].size(); j++) {
                    weights[d] += WEIGHT[paths[d][j] - 1];
                    cout << paths[d][j] << " -> ";
                }
                cout << "END";
                cout << "\tParcel Weights: " << weights[d];

                if (weights[d] > CAPACITY) {
                    totalPenalty += PENALTY;
                }

                // Calculate distance starting from the starting point (index 0)
                distances[d] = map[0][paths[d][0] - 1];

                // Calculate distances between subsequent locations in the path
                for (size_t j = 0; j < paths[d].size() - 1; j++) {
                    int loc1 = paths[d][j];
                    int loc2 = paths[d][j + 1];
                    distances[d] += map[loc1 - 1][loc2 - 1];
                }
                // Calculate distance from the last point back to the starting point (index 0)
                distances[d] += map[paths[d].back() - 1][0];

                cout << "\tDistance Path " << DRIVERS[d] << ": " << distances[d] << endl;
                driverUsed++;
            }
        }

        cout << "\n\tNumber of drivers: " << driverUsed << "\tTotal location visited: " << numOfLocation << "\tPenalty applied: " << totalPenalty;

        // Calculate fitness value
        float totalDistance = 0.0;
        for (int d = 0; d < NUMOFDRIVERS; d++) {
            totalDistance += distances[d];
        }
        float fitnessValue = ((numOfLocation * 1.0 / (V-1)) * 0.4) + ((1.0 / driverUsed) * 0.3) + ((1.0 / totalDistance) * 0.3) - totalPenalty;
        fitness[indexParents[i]] = fitnessValue;

        cout << "\tFitness: " << fitnessValue << endl << endl;
    }
}

void createNewPopulation() {
    // Identify the two chromosomes with the lowest fitness values
    vector<int> fitnessIndices(POPSIZE);
    iota(fitnessIndices.begin(), fitnessIndices.end(), 0);
    sort(fitnessIndices.begin(), fitnessIndices.end(), [](int a, int b) {
        return fitness[a] < fitness[b];
        });

    int lowestFitnessIndex1 = fitnessIndices[0];
    int lowestFitnessIndex2 = fitnessIndices[1];

    // Replace the two chromosomes with the lowest fitness values with the new offspring

    population[lowestFitnessIndex1] = offsprings[0];
    population[lowestFitnessIndex2] = offsprings[1];

    // Display the replacement
    cout << "\n\tReplacing chromosome " << lowestFitnessIndex1 << " with Child 1" << endl;
    cout << "\tReplacing chromosome " << lowestFitnessIndex2 << " with Child 2" << endl;

    // Display the new population
    cout << "\nNew Population Generated:\n";
    for (int i = 0; i < POPSIZE; i++) {
        cout << "\tChromosome " << i << ": ";
        for (int j = 0; j < GENE; j++) {
            cout << population[i].chromosome1[j] << population[i].chromosome2[j] << " ";
        }
        cout << endl;
    }
}

void calculateAverageFitness() {
    //1. Declare a variable for totalFitness, initialize to 0
    float totalFitness = 0;
    //2. For every chromosome
    for (int i = 0; i < POPSIZE; i++)
    {
        //2.1 Accumulate the fitness into totalFitness
        totalFitness += fitness[i];
    }
    //3. Divide the totalFitness with population size
    avgFitness = totalFitness / POPSIZE;
    //4. Print out the average to the screen
    cout << "\n\tAverage Fitness: " << avgFitness; //to the screen
    //5. Print out the average to an output file that keep average fitness
    avgFitnessFile << avgFitness << endl; //to the file
}

void recordBestFitness() {
    //1. Declare the bestChromosome data structure //declare above
    Gene& bestChromosome = bestChromo;
    //2. For each chromosome
    for (int i = 0; i < POPSIZE; i++)
    {
        Gene& chromo = population[i];

        //2.1. if (fitness current chromosome better than bestFitness){
        if (fitness[i] > bestFitness)
        {
            //2.1.1. bestFitness = fitness for the current chromosome
            bestFitness = fitness[i];
            //2.1.2. copy the chromosome to bestChromosome
            for (int j = 0; j < GENE; j++)
            {
                bestChromosome.chromosome1[j] = chromo.chromosome1[j];
            }
            cout << endl;
            for (int j = 0; j < GENE; j++)
            {
                bestChromosome.chromosome2[j] = chromo.chromosome2[j];
            }
        }
    }
    //3. Print the bestFitness and bestChromosome to the screen
    cout << "\tBSF: " << bestFitness;
    cout << "\n\tBest Chromosome: ";
    for (int j = 0; j < GENE; j++)
    {
        cout << bestChromosome.chromosome1[j] << " ";
        bestChromoFile << bestChromosome.chromosome1[j];
    }
    cout << endl;
    cout << "\t\t\t  ";
    for (int j = 0; j < GENE; j++)
    {
        cout << bestChromosome.chromosome2[j] << " ";
        bestChromoFile << bestChromosome.chromosome2[j];
    }
    cout << endl;
    bestChromoFile << endl;
    //4. Print the bestFitness and bestChromosome to two separate files
    bestFitnessFile << bestFitness << endl;

}

int main() {
    avgFitnessFile.open("avgFitness.txt"); //open in notepad, use .txt
    bestFitnessFile.open("bestFitness.txt");
    bestChromoFile.open("bestChromo.txt");
    cout << "GA starts..." << endl;
    cout << "\tInitialization" << endl;
    initializePopulation();
    for (int gen = 0; gen < MAX_GENERATION; gen++) {
        cout << "\n\tGeneration " << gen + 1 << endl;
        printChromosome();
        evaluateChromosome(); // Initial evaluation of population
        calculateAverageFitness();
        recordBestFitness();

        tournamentSelection();
        cout << "\n\nCrossover : N-Point\n";
        if ((double)rand() / RAND_MAX < XOVER_PROB) {
            cout << "\n\tCrossover happened!" << endl;
            nPointCrossover(); // 3-point crossover

        }
        else {
            cout << "\n\tCrossover does not happen!" << endl;
        }

        // Apply mutation to offspring1 and offspring2
        cout << "\nMutation: Swap Mutation" << endl;
        if ((double)rand() / RAND_MAX < MUTATE_PROB) {
            cout << "\n\tMutation happened!" << endl;
            swapMutation(); // Mutate Children
        }
        else
        {
            cout << "\n\tMutation does not happen!" << endl;
        }


        // Evaluate mutated offspring
        evaluateMutatedChromosome();

        // Create new population with the mutated offspring
        createNewPopulation();
        // Evaluate new population
        evaluateChromosome();

    }
    avgFitnessFile.close();
    bestFitnessFile.close();
    bestChromoFile.close();
    return 0;
}
