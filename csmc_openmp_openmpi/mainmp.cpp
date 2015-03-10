/** \file mainmp.cpp
 This is the main file for the parallel version of the PMC software. Parallelization is implemented using both OpenMPI and OpenMP.
 Also, c++11 support is required and the g++ -std=c++11 flag must be set during compilation.
 
 It is assumed that the software will be run on a 64 bit x86 system. Thus, the unsigned data type used to represent walkers
 is defined as: "#define walker_type uint64_t".
 If this is not the case then this must be changed to reflect the word length of the host system. (e.g. uint32_t)
 
 An input file is required containing the parameters for the Monte Carlo evolution and the physical parameters needed
 to define the number of particles, levels, and various aspects of the pairing hamiltonian.
 
 Format of input file assumes:
 1st line: number of walkers (numwalkers), number of steps (numsteps), length of each step (steplength), max iterations (max_iterations)
 2nd line: number of pair orbitals for |n> (levels), number of pairs (NP)
 3rd line: s.p. energy of each orbital (ee)
 all other: i, j matrix indices and Gij values (Gij)
 
 A sample input file is included labeled 2levels.inpt.
 
 For the first iteration of the algorithm (equilibration step) the evolution will occur as follows. Every independent walker is initialized
 and then randomly walked for numsteps of length steplength. The variable steplength indicates how many random moves will take place during
 each numstep before the walker is saved into the wave function. After this occurs for every walker the entire wave function is then normalized
 and a new generation of walkers is selected via importance sampling from the previous generation. After the last numstep estimates of the
 ground state energy and occupation numbers are calculated.
 
 For iterations 2 to max_iterations, numstep is set to 1. Every independent walker is initialized via importance sampling from the result of the
 previous iteration. The walkers are then walked for steplength and additional estimates of the energy and occupation numbers are calculated.
 Every estimate is saved so that the averages and statistical errors can be computed.
 
 The final output of the software will be the average ground state energy and it's statistical error and the occupation numbers for every level.
 */


#define walker_type uint64_t //Here is the unsigned int type used for the walkers

#include <fstream>
#include <iostream>
#include <cmath>
#include <ctime>
#include <utility>
#include <vector>
#include <string>
#include <omp.h>
#include <mpi.h>
#include "pairinghamiltonian.h"
#include "wavefunction.h"
#include "walker.h"
#include "operators.h"

//Global function prototypes
double average(const std::vector<double>& estimates);
double statistical_error(double avg, const std::vector<double>& estimates);
void averageoccs(int, std::vector<double>&);

int main(int argc, char **argv){
    
    //Check for interaction file
    if(argc != 2){
        std::cerr << "Usage: " << argv[0] << " <interaction file> " << std::endl;
        return 1;
    }
    
    int available_procs = 4; //omp_get_num_procs();
    omp_set_num_threads(available_procs);
    
    int numprocs, rank;
    MPI_Status Stat;
    MPI_Comm comm = MPI_COMM_WORLD;
    //int required = MPI_THREAD_MULTIPLE;
    //int provided;
    
    int error = MPI_Init(&argc, &argv);
    // int error = MPI_Init_thread(&argc, &argv, required, &provided);
    if(error != MPI_SUCCESS){
        std::cout << "Error starting MPI program. Terminating." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, error);
    }
    
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    std::ifstream infile;
    unsigned levels = 0, NP = 0; //number of levels and single particle energies
    uint64_t numwalkers = 0;
    unsigned numsteps = 0; //save wf every numsteps
    unsigned steplength = 0; //length of each iteration
    unsigned max_iterations = 0; //number of estimates we want
    
    if(rank == 0){
        
        infile.open(argv[1]);
        infile >> numwalkers;
        infile >> numsteps;
        infile >> steplength;
        infile >> max_iterations;
        infile >> levels;
        infile >> NP;
        
    }
    
    MPI_Bcast(&numwalkers, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&steplength, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&levels, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NP, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    std::vector<double> ee(levels, 0);  //single particle energies
    std::vector< std::vector<double> > Gij(levels, ee); //G matrix
    
    MPI_Barrier(comm);
    
    if(rank == 0) {
        
        double tmp2;
        for(size_t i = 0; i < levels; ++i){
            infile >> tmp2;
            ee[i] = tmp2;
        }
        
        for(size_t i = 0; i < levels; ++i){
            size_t x, y;
            for(size_t j = 0; j < levels; ++j){
                infile >> x;
                infile >> y;
                infile >> Gij[i][j];
            }
        }
        
        infile.close();
        
        for(size_t dest = 1; dest < unsigned(numprocs); ++dest){
            for(size_t i = 0; i < levels; ++i){
                for(size_t j = 0; j < levels; ++j){
                    MPI_Send(&Gij[i][j], 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
                }
            }
        }
        
    }
    else {
        
        for(size_t i = 0; i < levels; ++i){
            for(size_t j = 0; j < levels; ++j){
                double tmpv = 0;
                MPI_Recv(&tmpv, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &Stat);
                Gij[i][j] = tmpv;
            }
        }
        
    }
    
    MPI_Bcast(&ee[0], levels, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(comm);
    
    max_iterations /= numprocs + 1; // all processes produce at least one estimate of the ground state
    
    //At this point everyone should have all the imput data
    CRandomMersenne Randnum(time(NULL) * rank + time(NULL)); //this is going to be used to seed each thread
    double seed = Randnum.Random() * Randnum.IRandomX(0, 2147483647); //seed for the random number generator
    
    //Class constructors here
    Hamiltonian<walker_type> H(levels, NP, ee, Gij);
    Wavefunction<walker_type> wf(available_procs, NP, levels, seed);
    Walker<walker_type> ww(available_procs, levels, NP, seed, wf, H);
    Operators<walker_type> O(levels, NP, ee, Gij, wf);
    
    ee.clear();
    Gij.clear();
    
    //Now for the fun stuff!
    walker_type tmp;
    size_t nbits = 8 * sizeof(tmp); //set the number of bits
    size_t size_needed;
    if(levels == nbits) { size_needed = 1; }
    else if(levels != nbits && levels % nbits == 0) { size_needed = unsigned(double(levels) / double(nbits)); }
    else { size_needed = unsigned(double(levels) / double(nbits)) + 1; } //size of the state needed depends on the number of bits in our representation
    
    std::vector<walker_type> tmpstate(size_needed, 0);
    std::pair< std::vector<walker_type>, double> tmpwalker(tmpstate, 0);
    std::vector< std::pair<std::vector<walker_type>, double> > walkers(numwalkers, tmpwalker);
    
    ww.initialize_walkers(walkers);
    
    std::vector<double> Eestimates(max_iterations * numprocs, 0);
    std::vector<double> occupationn; //occupations
    std::vector< std::vector<double> > ioccs;
    size_t iteration = 0;
    
    do{

        for(size_t s = 0; s < numsteps; ++s){
    
            if(iteration == 0 && s != 0){
                wf.weighted_choice(walkers);
                wf.clear();
            }
            else if(iteration != 0){
                wf.weighted_choice(walkers);
                wf.clear();
            }
     
            ww.random_walk(steplength, walkers);
            wf.insert(walkers);
     
        }
   
        Eestimates[iteration] = O.Egs();
        
        occupationn.clear();
        occupationn.resize(levels, 0);
        O.occupations(occupationn);
        ioccs.push_back(occupationn);
        
        ++iteration;
        numsteps = 1; // After the first long walk, do short walks.
        
    } while(iteration < max_iterations);
    
    MPI_Barrier(comm);
    
    MPI_Gather(&Eestimates[0], max_iterations, MPI_DOUBLE, &Eestimates[rank * max_iterations], max_iterations, MPI_DOUBLE, 0, comm);
    
    std::vector<double> osums(levels, 0);
    for(size_t i = 0; i < ioccs.size(); ++i){
        for(size_t j = 0; j < levels; ++j){
            osums[j] += ioccs[i][j];
        }
    }
    
    ioccs.clear();
    occupationn.clear();
    occupationn.resize(levels, 0);
    
    for(size_t i = 0; i < levels; ++i){
        MPI_Reduce(&osums[i], &occupationn[i], 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    }
    
    if(rank == 0){
        
        averageoccs(max_iterations * numprocs, occupationn);
        
        double EE = average(Eestimates);
        
        std::cout << "E = " << EE << " +/- " << statistical_error(EE, Eestimates) << std::endl;
      
        std::cout << "occupations " << std::endl;
        for(size_t i = 0; i < occupationn.size(); ++i){
            std::cout << i << " " << occupationn[i] << std::endl;
        } std::cout << std::endl;
        
    }
    
    MPI_Finalize();
    
    return 0;
}

double average(const std::vector<double>& estimates){
    
    double avg = 0;
    for(size_t i = 0; i < estimates.size(); ++i) { avg += estimates[i]; }
    avg /= double(estimates.size());
    return avg;
    
}

double statistical_error(double avg, const std::vector<double>& estimates){
    
    double staterr = 0;
    for(size_t i = 0; i < estimates.size(); ++i) { staterr += pow(estimates[i] - avg, 2.0); }
    staterr = sqrt((1.0 / (double(estimates.size()) * (double(estimates.size()) - 1.0))) * staterr);
    return staterr;
    
}

void averageoccs(int total_estimates, std::vector<double>& roccs){
    
    for(size_t i = 0; i < roccs.size(); ++i){
        roccs[i] /= double(total_estimates);
    }
    
}

