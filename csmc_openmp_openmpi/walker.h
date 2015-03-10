/** \file walker.h
    This is the Walker class which handles the actual random evolution of the ensemble.
 
    The class is templated on the unsigned integer data type chosen to represent a state in the configuration space. This can be uint32_t, uint64_t, or uint128_t depending on the architecture
    and compiler support of the host system.
*/

#ifndef WALKER_H
#define WALKER_H

#include <algorithm>
#include <vector>
#include <utility>
#include "randomc.h"
#include "mersenne.cpp"
#include "pairinghamiltonian.h"
#include "wavefunction.h"

/**
    Walker class. walker type must be an unsigned type like uint32_t or uint64_t!
 
    @param levels integer value for the number of levels.
    @param NP integer value for the number of pairs.
    @param nbits integer value for the number of bits in the walker data type.
    @param wf const reference to a Wavefunction class.
    @param H const reference to a Hamiltonian class.
    @param rgenerators an std::vector of mersenne generator objects. One per a thread.
*/
template <typename walker>
class Walker {
    
public:
    Walker(unsigned, unsigned, unsigned, double, const Wavefunction<walker>&, const Hamiltonian<walker>&);
    ~Walker();
    void initialize_walkers(std::vector< std::pair< std::vector<walker>, double> >&);
    void random_walk(unsigned, std::vector< std::pair<std::vector<walker>, double> >&); //walk all the walkers
                                
private:
    
    void random_mv(unsigned&, unsigned&, std::vector<walker>&); //do a random move in a state and change the bag
    void single_step(std::pair< std::vector<walker>, double >&); //do a single step with a walker, maybe useful if we want to save the history of a walker

    unsigned levels;
    unsigned NP;
    size_t nbits;
    
    const Wavefunction<walker>& wf;
    const Hamiltonian<walker>& H;
    
    std::vector<CRandomMersenne> rgenerators;

};

/**
    Constructor for the walker class.
 
    @param nthreads integer value indicating the number of threads.
    @param nlevels integer value indicating the total number of levels.
    @param np integer value indicating the total number of pairs.
    @param seed floating point initial seed value for the mersenne generator.
    @param rwf const reference to Wavefunction class.
    @param rH const reference to the Hamiltonian class.
*/
template <typename walker>
Walker<walker>::Walker(unsigned nthreads, unsigned nlevels, unsigned np, double seed, const Wavefunction<walker>& rwf, const Hamiltonian<walker>& rH)
: levels(nlevels), NP(np), wf(rwf), H(rH), rgenerators(nthreads, CRandomMersenne(seed)) {
    
    walker tmp;
    nbits = 8 * sizeof(tmp); // set the number of bits
   
    for(size_t i = 0; i < rgenerators.size(); ++i) {
        CRandomMersenne Randnum(i*seed);
        rgenerators[i] = Randnum;
    }

}

template <typename walker>
Walker<walker>::~Walker() {}

/**
    This function loops through the entire ensemble of walkers, from 0 to numsteps, and performs the random walk.
 
    @param numsteps unsigned int containing the number of single steps that occur.
    @param walkers const reference to an std::vector of std::pair containing all of the current generation of walkers.
*/
template <typename walker>
void Walker<walker>::random_walk(unsigned numsteps, std::vector< std::pair< std::vector<walker>, double > >& walkers) {
    
    for(size_t ss = 0; ss < numsteps; ++ss) {
    
#pragma omp parallel for
        for(size_t w = 0; w < walkers.size(); ++w) {
            single_step(walkers[w]);
        }
        
    }
    
}

/**
    This function performs a single step for a walker. 
 
    @param rstate reference to an std::pair defining a random walker and it's amplitude.
*/
template <typename walker>
void Walker<walker>::single_step(std::pair< std::vector<walker>, double >& rstate){
    
    // Here are our two probabilities. Diagonal move and off diagonal moves.
    double Pdiagonal = 0;
    double Poffdiagonal = 0;
    
    if(NP == levels) {
        Pdiagonal = 1.0;
        Poffdiagonal = 0.0;
    }
    else {
        Pdiagonal = H.Probdiagonal(rstate.first);
        Poffdiagonal = H.Proboffdiagonal(rstate.first);
    }
   
    double randdoub = rgenerators[ omp_get_thread_num() ].Random();

    if(Pdiagonal > randdoub) { // diagonal transiton
        (rstate.second) *= (H.Diagmatrixelement(rstate.first) / Pdiagonal);
    }
    else{ // off diagonal
        unsigned i, j; // indices to track pair transfer
        random_mv(i, j, rstate.first); // random move in the state
        (rstate.second) *= (H.Offdiagmatrixelement(i, j) / Poffdiagonal);
    }
    
}

/**
    This function performs a random move within a single configuration. 
 
    @param ri reference to an unsigned integer defining the initial occupied level where a pair is removed.
    @param rj reference to an unsigned int defining to final level where a pair is moved to.
    @param rstate reference to an std::vector containing the configuration of the current walker.
*/
template <typename walker>
void Walker<walker>::random_mv(unsigned& ri, unsigned& rj, std::vector<walker>& rstate) { // inputs will be modified!
    
    walker bit = 0;
    unsigned index1, index2;
    
    do {
        ri = rgenerators[ omp_get_thread_num() ].IRandomX(0, levels - 1);
        index1 = unsigned(double(ri) / double(nbits));
        bit = rstate[index1] & (walker(1) << (ri - index1 * nbits));
    } while(bit == 0); // select initial occupied level
    
    do {
        rj = rgenerators[ omp_get_thread_num() ].IRandomX(0, levels - 1);
        index2 = unsigned(double(rj) / double(nbits));
        bit = rstate[index2] & (walker(1) << (rj - index2 * nbits));
    } while(bit != 0); // select final level that is unoccupied
    
    // This moves the pair from the init level to the fin level.
    rstate[index1] ^= walker(1) << (ri - index1 * nbits);
    rstate[index2] ^= walker(1) << (rj - index2 * nbits); 
    
    
}

/**
    This function is used to initialize the entire ensemble of walkers before the first step of the random walk.
    All walkers are initialized by placing all pairs into the lowest energy levels of the state.
 
    @param walkers reference to an std::vector of std::pairs containing all the walkers in the ensemble.
*/
template <typename walker>
void Walker<walker>::initialize_walkers(std::vector< std::pair< std::vector<walker>, double> >& walkers) { 
    
    unsigned int size = (walkers[0].first).size();
    
#pragma omp parallel for
    for(size_t added = 0; added < walkers.size(); ++added) {
        
        std::vector<walker> state(size, 0);
        unsigned index = 0;
        for(size_t placed = 0; placed < NP; ++placed) {
            state[index] += pow(2.0, placed);
            if(placed == nbits) { ++index; } //state[index] full so go to next one
        }
        std::pair<std::vector<walker>, double> tmp(state, 1.0);
        walkers[added] = tmp;
        
    }

    
}

#endif

