/** \file wavefunction.h
    This class represents the wave function for a paired, many body wave function in configuration space
    \f[ |\Psi \rangle = \sum _{i}^{\omega} \alpha_{i} |\textbf{n}_{i} \rangle \f]
    where \f$\omega\f$ is the total number of different configurations, \f$|\textbf{n}_{i} \rangle\f$ is
    the i-th configuration, and \f$\omega\f$ is the amplitude of that configuration.
 
    The underlying data structure of this class is an std::unordered_map that uses contains std::vector<walker>/double key/value pairs.
 
    The class is templated on the unsigned integer data type chosen to represent a state in the configuration space. This can be uint32_t, uint64_t, or uint128_t depending on the architecture
    and compiler support of the host system.
 
    The constructor requires as input the number of threads, the number of pairs (NP), the number of levels (numlevels), and the initial seed for the mersenne generator.
*/


#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <unordered_map>
#include <functional>
#include "randomc.h"
#include "am.h"

/**
    This class specializes the hashing function for the vector<walker> key and then hoists it into the std namespace.
    @param rstate const reference to the std::vector being hashed into std::unordered_map.
*/
namespace std {
    template <typename walker>
    class hash< vector<walker> > {
        public:
        
        size_t operator()(const std::vector<walker> &rstate) const {
            
            walker hash = 0;
            
            for(size_t i = 0; i < rstate.size(); i++) {
                hash ^= std::hash<walker>()(rstate[i]);
            }
            
            return hash;
            
        }

    };
}

/**
    Wavefunction class. walker type must be an unsigned type like uint32_t or uint64_t!
 
    @param NP integer value for the number of pairs.
    @param numlevels integer value for the number of levels.
    @param nbits nbits integer value for the number of bits in the walker data type.
    @param rgenerators an std::vector of mersenne generator objects. One per a thread.
    @param wave_function std::unordered_map container for the wave function.
*/
template <typename walker>
class Wavefunction {
    
public:
    Wavefunction(unsigned, unsigned);
    Wavefunction(unsigned, unsigned, unsigned, unsigned);
    ~Wavefunction();
    void random_choice(std::vector< std::pair< std::vector<walker>, double> >&); //select a (uniformly) random state from the wave function and its weight
    void weighted_choice(std::vector< std::pair< std::vector<walker>, double> >&);
    void weighted_choice_am(std::vector< std::pair< std::vector<walker>, double> >&);
    void insert(const std::vector< std::pair<std::vector<walker>, double> >&);
    double configuration_weight(const std::vector<walker>&) const; //return the weight from wave_function
    void clear();
    size_t size() const;
    void normalise_ln(); //normalise the entire wave function to the linear norm
    void normalise_sq(); //normalise the entire wave function to the square norm
    
    typedef typename std::unordered_map< std::vector<walker>, double >::const_iterator const_iterator;
    const_iterator begin() const { return wave_function.begin(); }
    const_iterator end() const { return wave_function.end(); }
    
private:
    
    unsigned NP;
    unsigned numlevels;
    std::size_t nbits;
    
    double factorial(int);
    
    std::vector<CRandomMersenne> rgenerators;
    std::unordered_map< std::vector<walker>, double > wave_function;
    
};

template <typename walker>
Wavefunction<walker>::Wavefunction(unsigned nthreads, unsigned seed) : NP(0), numlevels(0), rgenerators(nthreads, CRandomMersenne(seed)) {
  
    for(size_t i = 0; i < rgenerators.size(); ++i){
        CRandomMersenne Randnum(i*seed);
        rgenerators[i] = Randnum;
    }
						 
}

/**
    Constructor the Wavefunction class. 
 
    @param nthreads integer value representing the number of threads being used
    @param np integer value representing the number of pairs
    @param nlevels integer value representing the number of levels
    @param seed value representing the initial seed to be used to create the mersenne generator for each thread
*/
template <typename walker>
Wavefunction<walker>::Wavefunction(unsigned nthreads, unsigned np, unsigned nlevels, unsigned seed) : NP(np), numlevels(nlevels), rgenerators(nthreads, CRandomMersenne(seed)) {
    
    walker tmp;
    nbits = 8 * sizeof(tmp); //set the number of bits

    for(size_t i = 0; i < rgenerators.size(); ++i) {
         CRandomMersenne Randnum(i*seed);
         rgenerators[i] = Randnum;
    }
    
}

/**
    Destructor the Wavefunction class.
*/
template <typename walker>
Wavefunction<walker>::~Wavefunction() {}

/**
    This function normalizes the wave function following a linear normalization condition.
    \f[Norm = \sum_{i} A_{i}\f]
    where \f$A_{i}\f$ is the un-normalized amplitude of a state i in the wave function.
 
    The result is that each amplitude \f$\alpha_{i}\f$ in the normalized wave function represents 
    the classical probability to select the state via a random sampling.
*/
template <typename walker>
void Wavefunction<walker>::normalise_ln() {
    
    double norm = 0;
    
    for(auto it = wave_function.begin(); it != wave_function.end(); ++it){
        norm += (it -> second);
    }
    
    // use the value to normalise
    for(auto it = wave_function.begin(); it != wave_function.end(); ++it){
        (it -> second) /= norm;
    }
    
}

/**
    This function normalizes the wave function following a quadratic normalization condition.
    \f[Norm = \sqrt{\sum_{i} A_{i}^{2}}\f]
    where \f$A_{i}\f$ is the un-normalized amplitude of a state i in the wave function.
 
    The result is that each amplitude \f$\alpha_{i}\f$ in the normalized wave function represents
    the quantum probability to select the state via a random sampling.
*/
template <typename walker>
void Wavefunction<walker>::normalise_sq(){
    
    double norm = 0;
    
    for(auto it = wave_function.begin(); it != wave_function.end(); ++it) {
        norm += (it -> second) * (it -> second);
    }
    
    // use the value to normalise
    for(auto it = wave_function.begin(); it != wave_function.end(); ++it) {
        (it -> second) /= sqrt(norm);
    }
    
}

/**
    This function randomly selects, in parallel, a new ensemble of walkers from the previous wave_function.
    Each walker is represented by an std::pair containing an std::vector representing the states configuration
    and a double representing the amplitude. The random selection initializes the vector to a configuration 
    existing in the wave function and sets the walkers initial amplitude to the normalized amplitude from
    the wave function.
 
    @param walkers const reference to a two dimensional std::vector of size equal to the number of walkers.
*/
template <typename walker>
void Wavefunction<walker>::random_choice(std::vector< std::pair< std::vector<walker>, double> >& walkers){
    
    normalise_ln();
    
#pragma omp parallel for
    for(size_t added = 0; added < walkers.size(); ++added) {
        
        auto index = wave_function.begin();
        size_t randstate = rgenerators[ omp_get_thread_num() ].IRandomX(0, wave_function.size() - 1);
        advance(index, randstate); // advance it to a random state
        
        walkers[added].first = index -> first;
        walkers[added].second = index -> second;
        
    }
    
}

/**
    This function selects via the normalized probabilities, in parallel, a new ensemble of walkers from the previous wave_function.
    Each walker is represented by an std::pair containing an std::vector representing the states configuration
    and a double representing the amplitude. The weighted selection initializes the vector to a configuration
    existing in the wave function and sets the walkers initial amplitude to unity. Thus, most probable configurations are
    selected from the wave_function.
 
    @param walkers const reference to a two dimensional std::vector of size equal to the number of walkers.
*/
template <typename walker>
void Wavefunction<walker>::weighted_choice(std::vector< std::pair< std::vector<walker>, double> >& walkers) {
    
    normalise_ln();
    
#pragma omp parallel for
    for(size_t added = 0; added < walkers.size(); ++added) {
        
        auto index = wave_function.begin(); //define an index used to choose a state
        double num = rgenerators[ omp_get_thread_num() ].Random();
        while(num >= index -> second){
            num -= index -> second;
            ++index;
        }
        
        walkers[added].first = index -> first;
        walkers[added].second = 1.0;
        
    }
    
}

/**
 This function uses the alias method to select a new ensemble of walkers from the previous wave_function.
 Each walker is represented by an std::pair containing an std::vector representing the states configuration
 and a double representing the amplitude. The weighted selection initializes the vector to a configuration
 existing in the wave function and sets the walkers initial amplitude to unity. Thus, most probable configurations are
 selected from the wave_function.
 
 @param walkers const reference to a two dimensional std::vector of size equal to the number of walkers.
 */
template <typename walker>
void Wavefunction<walker>::weighted_choice_am(std::vector< std::pair< std::vector<walker>, double> >& walkers) {
    
    normalise_ln();
    
    std::vector<double> probabilities(wave_function.size(), 0);
    auto it2 = probabilities.begin();
    for(auto it = wave_function.begin(); it != wave_function.end(); ++it, ++it2) {
        *it2 = it -> second;
    }
    
    Alias_Method am(probabilities);
    
#pragma omp parallel for
    for(size_t added = 0; added < walkers.size(); ++added){
        
        auto index = wave_function.begin(); //define an index used to choose a state
        std::advance(index, am.next(rgenerators));
        
        walkers[added].first = index -> first;
        walkers[added].second = 1.0;
        
    }
    
}

/**
    This function inserts an ensemble of walkers in the wave function.
 
    @param walkers const reference to a two dimensional std::vector of size equal to the number of walkers.
*/
template <typename walker>
void Wavefunction<walker>::insert(const std::vector< std::pair< std::vector<walker>, double> >& walkers) {
    
    for(auto it = walkers.begin(); it != walkers.end(); ++it) {
        wave_function[ (*it).first ] += (*it).second;
    }
    
    //normalise_ln();
    
}

/**
    This function returns the current weight of a configuration in the wave function. If the configuration 
    exists then it's weight is returned otherwise 0 is returned. 
 
    @param crwalker a const reference to an std::vector representing a configuration \f$| \textbf{n}_{i} \rangle\f$
*/
template <typename walker>
double Wavefunction<walker>::configuration_weight(const std::vector<walker>& crwalker) const {
    
    auto it = wave_function.find(crwalker);
    
    if( it == wave_function.end() ) { return 0.0; }
    return it -> second;
    
}

/**
    The function returns the value of n!.
 
    @param n integer value to compute the factorial of
*/
template <typename walker>
double Wavefunction<walker>::factorial(int n){
    
    static std::vector<double> a(171, 0);
    static bool init = true;
    if(init) {
        init = false;
        a[0] = 1.0;
        for(int i = 1; i < 171; ++i) { a[i] = i * a[i - 1]; }
    }
    if (n < 0 || n > 170) { throw("factorial out of range"); }
    return a[n];
    
}

/**
    This function will clear the std::unordered_map containing the wave function.
*/
template <typename walker>
void Wavefunction<walker>::clear() {
    wave_function.clear();
}

/**
 This function will return the size of the std::unordered_map containing the wave function.
 */
template <typename walker>
size_t Wavefunction<walker>::size() const {
    return wave_function.size();
}

#endif
