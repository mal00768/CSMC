/** \file operators.h
    This class defines the functions necessary to calculate the ground state energy expectation values, \f$\sum_{i} \langle \textbf{n}_{i} | H | \Psi \rangle\f$ and
    \f$\langle \Psi | H | \Psi \rangle\f$, and the occupation numbers for all levels. It inherits public and protected members from the Hamiltonian<walker> class. This 
    is necessary to include functions defining how the hamiltonian acts on states.
 
    The class is templated on the unsigned integer data type chosen to represent a state in the configuration space. This can be uint32_t, uint64_t, or uint128_t depending on the architecture
    and compiler support of the host system.
    
    The constructor takes as input the number of levels (levels), the number of pairs (np), an std::vector containing the set of single particle energies \f$\epsilon_{k}\f$ (Esp), a
    two dimensional std::vector containing the set of pairing matrix elements \f$G_{ij}\f$ (Gij), and finally a const reference to the wave function (Wavefunction<walker>).
*/


#ifndef OPERATORS_H
#define OPERATORS_H

#include <unordered_map>
#include "pairinghamiltonian.h"
#include "wavefunction.h"

/**
    Operator class. Subclasses from Hamiltonian class. walker type must be an unsigned type like uint32_t or uint64_t!
 
    @param nbits unsigned integer indicating the number of bits in the walker type.
    @param levels unsigned integer containing the number of levels.
    @param np unsigned integer containing the number of pairs.
    @param wf const reference to the wave function class.
*/
template <typename walker>
class Operators : protected Hamiltonian<walker> {
    
public:
    Operators(unsigned, unsigned, const std::vector<double>&, const std::vector< std::vector<double> >&, const Wavefunction<walker>&);
    ~Operators();
    double Egs() const;
    double ExpectationH() const;
    double ExpHsp(typename std::unordered_map< std::vector<walker>, double >::const_iterator) const;
    double ExpHGij(typename std::unordered_map< std::vector<walker>, double >::const_iterator) const;
    double occupation_level(unsigned);
    void occupations(std::vector<double>&);
    
private:
    size_t nbits; // This is the number of bits we are using. The number of bits in the walker type!
    unsigned levels;
    unsigned np;
    const Wavefunction<walker>& wf;
    using Hamiltonian<walker>::Hsp;
    using Hamiltonian<walker>::HGij;
    
};

/**
 Constructor for the Operators class.
 
 @param nlevels integer value representing the number of levels.
 @param nnp integer value representing the number of pairs.
 @param rEsp const reference to an std::vector containing the single particle energies.
 @param rGij const reference to a 2-dimensional std::vector containing all the two body matrix elements Gij.
 @param rwf const reference to the wave function class.
 */
template <typename walker>
Operators<walker>::Operators(unsigned nlevels, unsigned nnp, const std::vector<double>& rEsp, const std::vector< std::vector<double> >& rGij, const Wavefunction<walker>& rwf)
: Hamiltonian<walker>(nlevels, nnp, rEsp, rGij), levels(nlevels), np(nnp), wf(rwf) {
    
    walker tmp;
    nbits = 8 * sizeof(tmp); // set the number of bits
    
}

template <typename walker>
Operators<walker>::~Operators() {}

/**
    This function computes the occupation number for the level k. The wave function should be normalized using the square norm first.
 
    @param k unsigned integer indicating the desired level
*/
template <typename walker>
double Operators<walker>::occupation_level(unsigned k){
 
    unsigned index = unsigned(double(k) / double(nbits));
    double O = 0;
    
    for(auto itr = wf.begin(); itr != wf.end(); ++itr) {
        unsigned kk = k;
        if(k >= nbits) { kk -= nbits * ( (itr -> first).size() - 1 ); }
        walker bit = (itr -> first)[index] & (walker(1) << kk);
        if(bit != 0){ //if a level is occupied then add the amplitude
            O += (itr -> second) * (itr -> second);
        }
    }
    
    return 2.0 * O;

}

/**
    This function calls the occupation_level() function for each level andsaves the result into the std::vector occupationn.
 
    @param occupationn reference to an std::vector used to save the occupation numbers for each level.
*/
template <typename walker>
void Operators<walker>::occupations(std::vector<double>& occupationn) {
    
    const_cast<Wavefunction<walker>&> (wf).normalise_sq();
    
#pragma omp parallel for
    for(size_t l = 0; l < occupationn.size(); ++l) {
        occupationn[l] = occupation_level(l);
    }
    
}

/**
    This function returns the ground state energy through the expectation value \f$\sum_{i} \langle \textbf{n}_{i} | H | \Psi \rangle\f$.
    This uses the non-standard linear normalization of the wave function. It does not require the wave function to be normalized prior
    to the function call.
*/
template <typename walker>
double Operators<walker>::Egs() const {
    
    const_cast<Wavefunction<walker>&> (wf).normalise_ln(); // ensure proper normalization
    
    double Egs = 0;
    for(typename Wavefunction<walker>::const_iterator it = wf.begin(); it != wf.end(); ++it) {
        Egs += (it->second) * (Hsp(it->first) - HGij(it->first));
    }
    
    return Egs;
    
}

/**
    This function returns the ground state energy by calculating \f$\langle \Psi | H | \Psi \rangle\f$. This uses the quantum mechanical square 
    normalization. It does not require the wave function to be normalized prior to the function call.
*/
template <typename walker>
double Operators<walker>::ExpectationH() const {
    
    const_cast<Wavefunction<walker>&> (wf).normalise_sq(); // ensure proper normalization
    
    double sum = 0;
    for(typename Wavefunction<walker>::const_iterator it = wf.begin(); it != wf.end(); ++it) {
        sum += (it->second) * (it->second);
    }
    
    double E = 0;
    for(typename Wavefunction<walker>::const_iterator it = wf.begin(); it != wf.end(); ++it) {
        E += (ExpHsp(it) - ExpHGij(it));
    }
    
    return E / sum;
}

/**
    This function returns the expectation value of the one body part of the hamiltonian in configuration space. 
    \f[ \alpha_{i}^{2} \langle \textbf{n}_{i}|H_{1b}|\textbf{n}_{i} \rangle \f]
    \f[ = \alpha_{i}^{2} \langle \textbf{n}_{i}| 2\sum_{i>0}^{\omega }\epsilon_{i}(a^{\dagger}_{i}a_{i} + a^{\dagger}_{\tilde{i}}a_{\tilde{i}}) |\textbf{n}_{i} \rangle \f]
 
    @param it const iterator pointing to a state in the wave function.
*/
template <typename walker>
double Operators<walker>::ExpHsp(typename std::unordered_map< std::vector<walker>, double >::const_iterator it) const {
    
    const std::vector<walker>& rstate = it -> first;
    double sum = 0;
    
    if(rstate.size() == 1) { // we don't need to iterate over rstate
        
        for(size_t i = 0; i < levels; ++i) {
            walker bit = rstate[0] & (walker(1) << i);
            if(bit != 0) { sum += this -> Esp[i]; }
        }
        
    }
    else { // need to iterate!
        
        for(size_t i = 0; i < rstate.size(); ++i) {
            
            if(i != rstate.size() - 1) { // Not at the end so need to iterate over all the bits in this index.
                
                for(size_t j = 0; j < nbits; ++j) {
                    walker bit = rstate[i] & (walker(1) << j);
                    if(bit != 0) { sum += this -> Esp[i * nbits + j]; }
                }
                
            }
            else { //A t the end. We might not need to iterate over all the bits on the last index.
                
                for(size_t j = 0; j < (levels - nbits * (rstate.size() - 1)); ++j) {
                    walker bit = rstate[i] & (walker(1) << j);
                    if(bit != 0) { sum += this -> Esp[i * nbits + j]; }
                }
                
            }
        }
        
    }
    
    return 2.0 * pow(it -> second, 2.0) * sum;
    
}

/**
    This function returns the expectation value of the two body part of the hamiltonian in configuration space.
    \f[ \alpha_{i}\alpha_{j} \langle \textbf{n}_{j}|H_{2b}|\textbf{n}_{i} \rangle \f]
    \f[ = \alpha_{i}\alpha_{j} \langle \textbf{n}_{j}| \sum_{i,j>0}^{\omega}G_{ij}a^{\dagger}_{j}a^{\dagger}_{\tilde{j}}a_{i}a_{\tilde{i}} |\textbf{n}_{i} \rangle \f]
 
    @param it const iterator pointing to a state in the wave function.
*/
template <typename walker>
double Operators<walker>::ExpHGij(typename std::unordered_map< std::vector<walker>, double >::const_iterator it) const {
    
    const std::vector<walker>& rstate = it -> first;
    double sum = 0;
    
    if(rstate.size() == 1) { // we don't need to iterate over rstate
        
        for(size_t i = 0; i < levels; ++i) {
            
            std::vector<walker> tmpstate = it -> first;
            walker bit1 = rstate[0] & (walker(1) << i); // should only procede into loop if level is occupied
            
            if(bit1 != 0) {
                for(size_t j = 0; j < levels; j++) {
                    walker bit2 = tmpstate[0] & (walker(1) << j);
                    
                    if(i != j && bit2 == 0) {
                        tmpstate[0] ^= walker(1) << i;
                        tmpstate[0] ^= walker(1) << j;
                        sum += (it -> second) * (wf.configuration_weight(tmpstate)) * (this -> Gij[i][j]);
                    }
                    else if(i == j){
                        sum += (it->second) * (it -> second) * (this -> Gij[i][j]);
                    }
                    
                    tmpstate = rstate; // reset the tmpstate
                }
            }
            
        }
        
    }
    else { // need to iterate!
        
        for(size_t i = 0; i < rstate.size(); ++i) {
            
            if(i != rstate.size() - 1){ // Not at the end so need to iterate over all the bits in this index
                
                for(size_t j = 0; j < nbits; ++j){
                    
                    std::vector<walker> tmpstate = it -> first;
                    walker bit1 = rstate[i] & (walker(1) << j);
                    
                    if(bit1 != 0) {
                        for(size_t k = 0; k < rstate.size(); ++k) {
                            if(k != rstate.size() - 1) { //Not at the end so need to iterate over all the bits in this index
                                for(size_t l = 0; l < nbits; ++l){
                                    walker bit2 = tmpstate[k] & (walker(1) << l);
                                    
                                    if( (i * nbits + j) != (k * nbits + l) && bit2 == 0) {
                                        tmpstate[i] ^= walker(1) << j;
                                        tmpstate[k] ^= walker(1) << l;
                                        sum += (it -> second) * (wf.configuration_weight(tmpstate)) * (this -> Gij[i * nbits + j][k * nbits + l]);
                                    }
                                    else if( (i * nbits + j) == (k * nbits + l) ) {
                                        sum += (it -> second) * (it -> second) * (this -> Gij[i * nbits + j][k * nbits + l]);
                                    }
                                    
                                    tmpstate = rstate; //reset the tmpstate
                                    
                                }
                            }
                            else {
                                for(size_t l = 0; l < (levels - nbits * (rstate.size() - 1)); ++l) {
                                    walker bit2 = tmpstate[k] & (walker(1) << l);
                                    
                                    if( (i * nbits + j) != (k * nbits + l) && bit2 == 0) {
                                        tmpstate[i] ^= walker(1) << j;
                                        tmpstate[k] ^= walker(1) << l;
                                        sum += (it -> second) * (wf.configuration_weight(tmpstate)) * (this -> Gij[i * nbits + j][k * nbits + l]);
                                    }
                                    else if( (i * nbits + j) == (k * nbits + l) ) {
                                        sum += (it -> second) * (it -> second) * (this -> Gij[i * nbits + j][k * nbits + l]);
                                    }
                                    
                                    tmpstate = rstate; //reset the tmpstate
                                    
                                }
                            }
                            
                        }
                    }
                }
                
            }
            else { // At the end. We might not need to iterate over all the bits on the last index.
            
                for(size_t j = 0; j < (levels - nbits * (rstate.size() - 1)); ++j) {
                    
                    std::vector<walker> tmpstate = it -> first;
                    walker bit1 = rstate[i] & (walker(1) << j);
                    
                    if(bit1 != 0) {
                        for(size_t k = 0; k < rstate.size(); ++k){
                            if(k != rstate.size() - 1){ //Not at the end so need to iterate over all the bits in this index
                                for(size_t l = 0; l < nbits; ++l) {
                                    walker bit2 = tmpstate[k] & (walker(1) << l);
                                    
                                    if( (i * nbits + j) != (k * nbits + l) && bit2 == 0) {
                                        tmpstate[i] ^= walker(1) << j;
                                        tmpstate[k] ^= walker(1) << l;
                                        sum += (it -> second) * (wf.configuration_weight(tmpstate)) * (this -> Gij[i * nbits + j][k * nbits + l]);
                                    }
                                    else if( (i * nbits + j) == (k * nbits + l) ) {
                                        sum += (it -> second) * (it -> second) * (this -> Gij[i * nbits + j][k * nbits + l]);
                                    }
                                    
                                    tmpstate = rstate; //reset the tmpstate
                                    
                                }
                            }
                            else {
                                for(size_t l = 0; l < (levels - nbits * (rstate.size() - 1)); ++l) {
                                    walker bit2 = rstate[k] & (walker(1) << l);
                                    
                                    if( (i * nbits + j) != (k * nbits + l) && bit2 == 0){
                                        tmpstate[i] ^= walker(1) << j;
                                        tmpstate[k] ^= walker(1) << l;
                                        sum += (it -> second) * (wf.configuration_weight(tmpstate)) * (this -> Gij[i * nbits + j][k * nbits + l]);
                                    }
                                    else if( (i*nbits + j) == (k * nbits + l) ) {
                                        sum += (it -> second) * (it -> second) * (this -> Gij[i * nbits + j][k * nbits + l]);
                                    }
                                    
                                    tmpstate = rstate; //reset the tmpstate
                                    
                                }
                                
                            }
                            
                        }
                    }
                }
                
            }
            
        }
        
    }
    
    return sum;
    
}


#endif
