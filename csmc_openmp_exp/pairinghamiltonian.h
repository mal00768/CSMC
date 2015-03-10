/** \file pairinghamiltonian.h
    This class describes the pairing hamiltonian defined as
    \f[H_{pairing} = H_{1b} + H_{2b} = 2\sum_{i>0}^{\omega }\epsilon_{i}(a^{\dagger}_{i}a_{i} + a^{\dagger}_{\tilde{i}}a_{\tilde{i}}) - \sum_{i,j>0}^{\omega}G_{ij}a^{\dagger}_{j}a^{\dagger}_{\tilde{j}}a_{i}a_{\tilde{i}}\f]
 
    The class is templated on the unsigned integer data type chosen to represent a state in the configuration space. This can be uint32_t, uint64_t, or uint128_t depending on the architecture
    and compiler support of the host system.
 
    The constructor requires as input the number of levels (levels), the number of pairs (NP), an std::vector containing the set of single particle energies \f$\epsilon_{k}\f$ (Esp), and lastly
    a two dimensional std::vector containing the set of pairing matrix elements \f$G_{ij}\f$ (Gij).

*/


#ifndef PAIRINGHAMILTONIAN_H
#define PAIRINGHAMILTONIAN_H

#include <vector>
#include <algorithm>
#include "wavefunction.h"

/**
    Hamiltonian class. walker type must be an unsigned type like uint32_t or uint64_t!
 
    @param NP integer value for the number of pairs.
    @param levels integer value for the number of levels.
    @param Esp std::vector containing the single particle energies.
    @param Gij 2-dimensional std::vector containing all the two body matrix elements Gij.
    @param nbits integer value for the number of bits in the walker data type.
    @param Eshift double value for the shift in the diagonal matrix elements.s
 
    The constructor requires as input the number of levels (levels), the number of pairs (NP), an std::vector containing the set of single particle energies \f$\epsilon_{k}\f$ (Esp), and lastly
    a two dimensional std::vector containing the set of pairing matrix elements \f$G_{ij}\f$ (Gij).
*/
template <typename walker>
class Hamiltonian {
    
public:
    Hamiltonian();
    Hamiltonian(unsigned, unsigned, const std::vector<double>&, const std::vector< std::vector<double> >&);
    ~Hamiltonian();
    double Hsp(const std::vector<walker>&) const; // action of H1b on a state
    double HGij(const std::vector<walker>&) const; // action of H2b on a state
    double Diagmatrixelement(const std::vector<walker>&) const;
    double Offdiagmatrixelement(unsigned, unsigned) const;
    double Probdiagonal(const std::vector<walker>&) const;
    double Proboffdiagonal(const std::vector<walker>&) const;
    double shift() const;
    
protected:
    unsigned levels; // number of levels
    unsigned NP; // number of pairs
    
    std::vector<double> Esp; // single particle energies
    std::vector< std::vector<double> > Gij; // matrix of G values
    
private:
    void determineshift();
    
    size_t nbits; // This is the number of bits we are using. The number of bits in the walker type!
    double Eshift; // diagonal shift
    
};

template <typename walker>
Hamiltonian<walker>::Hamiltonian() { }

/**
    Constructor for the Hamiltonian class.
 
    @param nlevels integer value representing the number of levels.
    @param np integer value representing the number of pairs.
    @param rEsp const reference to an std::vector containing the single particle energies.
    @param rGij const reference to a 2-dimensional std::vector containing all the two body matrix elements Gij.
*/
template <typename walker>
Hamiltonian<walker>::Hamiltonian(unsigned nlevels, unsigned np, const std::vector<double>& rEsp, const std::vector< std::vector<double> >& rGij)
: levels(nlevels), NP(np), Esp(rEsp), Gij(rGij) {
    
    determineshift(); //determine the diagonal shift
    walker tmp;
    nbits = 8 * sizeof(tmp); //set the number of bits
    
}

template <typename walker>
Hamiltonian<walker>::~Hamiltonian() { }

/**
    This function returns the result of the operation of the one-body part of the hamiltonian on a state in configuration space.
    \f[H_{1b}|\textbf{n} \rangle = 2\sum_{i>0}^{\omega }\epsilon_{i}(a^{\dagger}_{i}a_{i} + a^{\dagger}_{\tilde{i}}a_{\tilde{i}}) |\textbf{n} \rangle \f]
 
    @param rstate const reference to an std::vector representing a state in configuration space.
*/
template <typename walker>
double Hamiltonian<walker>::Hsp(const std::vector<walker>& rstate) const {
    
    double sum = 0;

    if(rstate.size() == 1) { //we don't need to iterate over rstate
        
        for(size_t i = 0; i < levels; ++i) {
            walker bit = rstate[0] & (walker(1) << i);
            if(bit != 0) { sum += Esp[i]; }
        }

    }
    else { //need to iterate!
        
        for(size_t i = 0; i < rstate.size(); ++i) {
            
            if(i != rstate.size() - 1) { //Not at the end so need to iterate over all the bits in this index
                
                for(size_t j = 0; j < nbits; ++j) {
                    walker bit = rstate[i] & (walker(1) << j);
                    if(bit != 0) { sum += Esp[i * nbits + j]; } // CHECK!
                }
                
            }
            else{ //At the end. We might not need to iterate over all the bits on the last index.
                
                for(size_t j = 0; j < (levels - nbits * (rstate.size() - 1)); ++j) {
                    walker bit = rstate[i] & (walker(1) << j);
                    if(bit != 0) { sum += Esp[i * nbits + j]; } // CHECK!
                }
                
            }
        }

    }
    
    return 2.0 * sum;
    
}

/**
    This function returns the result of the operation of the two-body part of the hamiltonian on a state in configuration space.
    \f[H_{2b}|\textbf{n} \rangle = \sum_{i,j>0}^{\omega}G_{ij}a^{\dagger}_{j}a^{\dagger}_{\tilde{j}}a_{i}a_{\tilde{i}} |\textbf{n} \rangle \f]
 
    @param rstate const reference to an std::vector representing a state in configuration space.
*/
template <typename walker>
double Hamiltonian<walker>::HGij(const std::vector<walker>& rstate) const {
    
    double sum = 0;
    
    if(rstate.size() == 1) { //we don't need to iterate over rstate
        
        for(size_t i = 0; i < levels; ++i) {
            walker bit1 = rstate[0] & (walker(1) << i);
            if(bit1 != 0) {
                for(size_t j = 0; j < levels; ++j){
                    walker bit2 = rstate[0] & (walker(1) << j);
                    if(bit2 == 0) { sum += Gij[i][j]; }
                    else if(i == j) { sum += Gij[i][j]; }
                }
            }
        }
        
    }
    else { //need to iterate!
        
        for(size_t i = 0; i < rstate.size(); ++i) {
            
            if(i != rstate.size() - 1) { //Not at the end so need to iterate over all the bits in this index
                
                for(size_t j = 0; j < nbits; ++j) {
                    walker bit1 = rstate[i] & (walker(1) << j);
                    if(bit1 != 0) {
                        for(size_t k = 0; k < rstate.size(); ++k) {
                            if(k != rstate.size() - 1) { //Not at the end so need to iterate over all the bits in this index
                                for(size_t l = 0; l < nbits; ++l) {
                                    walker bit2 = rstate[k] & (walker(1) << l);
                                    if(bit2 == 0) { sum += Gij[i * nbits + j][k * nbits + l]; }
                                    else if(i == k && j == l) { sum += Gij[i * nbits + j][k * nbits + l]; }
                                }
                            }
                            else {
                                for(size_t l = 0; l < (levels - nbits * (rstate.size() - 1)); ++l) {
                                    walker bit2 = rstate[k] & (walker(1) << l);
                                    if(bit2 == 0) { sum += Gij[i * nbits + j][k * nbits + l]; }
                                    else if(i == k && j == l) { sum += Gij[i * nbits + j][k * nbits + l]; }
                                }
                            }
                        }
                    }
                }
                
            }
            else{ //At the end. We might not need to iterate over all the bits on the last index.
                
                for(size_t j = 0; j < (levels - nbits * (rstate.size() - 1)); ++j) {
                    walker bit1 = rstate[i] & (walker(1) << j);
                    if(bit1 != 0) {
                        for(size_t k = 0; k < rstate.size(); ++k) {
                            if(k != rstate.size() - 1) { //Not at the end so need to iterate over all the bits in this index
                                for(size_t l = 0; l < nbits; ++l) {
                                    walker bit2 = rstate[k] & (walker(1) << l);
                                    if(bit2 == 0) { sum += Gij[i * nbits + j][k * nbits + l]; }
                                    else if(i == k && j == l) { sum += Gij[i * nbits + j][k * nbits + l]; }
                                }
                            }
                            else{
                                for(size_t l = 0; l < (levels - nbits * (rstate.size() - 1)); ++l){
                                    walker bit2 = rstate[k] & (walker(1) << l);
                                    if(bit2 == 0) { sum += Gij[i * nbits + j][k * nbits + l]; }
                                    else if(i == k && j == l) { sum += Gij[i * nbits + j][k * nbits + l]; }
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

/**
    This function determines Eshift. This is the diagonal shift necessary to maintain a positive sign during the random walk of each state.
*/
template <typename walker>
void Hamiltonian<walker>::determineshift() { //determine the diagonal shift
    
    Eshift = 0;
    
    if(NP == levels) {
        
        for(size_t i = 0; i < levels; ++i) { //add all the largest energy levels into the shift
            Eshift += -2.0 * Esp[i];
        }
        
        for(size_t i = 0; i < levels; ++i) { //add in the contribution for G transfers along the diagonal
            Eshift -= Gij[i][i];
        }
        
    }
    else{
        
        std::vector<double> Esptmp = Esp; //make a temp copy of the single particle energies
        std::sort(Esptmp.begin(), Esptmp.end()); //sort the temp copy
        
        std::vector<double> Gijtmp; //now make a copy of all the diagonal G
        for(size_t i = 0; i < Gij.size(); ++i) {
            Gijtmp.push_back(Gij[i][i]);
        }
        std::sort(Gijtmp.begin(), Gijtmp.end()); //sort the temp copy
        
        for(size_t i = Esptmp.size() - 1; i >= Esptmp.size() - NP; --i) { //add all the largest energy levels into the shift
            Eshift += -2.0*Esptmp[i];
        }
        
        for(size_t i = Gijtmp.size() - 1; i >= Gijtmp.size() - NP; --i) { //add in the contribution for G transfers along the diagonal
            Eshift -= Gijtmp[i];
        }
        
    }
    
}

/**
    This function returns a diagonal matrix element of the hamiltonian matrix i.e. \f$\langle \textbf{n}_{i} |H| \textbf{n}_{i} \rangle\f$.
    @param rstate const reference to an std::vector representing a state in configuration space.
*/
template <typename walker>
double Hamiltonian<walker>::Diagmatrixelement(const std::vector<walker>& rstate) const {
    
    double sum = -Hsp(rstate) - Eshift;
    
    if(rstate.size() == 1) { //we don't need to iterate over rstate
        
        for(size_t i = 0; i < levels; ++i){
            walker bit = rstate[0] & (walker(1) << i);
            if(bit != 0) { sum += Gij[i][i]; }
        }
        
    }
    else { //need to iterate!
        
        for(size_t i = 0; i < rstate.size(); ++i) {
            
            if(i != rstate.size() - 1) { //Not at the end so need to iterate over all the bits in this index
                
                for(size_t j = 0; j < nbits; ++j){
                    walker bit = rstate[i] & (walker(1) << j);
                    if(bit != 0) { sum += Gij[i*nbits + j][i*nbits + j]; } // CHECK!
                }
                
            }
            else { //At the end. We might not need to iterate over all the bits on the last index.
                
                for(size_t j = 0; j < (levels - nbits * (rstate.size() - 1)); ++j) {
                    walker bit = rstate[i] & (walker(1) << j);
                    if(bit != 0) { sum += Gij[i * nbits + j][i * nbits + j]; } // CHECK!
                }
                
            }
        }
        
    }

    return sum;
    
}

/**
    This function returns an off-diagonal matrix element of the hamiltonian matrix i.e. \f$\langle \textbf{n}_{j} |H| \textbf{n}_{i} \rangle\f$.
    @param rstate const reference to an std::vector representing a state in configuration space.
*/
template <typename walker>
double Hamiltonian<walker>::Offdiagmatrixelement(unsigned i, unsigned j) const {
    return Gij[i][j];
}

/**
    This function returns the probability of making a diagonal transition in configuration space i.e. where a state transitions into the same state and no pairs are moved.
 
    In general the set of transition probabilites are proportional to the normalized matrix elements of all states connected to the current configuration.
    \f[\mathcal{P}_{ij} = \frac{\langle \textbf{n}_{j} |H| \textbf{n}_{i} \rangle}{d_{ij}}\f]
    where \f$d_{ij} = \sum_{j} \langle \textbf{n}_{j} |H| \textbf{n}_{i} \rangle\f$ following the linear normalization condition.
 
    @param rstate const reference to an std::vector representing a state in configuration space.
*/
template <typename walker>
double Hamiltonian<walker>::Probdiagonal(const std::vector<walker>& rstate) const {    
    return abs(Diagmatrixelement(rstate)) / (-Hsp(rstate) - Eshift + HGij(rstate));
}

/**
    This function returns the probability of making a off-diagonal transition in configuration space i.e. where a state transitions into a new state and
    a pair is moved from an occupied level to an unoccupied level.
 
    In general the set of transition probabilites are proportional to the normalized matrix elements of all states connected to the current configuration.
    \f[\mathcal{P}_{ij} = \frac{\langle \textbf{n}_{j} |H| \textbf{n}_{i} \rangle}{d_{ij}}\f]
    where \f$d_{ij} = \sum_{j} \langle \textbf{n}_{j} |H| \textbf{n}_{i} \rangle\f$ following the linear normalization condition.
 
    @param rstate const reference to an std::vector representing a state in configuration space.
*/
template <typename walker>
double Hamiltonian<walker>::Proboffdiagonal(const std::vector<walker>& rstate) const {
    return (1 - Probdiagonal(rstate)) / double(NP * (levels - NP));
}

/**
    This function returns the diagonal shift Eshift.
*/
template <typename walker>
double Hamiltonian<walker>::shift() const {
    return Eshift;
}

#endif
