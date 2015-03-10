// File: am.h
// Author: Mark Lingle
//
// This is based on the implementations:
//  http://ideone.com/gzFzI
//  http://www.keithschwarz.com/interesting/code/?dir=alias-method
//
// Walker's alias method implemented using Vose's algorithm.
// This method allows for efficient sampling of random values from a
// discrete probability distribution in O(1) time each after O(n) preprocessing time.
//

#ifndef _AM_H_
#define _AM_H_
#include <vector>
#include "randomc.h"

class Alias_Method {
    
private:
    // The probability and alias tables.
    std::vector<size_t> alias_;
    std::vector<double> probability_;
    
public:
    Alias_Method(const std::vector<double>&);
    size_t next(std::vector<CRandomMersenne>&);
    
};

Alias_Method::Alias_Method(const std::vector<double>& probability) : probability_(probability) {
    
    // Allocate space for the alias table.
    alias_.resize(probability_.size());
    
    // Compute the average probability and store it for later.
    const double average = 1.0 / double(probability_.size());
    
    // two worklists used to populate the tables
    std::vector<size_t> small, large;
    
    // Populate the stacks with the input probabilities.
    for(size_t i = 0; i < probability_.size(); ++i) {
        // If the probability is below the average probability, then we add
        // it to the small list; otherwise we add it to the large list.
        if(probability_[i] >= average) { large.push_back(i); }
        else { small.push_back(i); }
    }
    
    // As a note: in the mathematical specification of the algorithm, we
    // will always exhaust the small list before the big list.  However,
    // due to floating point inaccuracies, this is not necessarily true.
    // Consequently, this inner loop (which tries to pair small and large
    // elements) will have to check that both lists aren't empty.
    while(!small.empty() && !large.empty()) {
        // Get the index of the small and the large probabilities.
        size_t less = small.back();
        small.pop_back();
        size_t more = large.back();
        large.pop_back();
        
        alias_[less] = more;
        
        // Decrease the probability of the larger one by the appropriate
        // amount.
        probability_[more] = (probability_[more] + probability_[less]) - average;
        
        // If the new probability is less than the average, add it into the
        // small list; otherwise add it to the large list.
        if (probability_[more] >= average) { large.push_back(more); }
        else {small.push_back(more);}
    }
    
    // At this point, everything is in one list, which means that the
    // remaining probabilities should all be 1/n.  Based on this, set them
    // appropriately.  Due to numerical issues, we can't be sure which
    // stack will hold the entries, so we empty both.
    while(!small.empty()) {
        probability_[small.back()] = average;
        small.pop_back();
    }
    
    while(!large.empty()) {
        probability_[large.back()] = average;
        large.pop_back();
    }
    
    // These probabilities have not yet been scaled up to be such that
    // 1/n is given weight 1.0.  We do this here.
    size_t n = static_cast<size_t> (probability_.size());
    for(size_t i = 0; i < probability_.size(); i++) { probability_[i] *= double(n); }
    
}

/*
 Return a random value from the underlying distribution.
*/
size_t Alias_Method::next(std::vector<CRandomMersenne>& rgenerators) {
    // Generate a fair die roll to determine which column to inspect.
    size_t column = rgenerators[ omp_get_thread_num() ].IRandomX(0, probability_.size() - 1);
    
    // Generate a biased coin toss to determine which option to pick.
    bool coin_toss = rgenerators[ omp_get_thread_num() ].Random() < probability_[column];
    
    // Based on the outcome, return either the column or its alias.
    return coin_toss ? column : alias_[column];
}

#endif
