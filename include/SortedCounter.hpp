/**
 * Maintain decreasing sorted counter with 
 * increments in O(1).
 * 
 *
 */

#ifndef STREAMGRAPHS_SORTED_COUNTERS_HPP_
#define STREAMGRAPHS_SORTED_COUNTERS_HPP_

#include <vector>
#include <unordered_map>
#include "globals.hpp"
namespace StreamGraphs {
template <class T>
class SortedCounters {
    public:
    std::vector<uint64_t> values; ///< the counter values, stored in a decreasing order
    //vector<uint64_t> counter2pos; ///< store position of counter c in values
   
    //std::unordered_map<uint64_t, uint64_t> counter2pos; ///< store position of counter c in values
    std::unordered_map<T, uint64_t> counter2pos; ///< store position of counter c in values


    std::vector<T> pos2counter; ///< store name of counter at position i in values
    //vector<uint64_t> val2pos; ///< store the first index of a counter at given value
    std::unordered_map<uint64_t, uint64_t> val2pos; ///< store the first index of a counter at given value

    //vector<uint64_t> distrib; ///< the counter distribution, store the number of counter with value v
    std::unordered_map<uint64_t, uint64_t> distrib;
    SortedCounters();
    ~SortedCounters();
    uint64_t get_value(T counter);
    bool has_counter(T counter);

    void increase_counter(T counter);

    void decrease_counter(T counter);

    void add_counter(T counter);
    void remove_counter(T counter);

};

/* Sorted Counters implementation */
template<class T>
SortedCounters<T>::SortedCounters(){}

template<class T>
SortedCounters<T>::~SortedCounters(){}

template<class T>
uint64_t SortedCounters<T>::get_value(T counter) {
    if (counter2pos.find(counter) != counter2pos.end()) {
        uint64_t counter_idx = counter2pos[counter];
        return values[counter_idx];
    } else {
        return none;
    }
}

template<class T>
void SortedCounters<T>::add_counter(T counter) {

    if (distrib.find(0) == distrib.end()) {
        distrib[0] = 0;
        val2pos[0] = values.size();
    }
    distrib[0] += 1;
    counter2pos[counter] = values.size();
    pos2counter.push_back(counter);
    values.push_back(0);
}

template<class T>
void SortedCounters<T>::increase_counter(T counter) {
    uint64_t val = get_value(counter);
    uint64_t val_idx = val2pos[val];
    T first_counter = pos2counter[val_idx];

    // swap counters
    if (first_counter != counter) {

        uint64_t counter_idx = counter2pos[counter];
        uint64_t first_counter_idx = counter2pos[first_counter];

        counter2pos[counter] = first_counter_idx;
        counter2pos[first_counter] = counter_idx;
        pos2counter[counter_idx] = first_counter;
        pos2counter[first_counter_idx] = counter;

    }

    // update distribution
    if (distrib.find(val+1) == distrib.end()) {
        distrib[val + 1] = 0;
        val2pos[val + 1] = val2pos[val];
    }
    distrib[val + 1] += 1;
    distrib[val] -= 1;
    val2pos[val] += 1;

    if (distrib[val] == 0) {
        // check removal
        distrib.erase(val);
        val2pos.erase(val);
    }

    values[counter2pos[counter]] += 1;
}

template<class T>
void SortedCounters<T>::decrease_counter(T counter) {
    uint64_t val = get_value(counter);
    uint64_t val_idx = val2pos[val];
    uint64_t N_val = distrib[val] - 1;
    T last_counter = pos2counter[val_idx + N_val];

    // swap counters
    if (last_counter != counter) { // TODO make function

        uint64_t counter_idx = counter2pos[counter];
        uint64_t last_counter_idx = counter2pos[last_counter];
        counter2pos[counter] = last_counter_idx;
        counter2pos[last_counter] = counter_idx;
        pos2counter[counter_idx] = last_counter;
        pos2counter[last_counter_idx] = counter;

    }

    // update distribution
    if (distrib.find(val-1) == distrib.end()) {// TODO 
        distrib[val - 1] = 0;
        val2pos[val - 1] = val2pos[val] + distrib[val];
    }
    distrib[val - 1] += 1;
    distrib[val] -= 1;
    val2pos[val - 1] -= 1;

    if (distrib[val] == 0) {
        distrib.erase(val);
        val2pos.erase(val);

    }
    values[counter2pos[counter]] -= 1;

}

template<class T>
void SortedCounters<T>::remove_counter(T counter){

    distrib[0] -= 1;
    if (distrib[0] == 0) {
        distrib.erase(0);
        val2pos.erase(0);
    }
    uint64_t val = get_value(counter);
    T last_counter = pos2counter[values.size() - 1];

    // swap counters
    if (last_counter != counter) { // TODO make function
        uint64_t counter_idx = counter2pos[counter];
        uint64_t last_counter_idx = counter2pos[last_counter];
        counter2pos[counter] = last_counter_idx;
        counter2pos[last_counter] = counter_idx;
        pos2counter[counter_idx] = last_counter;
        pos2counter[last_counter_idx] = counter;

    }
    counter2pos.erase(counter);
    pos2counter.pop_back();
    values.pop_back();
}



} /* namespace Streamgraphs */
#endif
