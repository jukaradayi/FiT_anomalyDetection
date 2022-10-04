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

namespace StreamGraphs {
class SortedCounters {
    public:
    std::vector<uint64_t> values; ///< the counter values, stored in a decreasing order
    //vector<uint64_t> counter2pos; ///< store position of counter c in values
   
    std::unordered_map<uint64_t, uint64_t> counter2pos; ///< store position of counter c in values

    std::vector<uint64_t> pos2counter; ///< store name of counter at position i in values
    //vector<uint64_t> val2pos; ///< store the first index of a counter at given value
    std::unordered_map<uint64_t, uint64_t> val2pos; ///< store the first index of a counter at given value

    //vector<uint64_t> distrib; ///< the counter distribution, store the number of counter with value v
    std::unordered_map<uint64_t, uint64_t> distrib;
    SortedCounters();
    ~SortedCounters();
    uint64_t get_value(uint64_t counter);
    bool has_counter(uint64_t counter);

    void increase_counter(uint64_t counter);

    void decrease_counter(uint64_t counter);

    void add_counter(uint64_t counter);
    void remove_counter(uint64_t counter);

};

}
#endif
