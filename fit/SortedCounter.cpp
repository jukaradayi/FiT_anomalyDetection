#include "globals.hpp"
#include "SortedCounter.hpp"
#include <iostream>

namespace StreamGraphs {
using namespace StreamGraphs;

SortedCounters::SortedCounters(){}

SortedCounters::~SortedCounters(){}

uint64_t SortedCounters::get_value(uint64_t counter) {
    if (counter2pos.find(counter) != counter2pos.end()) {
        uint64_t counter_idx = counter2pos[counter];
        return values[counter_idx];
    } else {
        return none;
    }
}


void SortedCounters::add_counter(uint64_t counter) {

    if (distrib.find(0) == distrib.end()) {
        distrib[0] = 0;
        val2pos[0] = values.size();
    }
    distrib[0] += 1;
    counter2pos[counter] = values.size();
    pos2counter.push_back(counter);
    values.push_back(0);
}

void SortedCounters::increase_counter(uint64_t counter) {
    uint64_t val = get_value(counter);
    uint64_t val_idx = val2pos[val];
    uint64_t first_counter = pos2counter[val_idx];

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

void SortedCounters::decrease_counter(uint64_t counter) {
    uint64_t val = get_value(counter);
    uint64_t val_idx = val2pos[val];
    uint64_t N_val = distrib[val] - 1;
    uint64_t last_counter = pos2counter[val_idx + N_val];

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
void SortedCounters::remove_counter(uint64_t counter){

    distrib[0] -= 1;
    if (distrib[0] == 0) {
        distrib.erase(0);
        val2pos.erase(0);
    }
    uint64_t val = get_value(counter);
    uint64_t last_counter = pos2counter[values.size() - 1];

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
}

