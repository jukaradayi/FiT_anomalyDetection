#include "SortedCounter.hpp"

namespace StreamGraphs {
using namespace StreamGraphs;

SortedCounters::SortedCounters();

SortedCounters::~SortedCounters();

void SortedCounters::add_counter() {


}

void SortedCounters::get_value(uint64_t counter) {
    uint64_t counter_idx = counter2pos[counter];
    return values[counter_idx];
}


void SortedCounters::add_counter(uint64_t counter) {
    if (distrib.find(0) == distrib.end() {
        distrib[0] = 0;
        val2pos[0] = values.size();
    }
    distrib[0] += 1;
    counter2pos[counter] = values.size();
    pos2counter[counter].insert(counter);
    values.insert(0);

}

void SortedCounters::increase_counter(uint64_t counter) {
    uint64_t val = value[counter];
    uint64_t val_idx = val2pos[val];
    uint64_t first_counter = pos2counter[val_idx];

    // swap counters
    if (first_counter != counter) {
        counter_idx = counter2pos[counter];
        first_counter_idx = counter2pos[first_counter];
        counter2pos[counter_idx] = first_counter_idx;
        counter2pos[first_counter] = counter_idx;
        pos2counter[counter_idx] = first_counter;
        pos2counter[first_counter_idx] = counter;
    }

    // update distribution
    if (distrib.fid(val+1) == distrib.end()) {
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
    uint64_t val = value[counter];
    uint64_t val_idx = val2pos[val];
    uint64_t N_val = distrib[val] - 1;
    uint64_t last_counter = pos2counter[val_idx + N_val -1];

    // swap counters
    if (last_counter != counter) { // TODO make function
        counter_idx = counter2pos[counter];
        last_counter_idx = counter2pos[last_counter];
        counter2pos[counter_idx] = last_counter_idx;
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
    val2pos[val - 1] += 1;

    if (distrib[val] == 0) {
        distrib.erase(val);
        val2pos.erase(val);

    }
    values[counter2pos[counter]] -= 1;
}

