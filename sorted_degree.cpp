//
// Created by Apple on 30/05/2022.
//

#include "sorted_degree.hpp"
namespace StreamGraphs {
    ConstantTimeMax::ConstantTimeMax() {
        this->deg_seq = new uint64_t[N] {none};
    }
}