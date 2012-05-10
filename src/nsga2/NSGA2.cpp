#include "nsga2/NSGA2.h"

#include <iostream>


using namespace nsga2;

NSGA2::NSGA2() :
    nreal(0),
    nbin(0),
    nobj(0),
    ncon(0),
    popsize(0),
    pcross_real(0),
    pcross_bin(0),
    pmut_real(0),
    pmut_bin(0),
    eta_c(0),
    eta_m(0),
    ngen(0),
    nbinmut(0),
    nrealmut(0),
    nbincross(0),
    nrealcross(0),
    nbits(0),
    limits_realvar(0),
    limits_binvar(0),
    bitlength(0),
    choice(0),
    obj1(0),
    obj2(0),
    obj3(0),
    angle1(0),
    angle2(0),
    parent_pop(0),
    child_pop(0),
    mixed_pop(0) {
}

NSGA2::~NSGA2() {
    if (parent_pop) {
        delete parent_pop;
        parent_pop = 0;
    }
    if (child_pop) {
        delete child_pop;
        child_pop = 0;
    }
    if (mixed_pop) {
        delete mixed_pop;
        mixed_pop = 0;
    }
}
