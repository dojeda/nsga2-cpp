#include "nsga2/NSGA2.h"

#include <iostream>

using namespace nsga2;
using namespace std;

NSGA2::NSGA2() :
    nreal(-1),
    nbin(-1),
    nobj(-1),
    ncon(-1),
    popsize(-1),
    ngen(-1),
    pcross_real(-1),
    pcross_bin(-1),
    pmut_real(-1),
    pmut_bin(-1),
    eta_c(-1),
    eta_m(-1),
    nbits(0),
    limits_realvar(0),
    limits_binvar(0),
    // choice(0),
    // obj1(0),
    // obj2(0),
    // obj3(0),
    // angle1(0),
    // angle2(0),
    nbinmut(0),
    nrealmut(0),
    nbincross(0),
    nrealcross(0),
    bitlength(0),
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


void NSGA2::initialize() throw (nsga2exception) {

    cout << "Initializing NSGA-II. Checking configuration" << endl;
    
    if (nreal < 0)
        throw nsga2exception("Invalid number of real variables");
    if (nbin < 0)
        throw nsga2exception("Invalid number of binary variables");
    if (nreal==0 && nbin==0)
        throw nsga2exception("Zero real and binary variables");
    if (nobj < 1)
        throw nsga2exception("Invalid number of objective functions");
    if (ncon < 0) 
        throw nsga2exception("Invalid number of constraints");
    if (popsize<4 || (popsize%4)!= 0)
        throw nsga2exception("Invalid size of population");
    if (pcross_real<0.0 || pcross_real>1.0)
        throw nsga2exception("Invalid probability of real crossover");
    if (pmut_real<0.0 || pmut_real>1.0)
        throw nsga2exception("Invalid probability of real mutation");
    if (pcross_bin<0.0 || pcross_bin>1.0)
        throw nsga2exception("Invalid probability of binary crossover");
    if (pmut_bin<0.0 || pmut_bin>1.0)
        throw nsga2exception("Invalid probability of binary mutation");
    if (eta_c<=0)
        throw nsga2exception("Invalid distribution index for crossover");
    if (eta_m<=0)
        throw nsga2exception("Invalid distribution index for mutation");
    if (ngen < 1)
        throw nsga2exception("Invalid number of generations");
    if (nbin != 0 && nbits.size() == 0)
        throw nsga2exception("Invalid number of bits for binary variables");
    if (limits_realvar.size() != nreal)
        throw nsga2exception("Invalid number of real variable limits");
    if (limits_binvar.size() != nbin)
        throw nsga2exception("Invalid number of binary variable limits");

}
