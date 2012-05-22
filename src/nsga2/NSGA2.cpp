#include "nsga2/NSGA2.h"
#include "rand.h"

#include <iostream>
#include <numeric>
#include <string>

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


    init_streams();
    report_parameters(fpt5);
    
    nbinmut    = 0;
    nrealmut   = 0;
    nbincross  = 0;
    nrealcross = 0;
    bitlength  = std::accumulate(nbits.begin(), nbits.end(), 0);

    parent_pop = new population(popsize,
                                nreal,
                                nbin,
                                ncon,
                                nbits,
                                limits_realvar,
                                limits_binvar,
                                nobj);
    child_pop  = new population(popsize,
                                nreal,
                                nbin,
                                ncon,
                                nbits,
                                limits_realvar,
                                limits_binvar,
                                nobj);
    mixed_pop  = new population(popsize*2,
                                nreal,
                                nbin,
                                ncon,
                                nbits,
                                limits_realvar,
                                limits_binvar,
                                nobj);

    randomize();
    parent_pop->initialize();
    cout << "Initialization done, now performing first generation" << endl;

    parent_pop->decode();
    parent_pop->evaluate();
    //parent_pop->assign_rank_and_crowding_distance();
    // assign_rank_and_crowding_distance (parent_pop);

    report_pop(*parent_pop,fpt1);
    
}

void NSGA2::init_streams() {
    fpt1.open("nsga2_initial_pop.out" , ios::out | ios::trunc);
    fpt2.open("nsga2_final_pop.out"   , ios::out | ios::trunc);
    fpt3.open("nsga2_best_pop.out"    , ios::out | ios::trunc);
    fpt4.open("nsga2_all_pop.out"     , ios::out | ios::trunc);
    fpt5.open("nsga2_params.out"      , ios::out | ios::trunc);

    fpt1.setf(ios::scientific);
    fpt2.setf(ios::scientific);
    fpt3.setf(ios::scientific);
    fpt4.setf(ios::scientific);
    fpt5.setf(ios::scientific);

    fpt1 << "# This file contains the data of initial population\n";
    fpt2 << "# This file contains the data of final population\n";
    fpt3 << "# This file contains the data of final feasible population (if found)\n";
    fpt4 << "# This file contains the data of all generations\n";
    fpt5 << "# This file contains information about inputs as read by the program\n";

    fpt1 << "# of objectives = "    << nobj
         << ", # of constraints = " << ncon
         << ", # of real_var = "    << nreal
         << ", # of bits of bin_var = " << bitlength
         << ", constr_violation, rank, crowding_distance\n";
    fpt2 << "# of objectives = "    << nobj
         << ", # of constraints = " << ncon
         << ", # of real_var = "    << nreal
         << ", # of bits of bin_var = " << bitlength
         << ", constr_violation, rank, crowding_distance\n";
    fpt3 << "# of objectives = "    << nobj
         << ", # of constraints = " << ncon
         << ", # of real_var = "    << nreal
         << ", # of bits of bin_var = " << bitlength
         << ", constr_violation, rank, crowding_distance\n";
    fpt4 << "# of objectives = "    << nobj
         << ", # of constraints = " << ncon
         << ", # of real_var = "    << nreal
         << ", # of bits of bin_var = " << bitlength
         << ", constr_violation, rank, crowding_distance\n";
}

void NSGA2::report_parameters(std::ostream& os) const {
    os << "Population size = " << popsize
       << "\nNumber of generations = " << ngen
       << "\nNumber of objective functions = " << nobj
       << "\nNumber of constraints = " << ncon
       << "\nNumber of real variables = " << nreal;

    if (nreal != 0) {
        for (int i = 0; i<nreal; ++i) {
            os << "\nLower limit of real variable " << (i+1)
               << " = " << limits_realvar[i].first;
            os << "\nUpper limit of real variable " << (i+1)
               << " = " << limits_realvar[i].second;
        }
        os << "\nProbability of crossover of real variable = " << pcross_real;
        os << "\nProbability of mutation of real variable = " << pmut_real;
        os << "\nDistribution index for crossover = " << eta_c;
        os << "\nDistribution index for mutation = " << eta_m;
    }
    
    os << "\nNumber of binary variables = " << nbin;
    if (nbin != 0) {
        for (int i = 0; i<nbin; ++i) {
            os << "\nNumber of bits for binary variable " << (i+1)
               << " = " << nbits[i];
            os << "\nLower limit of real variable " << (i+1)
               << " = " << limits_binvar[i].first;
            os << "\nUpper limit of real variable " << (i+1)
               << " = " << limits_binvar[i].second;
        }
        os << "Probability of crossover of binary variable = " << pcross_bin;
        os << "Probability of mutation of binary variable = " << pmut_bin;
    }
    os << "\nSeed for random number generator = " << seed << endl;
}

void NSGA2::report_pop(const population& pop, std::ostream& os) const {

    pop.report(os);
    
}

