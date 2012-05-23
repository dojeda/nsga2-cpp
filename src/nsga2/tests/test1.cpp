#include <nsga2/global.h>
#include <nsga2/NSGA2.h>

#include "rand.h"

#include <vector>
#include <iostream>
#include <utility>

using namespace std;


int main(int argc, char *argv[]) {
    std::vector<int> nbits(3,5);

    nsga2::individual_config conf;
    conf.nreal = 3;
    conf.nbin  = 3;
    conf.ncon  = 2;
    conf.nbits = nbits;
    conf.nobj  = 2;
    conf.limits_realvar.push_back(make_pair(0.1,1.2));
    conf.limits_realvar.push_back(make_pair(0.1,1.2));
    conf.limits_realvar.push_back(make_pair(0.1,1.2));
    conf.limits_binvar.push_back(make_pair(0.0,2.5));
    conf.limits_binvar.push_back(make_pair(0.0,2.5));
    conf.limits_binvar.push_back(make_pair(0.0,2.5));
    
    nsga2::individual ind(conf);

    randomize();
    ind.initialize();
    
    cout << "Individual: " << ind << endl;

    std::vector< std::pair<double,double> > limr = conf.limits_realvar;
    std::vector< std::pair<double,double> > limb = conf.limits_binvar;
    nsga2::population pop(4,3,3,2,nbits,limr,limb,2,0.7,0.1,10,NULL);

    pop.initialize();
    cout << "Population: " << pop << endl;

    nsga2::NSGA2 ga;

    cout << "GA: " << &ga << endl;
    
    return 0;
}
