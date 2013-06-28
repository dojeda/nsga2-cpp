#include <nsga2/global.h>
#include <nsga2/NSGA2.h>

#include "rand.h"

#include <vector>
#include <iostream>
#include <utility>

using namespace std;


int main(int argc, char *argv[]) {

    nsga2::individual_config conf;
    conf.nreal = 2;
    conf.nbin  = 0;
    conf.ncon  = 0;
    //conf.nbits = 0;
    conf.nobj  = 1;
    conf.pmut_real = 1;
    conf.eta_m = 10;
    conf.limits_realvar.push_back(make_pair(-10,10));
    conf.limits_realvar.push_back(make_pair(-10,10));
    
    nsga2::NSGA2 nsga2;
    nsga2.set_seed(time(0));
    nsga2.set_nreal(2);
    nsga2.set_nbin(0);
    nsga2.set_nobj(1);
    nsga2.set_ncon(0); // add a constraint due to possible simulation failures
    nsga2.set_popsize(4);
    nsga2.set_ngen(10);
    nsga2.set_pcross_real(0.0);
    nsga2.set_pcross_bin(0.0);
    nsga2.set_pmut_real(1.0);
    nsga2.set_pmut_bin(0.0);
    nsga2.set_eta_c(2);
    nsga2.set_eta_m(10);
    //nsga2.set_nbits(nbits);
    nsga2.set_limits_realvar(conf.limits_realvar);
    // nsga2.set_limits_binvar(limits_binvar);
    // nsga2.set_function(&dummy);
    // nsga2.set_crowdobj(false); // crowd over the parameters, not the objective functions
    // //nsga2.set_crowdobj(true); // crowd over objective function
    // nsga2.set_seed(seed);
    // nsga2.set_popfunction(&evaluate_population);
    // nsga2.set_custom_report_function(&update_generation);
    // nsga2.set_nreport(10);
    nsga2.set_backup_filename(""); // no backup


    nsga2::individual p1(conf); p1.initialize(); p1.xreal[0] = 2; p1.xreal[1] = 2;
    nsga2::individual p2(conf); p2.initialize(); p2.xreal[0] = 8; p2.xreal[1] = 8;
    //cout << "c1x,c1y,c2x,c2y" << endl;
    cout << "c1x,c1y" << endl;
    for (int i = 0; i < 10000; ++i) {
    	// nsga2::individual ind1(conf);
    	// ind1.initialize();
    	// cout << ind1.xreal[0] << endl;

	// nsga2::individual c1(conf);
	// nsga2::individual c2(conf);

	nsga2::individual c1(conf); c1.initialize(); c1.xreal[0] = 0; c1.xreal[1] = 0;
	c1.real_mutate();

	cout << c1.xreal[0] << ',' << c1.xreal[1] << endl;
	

	// nsga2.crossover(p1,p2,c1,c2);
	// cout << c1.xreal[0] << ',' << c1.xreal[1] << ',' << c2.xreal[0] << ',' << c2.xreal[1] << endl;
    }
    // nsga2::individual ind1(conf);
    // nsga2::individual ind2(conf);

    // //randomize();
    // ind1.initialize();
    // ind2.initialize();
    
    // cout << "P1: " << ind1 << endl;
    // cout << "P2: " << ind2 << endl;
    
    return 0;
}
