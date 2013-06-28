#include <nsga2/global.h>
#include <nsga2/NSGA2.h>

#include "rand.h"

#include <vector>
#include <iostream>
#include <utility>
#include <cmath>

#define PI 3.14159265358979323846264338327950288419716939937510582097494

using namespace std;

void evaluate_model(double *xreal,
                    double *xbin,
                    int **gene,
                    double *obj,
                    double *constr) {
    double x1 = xreal[0];//(xreal[0]-0.5)*2*5.12;
    double x2 = xreal[1];//(xreal[1]-0.5)*2*5.12;
    x1 = floor(x1*10)/10.0;
    x2 = floor(x2*10)/10.0;
    double A = 10;
    obj[0] = -(A*2 + ( (x1*x1 - A * std::cos(2*PI*x1)) +
		      (x2*x2 - A * std::cos(2*PI*x2))));
    return;
}


int main(int argc, char *argv[]) {

    nsga2::individual_config conf;
    // conf.nreal = 2;
    // conf.nbin  = 0;
    // conf.ncon  = 0;
    // //conf.nbits = 0;
    // conf.nobj  = 1;
    // conf.pmut_real = 1;
    // conf.eta_m = 10;
    conf.limits_realvar.push_back(make_pair(-5.12,5.12));
    conf.limits_realvar.push_back(make_pair(-5.12,5.12));
    
    int seed = time(0);
    cout << "Using seed " << seed << endl;

    nsga2::NSGA2 nsga2;
    nsga2.set_seed(seed);
    nsga2.set_nreal(2);
    nsga2.set_nbin(0);
    nsga2.set_nobj(1);
    nsga2.set_ncon(0); // add a constraint due to possible simulation failures
    nsga2.set_popsize(1000);
    nsga2.set_ngen(100);
    nsga2.set_pcross_real(1.0);
    nsga2.set_pcross_bin(0.0);
    nsga2.set_pmut_real(0.3);
    nsga2.set_pmut_bin(0.0);
    nsga2.set_eta_c(10);
    nsga2.set_eta_m(10);
    nsga2.set_epsilon_c(1e-14);
    //nsga2.set_nbits(0);
    nsga2.set_limits_realvar(conf.limits_realvar);
    // nsga2.set_limits_binvar(limits_binvar);
    nsga2.set_function(&evaluate_model);
    nsga2.set_crowdobj(false); // crowd over the parameters, not the objective functions
    //nsga2.set_crowdobj(true); // crowd over objective function
    // nsga2.set_seed(seed);
    // nsga2.set_popfunction(&evaluate_population);
    // nsga2.set_custom_report_function(&update_generation);
    // nsga2.set_nreport(10);
    nsga2.set_backup_filename(""); // no backup

    
    nsga2.initialize();
    nsga2.evolve();
    
    return 0;
}
