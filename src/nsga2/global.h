#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <vector>
#include <ostream>
#include <string>
#include <utility>

#include "nsga2/exception.h"

#define EPS 1e-14

namespace nsga2 {

class NSGA2;

struct individual_config {
    int nreal;
    int nbin;
    int nobj;
    int ncon;
    // int popsize;
    // int ngen;
    // double pcross_real;
    // double pcross_bin;
    double pmut_real;
    double pmut_bin;
    // double eta_c;
    double eta_m;
    std::vector<int> nbits;
    std::vector< std::pair<double,double> > limits_realvar;
    std::vector< std::pair<double,double> > limits_binvar;

};
    
struct individual {

    individual() throw (); // needed for std::vector<individual> allocator
        
    individual(const individual_config& c) throw (nsga2::nsga2exception);
    virtual ~individual();

    // individual& operator=(const individual& ind);
    
    void initialize() throw (nsga2::nsga2exception);

    void decode();
    void evaluate();
    std::pair<int,int> mutate();
    int real_mutate();
    int bin_mutate();

    int check_dominance(const individual& b) const;
    
    int rank;
    double constr_violation;
    std::vector<double> xreal;
    std::vector< std::vector<int> > gene;
    std::vector<double> xbin;
    std::vector<double> obj;
    std::vector<double> constr;
    double crowd_dist;

    int dcounter; // domination counter n_p
    std::vector<int> dominated;

private:
    const individual_config* config;
    friend std::ostream& operator<< (std::ostream& os, const individual& ind);
};

std::ostream& operator<< (std::ostream& os, const individual& ind);

struct population {
        
    population(const int size,
               const int nreal,
               const int nbin,
               const int ncon,
               const std::vector<int>& nbits,
               const std::vector< std::pair<double,double> >& limreal,
               const std::vector< std::pair<double,double> >& limbin,
               const int nobj,
               const double pmut_real,
               const double pmut_bin,
               const double eta_m) throw (nsga2::nsga2exception);
    virtual ~population();

    void initialize() throw (nsga2::nsga2exception);

    void decode();
    void evaluate();
    void fast_nds();
    std::pair<int,int> mutate();

    void report(std::ostream& os) const;

    int size() const {
        return ind.size();
    };
        
    std::vector<individual> ind;

private:
    individual_config ind_config;
    friend std::ostream& operator<< (std::ostream& os, const population& pop);
};

std::ostream& operator<< (std::ostream& os, const population& pop);
}

#endif /* GLOBAL_H_ */
