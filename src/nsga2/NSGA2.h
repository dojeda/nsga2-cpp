#ifndef NSGA2_H_
#define NSGA2_H_

#include <nsga2/global.h>

#include <vector>
#include <utility>

namespace nsga2 {

    class NSGA2 {
    public:
        NSGA2();
        virtual ~NSGA2();
        
    private:
        // Parameters to be defined by the user
        int nreal;
        int nbin;
        int nobj;
        int ncon;
        int popsize;
        double pcross_real;
        double pcross_bin;
        double pmut_real;
        double pmut_bin;
        double eta_c;
        double eta_m;
        int ngen;
        int nbinmut;
        int nrealmut;
        int nbincross;
        int nrealcross;
        std::vector<int> nbits;
        std::vector< std::pair<double,double> > limits_realvar;
        // std::vector<double> min_realvar;
        // std::vector<double> max_realvar;
        // double *min_binvar;
        // double *max_binvar;
        std::vector< std::pair<double,double> > limits_binvar;
        // int choice; // to be added later, maybe.
        // int obj1;
        // int obj2;
        // int obj3;
        // int angle1;
        // int angle2;

    private:
        int bitlength;
        // random generator?
        // FILE *fpt1;
        // FILE *fpt2;
        // FILE *fpt3;
        // FILE *fpt4;
        // FILE *fpt5;
        // FILE *gp;
        // population *parent_pop;
        // population *child_pop;
        // population *mixed_pop;
        population* parent_pop;
        population* child_pop;
        population* mixed_pop;
    };
    
}

#endif /* NSGA2_H_ */
