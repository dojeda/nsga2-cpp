#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <vector>
#include <ostream>
#include <string>
#include <utility>

#include "nsga2/exception.h"

namespace nsga2 {
    
    struct individual {

        individual() throw ();
        
        individual(const unsigned int nreal,
                   const unsigned int nbin,
                   const unsigned int ncon,
                   const std::vector<int>& nbits,
                   const unsigned int nobj) throw (nsga2::nsga2exception);
        virtual ~individual();

        
        int rank;
        double constr_violation;
        std::vector<double> xreal;
        std::vector< std::vector<int> > gene;
        std::vector<double> xbin;
        std::vector<double> obj;
        std::vector<double> constr;
        double crowd_dist;

    private:
        friend std::ostream& operator<< (std::ostream& os, const individual& ind);

    };

    std::ostream& operator<< (std::ostream& os, const individual& ind);

    struct population {

        typedef population* ptr;
        
        population(const unsigned int size,
                   const unsigned int nreal,
                   const unsigned int nbin,
                   const unsigned int ncon,
                   const std::vector<int>& nbits,
                   const unsigned int nobj) throw (nsga2::nsga2exception);
        virtual ~population();

        std::vector<individual> ind;

    private:
        friend std::ostream& operator<< (std::ostream& os, const population& pop);
    };

    std::ostream& operator<< (std::ostream& os, const population& pop);
}

#endif /* GLOBAL_H_ */
