#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <vector>
#include <ostream>
#include <string>

#include "nsga2/exception.h"

namespace nsga2 {
    
    struct individual {

        individual(const unsigned int nreal,
                   const unsigned int nbin,
                   const unsigned int ncon,
                   const std::vector<int>& nbits,
                   const unsigned int nobj);
        virtual ~individual();

        
        int rank;
        double constr_violation;
        std::vector<double> xreal;
        std::vector< std::vector<int> > gene; // TODO: why double pointer? Investigate
        std::vector<double> xbin;
        std::vector<double> obj;
        std::vector<double> constr;
        double crowd_dist;

    private:
        friend std::ostream& operator<< (std::ostream& os, const individual& ind);

    };

    std::ostream& operator<< (std::ostream& os, const individual& ind);
}

#endif /* GLOBAL_H_ */
