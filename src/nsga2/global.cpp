#include "nsga2/global.h"

using namespace nsga2;

individual::individual(const unsigned int nreal,
                       const unsigned int nbin,
                       const unsigned int ncon,
                       const std::vector<int>& nbits,
                       const unsigned int nobj) :
    rank(0),
    constr_violation(0) {
    
    xreal.resize(nreal,0);
    xbin.resize(nbin,0);
    gene.resize(nbin);
    if (nbits.size() != nbin)
        throw nsga2::nsga2exception("nbits size != nbin");
    for (unsigned int j = 0; j < nbin; ++j) {
        gene[j].resize(nbits[j],0);
    }
    obj.resize(nobj,0);
    constr.resize(ncon,0);
}


individual::~individual() {
}

std::ostream& nsga2::operator<< (std::ostream& os, const individual& ind) {
    
    os << "{Individual rank=" << ind.rank
       << "\nconstr_violation=" << ind.constr_violation;

    os << "\nxreal=[";
    std::vector<double>::const_iterator it;
    for (it = ind.xreal.begin(); it != ind.xreal.end(); ++it) {
        os << *it;
        if (it+1 != ind.xreal.end())
            os << ",";
    }
    
    os << "]\ngene=";
    std::vector< std::vector<int> >::const_iterator it1;
    for (it1 = ind.gene.begin(); it1 != ind.gene.end(); ++it1) {
        const std::vector<int>& tmp = *it1;
        std::vector<int>::const_iterator it2;
        if (it1 != ind.gene.begin())
            os << "     "; // tab space
        for (it2 = tmp.begin(); it2 != tmp.end(); ++it2) {
            os << *it2;
            if (it2+1 != tmp.end())
                os << ',';
        }
        //       gene=
        os << '\n';
    }

    os << "xbin=";
    for (it = ind.xbin.begin(); it != ind.xbin.end(); ++it) {
        os << *it;
        if (it+1 != ind.xbin.end())
            os << ",";
    }

    os << "\nobj=";
    for (it = ind.obj.begin(); it != ind.obj.end(); ++it) {
        os << *it;
        if (it+1 != ind.obj.end())
            os << ",";
    }

    os << "\nconstr=";
    for (it = ind.constr.begin(); it != ind.constr.end(); ++it) {
        os << *it;
        if (it+1 != ind.constr.end())
            os << ",";
    }

    os << "\ncrowd_dist=" << ind.crowd_dist;
 
    os << " }";
    
    return os;
}
