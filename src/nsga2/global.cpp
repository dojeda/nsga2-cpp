#include "nsga2/global.h"
#include "rand.h"

using namespace nsga2;

individual::individual() throw () :
    rank(0),
    constr_violation(0),
    xreal(0),
    gene(0),
    xbin(0),
    obj(0),
    constr(0),
    crowd_dist(0),
    config(0) {}

individual::individual(const individual_config& c
                       /*const unsigned int nreal,
                       const unsigned int nbin,
                       const unsigned int ncon,
                       const std::vector<int>& nbits,
                       const unsigned int nobj*/) throw (nsga2::nsga2exception) :
    rank(0),
    constr_violation(0),
    xreal(0),
    gene(0),
    xbin(0),
    obj(0),
    constr(0),
    crowd_dist(0),
    config(&c) {
    
    xreal.resize(config->nreal,0);
    xbin.resize(config->nbin,0);
    gene.resize(config->nbin);
    if (config->nbits.size() != config->nbin)
        throw nsga2::nsga2exception("nbits size != nbin");
    for (int j = 0; j < config->nbin; ++j) {
        gene[j].resize(config->nbits[j],0);
    }
    obj.resize(config->nobj,0);
    constr.resize(config->ncon,0);
}


individual::~individual() {
}

void individual::initialize() throw (nsga2::nsga2exception) {
    if (!config)
        throw nsga2::nsga2exception("Individual not configured");
    
    for (int i = 0; i < config->nreal; ++i) {
        xreal[i] = rndreal(config->limits_realvar[i].first,
                           config->limits_realvar[i].second);
    }

    for (int i = 0; i < config->nbin; ++i) {
        for (int j = 0; j < config->nbits[i]; ++j) {
            gene[i][j] = randomperc() <= 0.5 ? 0 : 1;
        }
    }
}

void individual::decode() {
    int sum;
    for (int i = 0; i < config->nbin; ++i) {
        sum = 0;
        for (int j = 0; j < config->nbits[i]; ++j) {
            sum += (1 << (config->nbits[i]-1-j));  // TODO: check
        }

        xbin[i] = config->limits_binvar[i].first +
            (double)sum*( config->limits_binvar[i].second - config->limits_binvar[i].first) / (double)((1 << (config->nbits[i]))-1); // TODO: check
    }
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

population::population(const int size,
                       const int nreal,
                       const int nbin,
                       const int ncon,
                       const std::vector<int>& nbits,
                       const std::vector< std::pair<double,double> >& limreal,
                       const std::vector< std::pair<double,double> >& limbin,
                       const int nobj) throw (nsga2::nsga2exception) :
    ind_config() {

    ind_config.nreal = nreal;
    ind_config.nbin  = nbin;
    ind_config.ncon  = ncon;
    ind_config.nbits = nbits;
    ind_config.limits_realvar = limreal;
    ind_config.limits_binvar  = limbin;
    
    for (unsigned int i = 0; i < size; ++i) {
        ind.push_back(individual(ind_config));
    }
    
}

population::~population() {
}

void population::initialize() throw (nsga2::nsga2exception) {
    std::vector<individual>::iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {
        it->initialize();
    }
}

void population::decode() {
    std::vector<individual>::iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {
        it->decode();
    }
}

std::ostream& nsga2::operator<< (std::ostream& os, const population& pop) {
    os << "Population: {\n";
    std::vector<individual>::const_iterator it;
    for (it = pop.ind.begin(); it != pop.ind.end(); ++it) {
        os << *it;
    }
    os << '}';
    return os;
}
