#include "nsga2/global.h"
#include "rand.h"

#include <cmath> // only for test problem def
#include <iostream>

using namespace nsga2;
using namespace std;

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

void test_problem (double *xreal,
                   double *xbin,
                   int **gene,
                   double *obj,
                   double *constr) {
    // obj[0] = 4.0*(xreal[0]*xreal[0] + xreal[1]*xreal[1]);
    // obj[1] = std::pow((xreal[0]-5.0),2.0) + std::pow((xreal[1]-5.0),2.0);
    // constr[0] = 1.0 - (std::pow((xreal[0]-5.0),2.0) + xreal[1]*xreal[1])/25.0;
    // constr[1] = (std::pow((xreal[0]-8.0),2.0) + std::pow((xreal[1]+3.0),2.0))/7.7 - 1.0;

    obj[0] = pow(xreal[0],2.0);
    obj[1] = pow((xreal[0]-2.0),2.0);
    return;
}

void individual::evaluate() {

    // workaround to respect the signature of test_problem and its (int**)
    std::vector<int*> tmp(gene.size());
    for (unsigned i=0; i < gene.size(); ++i) {
        tmp[i] = &(gene[i][0]);
    }
    // TODO: change test_problem so one can change it to whatever
    test_problem (&xreal[0], &xbin[0], &tmp[0], &obj[0], &constr[0]);
    
    if (config->ncon) {
        for (int i = 0; i < config->ncon; ++i)
            if (constr[i] < 0.0)
                constr_violation += constr[i];
    } else {
        constr_violation = 0.0;
    }
}

// returns:  1 if this < b (this dominates b),
//          -1 if this > b (this is dominated by b),
//           0 if they are nondominated
int individual::check_dominance(const individual& b) const {
    
    if (constr_violation < 0 && b.constr_violation < 0) {
        // both have constraint violations
        
        if (constr_violation > b.constr_violation)
            return 1; // this violates less
        else if (constr_violation < b.constr_violation)
            return -1; // b violates less
        else
            return 0; // they both violate equally
        
    } else if (constr_violation < 0 && b.constr_violation == 0) {
        // this violates and b doesn't => b dominates
        
        return -1; 
        
    } else if (constr_violation == 0 && b.constr_violation < 0) {
        // this doesn't violate and b does => this dominates
        
        return 1; 
        
    } else {
        // no constraint violations
       
        int flag1 = 0, // to check if this has a smaller objective
            flag2 = 0; // to check if b    has a smaller objective
        
        for (int i=0; i<config->nobj; ++i) {
            if (obj[i] < b.obj[i]) {
                flag1 = 1;
            } else if (obj[i] > b.obj[i]) {
                flag2 = 1;
            }
        }
        
        if (flag1==1 && flag2==0) {
            // there is at least one smaller objective for this and none for b
            
            return 1;
            
        } else if (flag1==0 && flag2==1) {
            // there is at least one smaller objective for b and none for this
            
            return -1;
            
        } else {
            // no smaller objective or both have one smaller
            
            return 0;
        }
               
    }
}

// returns num_mut_real, num_mut_bin
std::pair<int,int> individual::mutate() {
    std::pair<int,int> num_mut = std::make_pair(0,0);
    if (config->nreal)
        num_mut.first  += real_mutate();
    if (config->nbin)
        num_mut.second += bin_mutate();
    return num_mut;
}

int individual::real_mutate() {
    int j;
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;
    int num_mut = 0;
    for (j=0; j<config->nreal; j++) {
        if (randomperc() <= config->pmut_real) {
            y = xreal[j];
            yl = config->limits_realvar[j].first;
            yu = config->limits_realvar[j].second;
            delta1 = (y-yl)/(yu-yl);
            delta2 = (yu-y)/(yu-yl);
            rnd = randomperc();
            mut_pow = 1.0/(config->eta_m+1.0);
            if (rnd <= 0.5) {
                xy = 1.0-delta1;
                val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(config->eta_m+1.0)));
                deltaq =  pow(val,mut_pow) - 1.0;
            } else {
                xy = 1.0-delta2;
                val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(config->eta_m+1.0)));
                deltaq = 1.0 - (pow(val,mut_pow));
            }
            y = y + deltaq*(yu-yl);
            if (y<yl)
                y = yl;
            if (y>yu)
                y = yu;
            xreal[j] = y;
            num_mut+=1;
        }
    }
    return num_mut;
}

int individual::bin_mutate() {
    int j, k;
    double prob;
    int num_mut = 0;
    for (j=0; j<config->nbin; j++) {
        for (k=0; k<config->nbits[j]; k++) {
            prob = randomperc();
            if (prob <=config->pmut_bin) {
                if (gene[j][k] == 0) {
                    gene[j][k] = 1;
                } else {
                    gene[j][k] = 0;
                }
                num_mut+=1;
            }
        }
    }
    return num_mut;
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
                       const int nobj,
                       const double pmut_real,
                       const double pmut_bin,
                       const double eta_m) throw (nsga2::nsga2exception) :
    ind_config() {

    ind_config.nreal          = nreal;
    ind_config.nbin           = nbin;
    ind_config.nobj           = nobj;
    ind_config.ncon           = ncon;
    ind_config.nbits          = nbits;
    ind_config.limits_realvar = limreal;
    ind_config.limits_binvar  = limbin;
    ind_config.pmut_real      = pmut_real;
    ind_config.pmut_bin       = pmut_bin;
    ind_config.eta_m          = eta_m;
    
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

void population::evaluate() {
    std::vector<individual>::iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {
        it->evaluate();
    }
}

void population::fast_nds() { // TODO: should return F?
    std::vector< std::vector<int> >  F(1);
    int i,j;
    for (i = 0; i < ind.size(); ++i) {
        
        individual& p = ind[i];
        p.dcounter  = 0;
        p.dominated.clear();

        for (j = 0; j < ind.size(); ++j) {

            individual& q = ind[j];
            
            int compare = p.check_dominance(q);
            if (compare == 1) { // p dominates q
                p.dominated.push_back(j);
            } else if (compare == -1) { // q dominates p
                p.dcounter += 1;
            }
        }

        if (p.dcounter == 0) {
            p.rank = 1;
            F[0].push_back(i);
        }
        
    }

    int fi = 1;
    while (F[fi-1].size() > 0) {

        std::vector<int>& Fi = F[fi-1];
        std::vector<int> Q;
        for (i = 0; i < Fi.size(); ++i) {
            
            individual& p = ind[Fi[i]];
            
            for (j = 0; j < p.dominated.size() ; ++j) {

                
                individual& q = ind[p.dominated[j]];
                q.dcounter -= 1;

                if (q.dcounter == 0) {
                    q.rank = fi+1;
                    Q.push_back(p.dominated[j]);
                }
            }
        }
        

        fi += 1;
        F.push_back(Q);
    }
    
}

void population::report(std::ostream& os) const {

    std::vector<individual>::const_iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {

        for (int j = 0; j < ind_config.nobj; ++j)
            os << it->obj[j] << '\t';
        for (int j = 0; j < ind_config.ncon; ++j)
            os << it->constr[j] << '\t';
        for (int j = 0; j < ind_config.nreal; ++j) 
            os << it->xreal[j] << '\t';
        for (int j = 0; j < ind_config.nbin; ++j)
            for (int k = 0; k < ind_config.nbits[j]; ++k) 
                os << it->gene[j][k] << '\t';

        os << it->constr_violation << '\t'
           << it->rank << '\t'
           << it->crowd_dist << '\n';
        
    }
}

std::pair<int,int> population::mutate() {
    std::pair<int,int>
        num_mut = std::make_pair(0,0),
        tmp = std::make_pair(0,0);
    std::vector<individual>::iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {
        tmp = it->mutate();
        num_mut.first  += tmp.first;
        num_mut.second += tmp.second;
    }
    return num_mut;
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
