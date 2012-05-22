#include <nsga2/global.h>
#include <nsga2/NSGA2.h>

#include "rand.h"

#include <vector>
#include <numeric>
#include <iostream>
#include <cstdlib>

using namespace std;


int main(int argc, char *argv[]) {

    if (argc<2) {
        cout << "Usage " << argv[0] << " random_seed" << endl;
        return 1;
    }
    
    seed = atof(argv[1]);
    if (seed<=0.0 || seed>=1.0) {
        cout << "Entered seed value is wrong, seed value must be in (0,1)" << endl;
        return 1;
    }

    // fpt1 = fopen("initial_pop.out","w");
    // fpt2 = fopen("final_pop.out","w");
    // fpt3 = fopen("best_pop.out","w");
    // fpt4 = fopen("all_pop.out","w");
    // fpt5 = fopen("params.out","w");
    // fprintf(fpt1,"# This file contains the data of initial population\n");
    // fprintf(fpt2,"# This file contains the data of final population\n");
    // fprintf(fpt3,"# This file contains the data of final feasible population (if found)\n");
    // fprintf(fpt4,"# This file contains the data of all generations\n");
    // fprintf(fpt5,"# This file contains information about inputs as read by the program\n");

    int popsize;
    cout << "Enter the problem relevant and algorithm relevant parameters ... " << endl;
    cout << "Enter the population size (a multiple of 4) : " << endl;
    cin >> popsize;
    if (popsize<4 || (popsize%4)!= 0) {
        cout << "population size read is : " << popsize << endl;
        cout << "Wrong population size entered, hence exiting" << endl;
        return 1;
    }

    int ngen;
    cout << "Enter the number of generations : " << endl;
    cin >> ngen;
    if (ngen<1) {
        cout << "number of generations read is : " << ngen << endl;
        cout << "Wrong nuber of generations entered, hence exiting" << endl;
        return 1;
    }

    int nobj;
    cout << "Enter the number of objectives : " << endl;
    cin >> nobj;
    if (nobj<1) {
        cout << "number of objectives entered is : " << nobj << endl;
        cout << "Wrong number of objectives entered, hence exiting" << endl;
        return 1;
    }

    int ncon;
    cout << "Enter the number of constraints : " << endl;
    cin >> ncon;
    if (ncon<0) {
        cout << "number of constraints entered is : " << ncon << endl;
        cout << "Wrong number of constraints enetered, hence exiting" << endl;
        return 1;
    }

    int nreal;
    cout << "Enter the number of real variables : " << endl;
    cin >> nreal;
    if (nreal<0) {
        cout << "number of real variables entered is : " << nreal << endl;
        cout << "Wrong number of variables entered, hence exiting" << endl;
        return 1;
    }

    std::vector< std::pair<double,double> > limits_realvar;
    double pcross_real = 0;
    double pmut_real   = 0;
    double eta_c       = 0;
    double eta_m       = 0;
    if (nreal != 0) {
        for (int i=0; i<nreal; i++) {
            double realmax, realmin;
            cout << "Enter the lower limit of real variable " << (i+1) << ':' << endl;
            cin >> realmin;
            cout << "Enter the upper limit of real variable " << (i+1) << ':' << endl;
            cin >> realmax;
            if (realmax <= realmin) {
                cout << "Wrong limits entered for the min and max bounds of real variable, hence exiting" << endl;
                return 1;
            }
            limits_realvar.push_back(std::make_pair(realmin,realmax));
        }
        cout << "Enter the probability of crossover of real variable (0.6-1.0) : " << endl;
        cin >> pcross_real;
        if (pcross_real<0.0 || pcross_real>1.0) {
            cout << "Probability of crossover entered is :" << pcross_real << endl;
            cout << "Entered value of probability of crossover of real variables is out of bounds, hence exiting" << endl;
            return 1;
        }
        
        cout << "Enter the probablity of mutation of real variables (1/nreal) : " << endl;
        cin >> pmut_real;
        if (pmut_real<0.0 || pmut_real>1.0) {
            cout << "Probability of mutation entered is : " << pmut_real << endl;
            cout << "Entered value of probability of mutation of real variables is out of bounds, hence exiting" << endl;
            return 1;
        }
        
        cout << "Enter the value of distribution index for crossover (5-20): " << endl;
        cin >> eta_c;
        if (eta_c<=0) {
            cout << "The value entered is : " << eta_c << endl;
            cout << "Wrong value of distribution index for crossover entered, hence exiting" << endl;
            return 1;
        }
        
        cout << "Enter the value of distribution index for mutation (5-50): " << endl;
        cin >> eta_m;
        if (eta_m<=0) {
            cout << "The value entered is : " << eta_m << endl;
            cout << "Wrong value of distribution index for mutation entered, hence exiting" << endl;
            return 1;
        }
    }

    int nbin;
    cout << "Enter the number of binary variables : " << endl;
    cin >> nbin;
    if (nbin<0) {
        cout << "number of binary variables entered is : " << nbin << endl;
        cout << "Wrong number of binary variables entered, hence exiting" << endl;
        return 1;
    }

    std::vector< int > nbits;
    std::vector< std::pair<double,double> > limits_binvar;
    double pcross_bin = 0;
    double pmut_bin   = 0;
    if (nbin != 0) {
        for (int i=0; i<nbin; i++) {
            double binmax, binmin;
            int numbits;
            
            cout << "Enter the number of bits for binary variable " << (i+1) << ':' << endl;
            cin >> numbits;
            if (numbits < 1) {
                cout << "Wrong number of bits for binary variable entered, hence exiting" << endl;
                return 1;
            }
            cout << "Enter the lower limit of binary variable " << (i+1) << ':' << endl;
            cin >> binmin;
            cout << "Enter the upper limit of binary variable %d : " << (i+1) << ':' << endl;
            cin >> binmax;
            if (binmax <= binmin) {
                cout << "Wrong limits entered for the min and max bounds of binary variable entered, hence exiting" << endl;
                return 1;
            }
            nbits.push_back(numbits);
            limits_binvar.push_back(std::make_pair(binmin,binmax));
        }
        
        cout << "Enter the probability of crossover of binary variable (0.6-1.0): " << endl;
        cin >> pcross_bin;
        if (pcross_bin<0.0 || pcross_bin>1.0) {
            cout << "Probability of crossover entered is : " << pcross_bin << endl;
            cout << "Entered value of probability of crossover of binary variables is out of bounds, hence exiting" << endl;
            return 1;
        }
        
        cout << "Enter the probability of mutation of binary variables (1/nbits): " << endl;
        cin >> pmut_bin;
        if (pmut_bin<0.0 || pmut_bin>1.0) {
            cout << "Probability of mutation entered is : " << pmut_bin << endl;
            cout << "Entered value of probability  of mutation of binary variables is out of bounds, hence exiting" << endl;
            return 1;
        }
    }

    if (nreal==0 && nbin==0) {
        cout << "Number of real as well as binary variables, both are zero, hence exiting" << endl;
        return 1;
    }

    //     choice=0;
    // printf("\n Do you want to use gnuplot to display the results realtime (0 for NO) (1 for yes) : ");
    // scanf("%d",&choice);
    // if (choice!=0 && choice!=1)
    // {
    //     printf("\n Entered the wrong choice, hence exiting, choice entered was %d\n",choice);
    //     exit(1);
    // }
    // if (choice==1)
    // {
    //     gp = popen(GNUPLOT_COMMAND,"w");
    //     if (gp==NULL)
    //     {
    //         printf("\n Could not open a pipe to gnuplot, check the definition of GNUPLOT_COMMAND in file global.h\n");
    //         printf("\n Edit the string to suit your system configuration and rerun the program\n");
    //         exit(1);
    //     }
    //     if (nobj==2)
    //     {
    //         printf("\n Enter the objective for X axis display : ");
    //         scanf("%d",&obj1);
    //         if (obj1<1 || obj1>nobj)
    //         {
    //             printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
    //             exit(1);
    //         }
    //         printf("\n Enter the objective for Y axis display : ");
    //         scanf("%d",&obj2);
    //         if (obj2<1 || obj2>nobj)
    //         {
    //             printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
    //             exit(1);
    //         }
    //         obj3 = -1;
    //     }
    //     else
    //     {
    //         printf("\n #obj > 2, 2D display or a 3D display ?, enter 2 for 2D and 3 for 3D :");
    //         scanf("%d",&choice);
    //         if (choice!=2 && choice!=3)
    //         {
    //             printf("\n Entered the wrong choice, hence exiting, choice entered was %d\n",choice);
    //             exit(1);
    //         }
    //         if (choice==2)
    //         {
    //             printf("\n Enter the objective for X axis display : ");
    //             scanf("%d",&obj1);
    //             if (obj1<1 || obj1>nobj)
    //             {
    //                 printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
    //                 exit(1);
    //             }
    //             printf("\n Enter the objective for Y axis display : ");
    //             scanf("%d",&obj2);
    //             if (obj2<1 || obj2>nobj)
    //             {
    //                 printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
    //                 exit(1);
    //             }
    //             obj3 = -1;
    //         }
    //         else
    //         {
    //             printf("\n Enter the objective for X axis display : ");
    //             scanf("%d",&obj1);
    //             if (obj1<1 || obj1>nobj)
    //             {
    //                 printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
    //                 exit(1);
    //             }
    //             printf("\n Enter the objective for Y axis display : ");
    //             scanf("%d",&obj2);
    //             if (obj2<1 || obj2>nobj)
    //             {
    //                 printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
    //                 exit(1);
    //             }
    //             printf("\n Enter the objective for Z axis display : ");
    //             scanf("%d",&obj3);
    //             if (obj3<1 || obj3>nobj)
    //             {
    //                 printf("\n Wrong value of Z objective entered, value entered was %d\n",obj3);
    //                 exit(1);
    //             }
    //             printf("\n You have chosen 3D display, hence location of eye required \n");
    //             printf("\n Enter the first angle (an integer in the range 0-180) (if not known, enter 60) :");
    //             scanf("%d",&angle1);
    //             if (angle1<0 || angle1>180)
    //             {
    //                 printf("\n Wrong value for first angle entered, hence exiting \n");
    //                 exit(1);
    //             }
    //             printf("\n Enter the second angle (an integer in the range 0-360) (if not known, enter 30) :");
    //             scanf("%d",&angle2);
    //             if (angle2<0 || angle2>360)
    //             {
    //                 printf("\n Wrong value for second angle entered, hence exiting \n");
    //                 exit(1);
    //             }
    //         }
    //     }
    // }

    cout << "Input data successfully entered, now performing initialization" << endl;
    // fprintf(fpt5,"\n Population size = %d",popsize);
    // fprintf(fpt5,"\n Number of generations = %d",ngen);
    // fprintf(fpt5,"\n Number of objective functions = %d",nobj);
    // fprintf(fpt5,"\n Number of constraints = %d",ncon);
    // fprintf(fpt5,"\n Number of real variables = %d",nreal);
    // if (nreal!=0) {
    //     for (i=0; i<nreal; i++) {
    //         fprintf(fpt5,"\n Lower limit of real variable %d = %e",i+1,min_realvar[i]);
    //         fprintf(fpt5,"\n Upper limit of real variable %d = %e",i+1,max_realvar[i]);
    //     }
    //     fprintf(fpt5,"\n Probability of crossover of real variable = %e",pcross_real);
    //     fprintf(fpt5,"\n Probability of mutation of real variable = %e",pmut_real);
    //     fprintf(fpt5,"\n Distribution index for crossover = %e",eta_c);
    //     fprintf(fpt5,"\n Distribution index for mutation = %e",eta_m);
    // }
    // fprintf(fpt5,"\n Number of binary variables = %d",nbin);
    // if (nbin!=0) {
    //     for (i=0; i<nbin; i++) {
    //         fprintf(fpt5,"\n Number of bits for binary variable %d = %d",i+1,nbits[i]);
    //         fprintf(fpt5,"\n Lower limit of binary variable %d = %e",i+1,min_binvar[i]);
    //         fprintf(fpt5,"\n Upper limit of binary variable %d = %e",i+1,max_binvar[i]);
    //     }
    //     fprintf(fpt5,"\n Probability of crossover of binary variable = %e",pcross_bin);
    //     fprintf(fpt5,"\n Probability of mutation of binary variable = %e",pmut_bin);
    // }
    // fprintf(fpt5,"\n Seed for random number generator = %e",seed);
    
    // int bitlength = std::accumulate(nbits.begin(), nbits.end(), 0);
    // fprintf(fpt1,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
    // fprintf(fpt2,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
    // fprintf(fpt3,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
    // fprintf(fpt4,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);

    nsga2::NSGA2 nsga2;
    
    nsga2.set_nreal(nreal);
    nsga2.set_nbin(nbin);
    nsga2.set_nobj(nobj);
    nsga2.set_ncon(ncon);
    nsga2.set_popsize(popsize);
    nsga2.set_ngen(ngen);
    nsga2.set_pcross_real(pcross_real);
    nsga2.set_pcross_bin(pcross_bin);
    nsga2.set_pmut_real(pmut_real);
    nsga2.set_pmut_bin(pmut_bin);
    nsga2.set_eta_c(eta_c);
    nsga2.set_eta_m(eta_m);
    nsga2.set_nbits(nbits);
    nsga2.set_limits_realvar(limits_realvar);
    nsga2.set_limits_binvar(limits_binvar);

    nsga2.initialize();
    
    return 0;
}
