/* Routines to display the population information using gnuplot */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <unistd.h>

# include "global.h"
# include "rand.h"

/* Function to display the current population for the subsequent generation */
void onthefly_display (population *pop, FILE *gp, int ii)
{
    int i;
    int flag;
    FILE *fpt;
    fpt = fopen("plot.out","w");
    flag = 0;
    for (i=0; i<popsize; i++)
    {
        if (pop->ind[i].constr_violation==0)
        {
            if (choice==4) {
                int j;
                for (j=0; j < nobj; ++j) {
                    fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
                }
                fprintf(fpt,"\n");
            } else if (choice!=3) {
                fprintf(fpt,"%e\t%e\n",pop->ind[i].obj[obj1-1],pop->ind[i].obj[obj2-1]);
            } else {
                fprintf(fpt,"%e\t%e\t%e\n",pop->ind[i].obj[obj1-1],pop->ind[i].obj[obj2-1],pop->ind[i].obj[obj3-1]);
            }
            fflush(fpt);
            flag = 1;
        }
    }
    if (flag==0)
    {
        printf("\n No feasible soln in this pop, hence no display");
    }
    else
    {
        if (choice==4) {
            fprintf(gp,
                    /* "set title 'Generation test #%d'\n" */
                    "set size 1.0, 1.0\n"
                    "set origin 0.0, 0.0\n"
                    "set multiplot\n"
                    /* "unset key\n " */
                    /* "plot 'plot.out' using 1:2 w points pointtype 1 pointsize 1\n" */
                    ,ii);
            int i,j;
            double step = 1.0/(nobj);
            for (i = 0; i < nobj; ++i) {
                for (j = 0; j < nobj; ++j) {
                    /* if (j==i) { */
                    /*     fprintf(gp, */
                    /*             "set size %.2f, %.2f\n" */
                    /*             "set origin %.2f, %.2f\n" */
                    /*             "unset key\n" */
                    /*             "set xrange [-1:1]\n" */
                    /*             "set yrange [-1:1]\n" */
                    /*             "set label \"hello at 0,0 center\"\n" */
                    /*             "plot 1/0\n" */
                    /*             , step, step */
                    /*             , 1.0-step*i, 1.0-step*j); */
                    /* } else { */
                        fprintf(gp,
                                "set title '%dG vs %dG (%.2f, %.2f)'\n"
                                "set size %.2f, %.2f\n"
                                "set origin %.2f, %.2f\n"
                                "unset key\n"
                                "plot 'plot.out' using %d:%d w p pt 1 ps 1\n"
                                , j, i, step*j, 1.0-step*(i+1)
                                , step, step
                                , step*j, 1.0-step*(i+1)
                                , j+1, i+1 );
                    /* } */
                }
            }
        }
        else if (choice!=3)
            fprintf(gp,"set title 'Generation #%d'\n unset key\n plot 'plot.out' w points pointtype 6 pointsize 1\n",ii);
        else
            fprintf(gp,"set title 'Generation #%d'\n set view %d,%d\n unset key\n splot 'plot.out' w points pointtype 6 pointsize 1\n",ii,angle1,angle2);
        fflush(gp);
    }
    fclose(fpt);
    return;
}
