/* Rank assignment routine */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to assign rank and crowding distance to a population of size pop_size*/
void assign_rank_and_crowding_distance (population *new_pop)
{
    int flag;
    int i;
    int end;
    int front_size;
    int rank=1;
    list *orig;
    list *cur;
    list *temp1, *temp2;
    orig = (list *)malloc(sizeof(list));
    cur = (list *)malloc(sizeof(list));
    front_size = 0;
    orig->index = -1;        // orig empty : -1, null
    orig->parent = NULL;
    orig->child = NULL;
    cur->index = -1;         // cur empty  : -1, null
    cur->parent = NULL;
    cur->child = NULL;
    temp1 = orig;
    for (i=0; i<popsize; i++) // orig will have: -1,0,1,2,...,popsize-1
    {
        insert (temp1,i);
        temp1 = temp1->child;
    }
    do // for all elements in orig (i.e. all population)
    {
        if (orig->child->child == NULL) // if list has only one item, (item 0)
        {
            new_pop->ind[orig->child->index].rank = rank;
            new_pop->ind[orig->child->index].crowd_dist = INF;
            break;
        }
        temp1 = orig->child;        // temp1 = first element in orig
        insert (cur, temp1->index); // put first element of orig in cur
        front_size = 1;
        temp2 = cur->child;         // temp2 = first element in cur
        temp1 = del (temp1);        // erase first element in orig
        temp1 = temp1->child;       // temp1 = next element in orig
        do // for every element left in orig
        { 
            temp2 = cur->child;
            do // for every element in cur or until one of them is < temp1
            {
                end = 0;
                flag = check_dominance (&(new_pop->ind[temp1->index]), &(new_pop->ind[temp2->index]));
                if (flag == 1) // temp1 < temp2 : temp1 dominates temp2
                {
                    insert (orig, temp2->index); // put non dominant in orig
                    temp2 = del (temp2);         // delete dominant from cur
                    front_size--;
                    temp2 = temp2->child;
                }
                if (flag == 0) // temp1 and temp2 nondominated
                {
                    temp2 = temp2->child; // continue to next element
                }
                if (flag == -1) // temp2 < temp1: temp2 dominates temp1
                { 
                    end = 1;              // quit iterating over cur
                }
            }
            while (end!=1 && temp2!=NULL);
            if (flag == 0 || flag == 1) // last comparison was (temp1 <= temp2)
            {
                insert (cur, temp1->index); // insert dominant 
                front_size++;
                temp1 = del (temp1);
            }
            temp1 = temp1->child;
        }
        while (temp1 != NULL);
        temp2 = cur->child;
        do
        {
            new_pop->ind[temp2->index].rank = rank;
            temp2 = temp2->child;
        }
        while (temp2 != NULL);
        assign_crowding_distance_list (new_pop, cur->child, front_size);
        temp2 = cur->child;
        do
        {
            temp2 = del (temp2);
            temp2 = temp2->child;
        }
        while (cur->child !=NULL);
        rank+=1;
    }
    while (orig->child!=NULL);
    free (orig);
    free (cur);
    return;
}
