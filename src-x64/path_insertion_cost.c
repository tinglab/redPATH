#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "matrix_pos.h"


/*
  * Calculate insertion cost for the insertion algorithms
*/

  SEXP path_insertion_cost(SEXP R_matrix, SEXP R_order, SEXP R_k) {

    int n = INTEGER(GET_DIM(R_matrix))[0];
    int length = LENGTH(R_order);
    int *order = INTEGER(R_order);
    int k = INTEGER(R_k)[0]-1;
    SEXP R_cost;
    PROTECT(R_cost = NEW_NUMERIC(length+1));

    double link_add1, link_add2, link_remove, link0, linkn;
    double *matrix = REAL(R_matrix);
    double *cost = REAL(R_cost);


    //if (length == 1) {
    //  cost[0] = matrix[M_POS(n, order[0]-1, k)];

    //}else{
      for (int i = 0; i < (length-1); i++) {

        link_add1 = matrix[M_POS(n, order[i]-1, k)];
        link_add2 = matrix[M_POS(n, k, order[i+1]-1)];
        link_remove = matrix[M_POS(n, order[i]-1, order[i+1]-1)];

        //Rprintf("Links: %f %f %f\n", link_add1, link_add2, link_remove);

        // check for infinity
        //if(link_add1 == R_NegInf
        //   || link_add2 == R_NegInf
        //   || link_remove == R_PosInf)
        //  cost[i+1] = R_NegInf;

        //else if(link_add1 == R_PosInf
        //        || link_add2 == R_PosInf
        //        || link_remove == R_NegInf)
        //  cost[i+1] = R_PosInf;

        //else
          cost[i+1] = link_add1 + link_add2 - link_remove;
      }
      link0 = matrix[M_POS(n, k, order[0]-1)];

      linkn = matrix[M_POS(n, order[length-1]-1, k)];
      //Rprintf("Links: %f %f\n", link0, linkn);
      cost[0] = link0;
      cost[length] = linkn;

    //}

    UNPROTECT(1);
    return R_cost;
  }
