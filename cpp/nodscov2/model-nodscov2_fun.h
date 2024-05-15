#ifndef MODEL_NODSCOV2_FUN__H
#define MODEL_NODSCOV2_FUN__H

#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;

// extern double beta;
// extern double epsilon;
// extern double nu;
// extern double mu;
// extern double tau;

RcppExport Rcpp::List List_encountered(
    Rcpp::String id,
    Rcpp::DataFrame interactions_t
);

RcppExport Rcpp::NumericVector Lambda_c (
    Rcpp::DataFrame lambda_tim1,
    Rcpp::DataFrame interaction_ti,
    Rcpp::DataFrame status_ti,
    const double beta,
    const double dt,
    const int ti
);

// Update ENV
// RcppExport DataFrame Update_environment(
//     DataFrame environment_tim1,
//     DataFrame interaction_tim1,
//     DataFrame status_tim1,
//     //DataFrame rooms,
//     double mu,
//     double nu,
//     double dt
//     //double t
// );



// List encountered

// FOI (lambda_c / lambda_e)

// Update STATUS

// UPDATE t --> GIVEN GLOBAL LISTS -> GIVE DATAFRAME FOR TIME t
#endif