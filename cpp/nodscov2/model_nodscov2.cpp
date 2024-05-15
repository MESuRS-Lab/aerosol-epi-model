#include "model-nodscov2_fun.h"
#include <Rcpp.h>


// [[Rcpp::plugins(cpp11)]]

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::DataFrame simulation(
    Rcpp::List global_interaction,
    Rcpp::List global_localization,
    Rcpp::List global_environment,
    Rcpp::List global_status,
    Rcpp::List global_lambda,
    double beta,
    double epsilon,
    double nu,
    double mu,
    double dt,
    Rcpp::DataFrame info_patient_HCW
) {
    return(info_patient_HCW);
}

    // Rcpp::List clusters,
    // Rcpp::DataFrame inferred_admission,
    // Rcpp::DataFrame HCW_interacting_id,
    // Rcpp::List double_rooms,
    // Rcpp::Datetime begin_date,
    // Rcpp::Datetime end_date,
    // int n_subdivisions,