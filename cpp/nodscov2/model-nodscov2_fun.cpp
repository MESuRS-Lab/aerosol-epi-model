#include "model-nodscov2_fun.h"
#include <Rcpp.h>
#include <iostream>

using namespace std;

// [[Rcpp::plugins(cpp11)]]

// R UNIQUE(X) FUNCTION
Rcpp::Environment base("package:base");
Function do_unique = base["unique"];


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List List_encountered(
    String id,
    DataFrame interactions_t
) {
  Rcpp::List list_id;
  Rcpp::CharacterVector from = interactions_t["from"];
  Rcpp::CharacterVector to = interactions_t["to"];

  for (int j = 0; j < interactions_t.nrows(); ++j) {
    if (from[j] == id) {
        String push = to[j];
        list_id.push_back(push);
    }
    if (to[j] == id) {
        String push = from[j];
        list_id.push_back(push);
    }
  }
  // Need to have unique list of ids
  list_id = do_unique(list_id);
  return list_id;
};


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::NumericVector Lambda_c (
    Rcpp::DataFrame lambda_tim1,
    Rcpp::DataFrame interaction_ti,
    Rcpp::DataFrame status_ti,
    const double beta,
    const double dt,
    const int ti
) {
    Rcpp::CharacterVector ids = lambda_tim1["id"];
    Rcpp::NumericVector temp_lambda_c = lambda_tim1["lambda_c"];
    Rcpp::NumericVector lambda_c_ti = clone(temp_lambda_c);
    Rcpp::IntegerVector status = status_ti["status"];
    Rcpp::String id;
    Rcpp::List liste_ind_r;
    int nb_inf_r; 

    for (int j = 0; j < lambda_tim1.nrows(); ++j){
        id = ids[j];
        nb_inf_r = 0;
        liste_ind_r = List_encountered(id, interaction_ti); // 
        for (int i = 0; i < liste_ind_r.size(); ++i){
            Rcpp::String identifiant_r = liste_ind_r[i];
            //search for the index of individual encountered in ids vector 
            int index_r = -1;
            for (int k = 0; k < ids.size(); ++k) {
                if (ids[k] == identifiant_r) {
                    index_r = k;
                    break;
                }
            }
            if (index_r != -1 && status[index_r] == 1) {
                nb_inf_r += 1;
            }
        }
        
        lambda_c_ti[j] = beta * dt * nb_inf_r;
    }

    return lambda_c_ti;
};

