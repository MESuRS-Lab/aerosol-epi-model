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
    Rcpp::List list_ind_r;
    int nb_inf_r; 

    for (int j = 0; j < lambda_tim1.nrows(); ++j){
        id = ids[j];
        nb_inf_r = 0;
        list_ind_r = List_encountered(id, interaction_ti); // 
        for (int i = 0; i < list_ind_r.size(); ++i){
            Rcpp::String id_r = list_ind_r[i];
            // Search for the index of individual encountered in ids vector 
            int index_r = -1;
            for (int k = 0; k < ids.size(); ++k) {
                if (ids[k] == id_r) {
                    index_r = k;
                    break;
                }
            }
            // if individual r is infected & we found its index (for safety)
            if (index_r != -1 && status[index_r] == 1) {
                nb_inf_r += 1;
            }
        }
        
        lambda_c_ti[j] = beta * dt * nb_inf_r;
    }

    return lambda_c_ti;
};


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::NumericVector Lambda_e (
    Rcpp::DataFrame lambda_tim1,
    Rcpp::DataFrame localisation_ti,
    Rcpp::DataFrame environment_ti,
    Rcpp::DataFrame rooms,
    Rcpp::IntegerVector info_patient_HCW, // "0" IF PATIENT, "1" IF HCW
    const double epsilon,
    const double dt
) {
    Rcpp::CharacterVector ids_lambda = lambda_tim1["id"];
    Rcpp::CharacterVector ids_ind_rooms = rooms["id"]; // IDS OF INDIVIDUALS
    Rcpp::IntegerVector ids_rooms = rooms["id_room"]; // IDS OF THE ROOMS
    Rcpp::CharacterVector ids_localization = localisation_ti["id"];

    Rcpp::NumericVector temp_lambda_e = lambda_tim1["lambda_e"];
    Rcpp::NumericVector lambda_e_ti = clone(temp_lambda_e);
    Rcpp::NumericVector environment = environment_ti["env"];

    Rcpp::IntegerVector localizations = localisation_ti["localization"];
    Rcpp::String id;

  for (int j = 0; j < lambda_tim1.nrows(); ++j){
    id = ids_lambda[j];
    // TWO CASES (PATIENTS AND HCWS)
    // CASE 1. IF INDIVIDUAL j IS A PATIENT --> ENVIRONMENT ACCORDING TO ITS ROOM
    if (info_patient_HCW[j] == 0){
        // Search for the index of patient's room
        int index_room = -1;
            for (int k = 0; k < ids_ind_rooms.size(); ++k) {
                if (ids_ind_rooms[k] == id) {
                    index_room = k;
                    break;
                }
            }
        // ENVIRONMENT AND ROOM HAVE THE SAME INDEX
        lambda_e_ti[j] = epsilon * dt * environment[index_room];
    }
    
    // CASE 2. IF INDIVIDUAL j IS A HCW --> ENVIRONMENT ACCORDING TO ITS LOCALIZATION
    // WARNING, LOCALIZATION DF INDEX != LAMBDA DF INDEX ETC
    if (info_patient_HCW[j] == 1){
        // Search for the room where the HCW is located
        int index_localisation = -1;
            for (int k = 0; k < ids_localization.size(); ++k) {
                if (ids_localization[k] == id) {
                    index_localisation = k;
                    break;
                }
            }
        int room_j = localizations[index_localisation];
        // Search for the index of this room
        int index_room = -1;
            for (int k = 0; k < ids_rooms.size(); ++k) {
                if (ids_rooms[k] == room_j) {
                    index_room = k;
                    break;
                }
            }
        // ENVIRONMENT AND ROOM HAVE THE SAME INDEX
        lambda_e_ti[j] = epsilon * dt * environment[index_room];


    }
  }
  return lambda_e_ti;
};
