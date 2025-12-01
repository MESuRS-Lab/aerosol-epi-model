#include "dev-model-nodscov2_fun.h"
#include <Rcpp.h>


#include <iostream>

// [[Rcpp::plugins(cpp11)]]

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::DataFrame simulation(
    Rcpp::List global_interaction,
    Rcpp::List global_environment,
    Rcpp::List global_data,
    Rcpp::DataFrame global_status,
    double beta_c,
    double beta_e,
    double B,
    double nu,
    double mu,
    String env_model,
    double dt,
    std::string intervention,
    double mu_int = 0,
    double nu_int = 0,
    double rel_trans_risk = 0
    
) {
    // INITIALIZATION
    double tau = 60 * 60 * 24; //seconds in a day
    double deltat = dt/(tau);

    Rcpp::DataFrame environment_ti = Get_t(global_environment, 0);
    Rcpp::DataFrame environment_tim1 = Get_t(global_environment, 0);
    
    Rcpp::CharacterVector ids_ti;    
    
    Rcpp::IntegerVector status_ti;
    Rcpp::IntegerVector status_tim1;

    Rcpp::DataFrame interaction_ti;
    Rcpp::DataFrame interaction_tim1;

    Rcpp::IntegerVector location_ti;
    Rcpp::IntegerVector location_tmi;

    Rcpp::IntegerVector interaction_with_patient_ti;
    Rcpp::IntegerVector interaction_with_patient_tmi;

    Rcpp::DataFrame global_data_t;
    Rcpp::IntegerVector info_ti;
    

    // INTERVENTIONS
    if (intervention == "None") {
        mu_int = mu;
        nu_int = nu;
        rel_trans_risk = 1.0;

    } else if (intervention == "Symptomatic masking") {
        mu_int = mu;
        if (nu == nu_int || rel_trans_risk == 1.0) stop("Symptomatic masking intervention needs a specific shedding and relative transmission risk rates.");

    } else if (intervention == "Universal masking") {
        mu_int = mu;
        if (nu == nu_int || rel_trans_risk == 1.0) stop("Universal masking intervention needs a specific shedding and relative transmission risk rates.");

    } else if (intervention == "Hand hygiene") {
        mu_int = mu;
        nu_int = nu;
        if (rel_trans_risk == 1.0) stop("Hand hygiene intervention needs a specific relative transmission risk rate.");

    } else if (intervention == "Improved ventilation patients" || intervention == "Improved ventilation hcws") {
        nu_int = nu;
        rel_trans_risk = 1.0;
        if (mu == mu_int) stop("Improved ventilation intervention needs a specific aerosol inactivation rate.");

    } else if (intervention == "Mixed1") {
        if (mu == mu_int || nu==nu_int || rel_trans_risk == 1.0) stop("Mixed intervention n°1 needs specific values for nu, mu, and relative transmission rates.");

    } else if (intervention == "Mixed2") {
        if (mu == mu_int || nu==nu_int || rel_trans_risk == 1.0) stop("Mixed intervention n°2 needs specific values for nu, mu, and relative transmission rates.");

    } else {
        stop("Intervention does not exist");
        
    }

    // Display simulation configuration
    Rcout << "-----------Simulation configuration-----------" << std::endl;
    Rcout << "Intervention: " << intervention << std::endl;
    Rcout << "Shedding rate: " << nu << std::endl;
    Rcout << "Modified shedding rate: " << nu_int << std::endl; 
    Rcout << "Inactivation rate: " << mu << std::endl; 
    Rcout << "Modified inactivation rate: " << mu_int << std::endl; 
    Rcout << "Relative transmission risk during close-contact: " << rel_trans_risk << std::endl;


    // Parameters for the discretized gamma distribution
    Rcpp::List params = param_gamma_discretized();

    ///////////////
    // R's t = 1 //
    ///////////////
    ids_ti = Get_t(global_data, 0)["id"];
    info_ti = Get_t(global_data, 0)["info"];
    location_ti = Get_t(global_data, 0)["location_ti"];
    status_tim1 = Get_status_t(global_status, ids_ti, 0); 
    interaction_with_patient_ti = Get_t(global_data, 0)["interaction_with_patient"];

    // Get id of index patient
    std::string index_patient;
    Rcpp::CharacterVector inf_by_t0 = global_status["inf_by"];
    Rcpp::CharacterVector ids_all = global_status["id"];
    for (size_t i=0; i<inf_by_t0.size(); i++) {
        if (inf_by_t0[i] == "INDEX") index_patient = ids_all[i];
    }
    Rcout << "Index patient: " << index_patient << std::endl;

    // Shedding of the index patient
    environment_ti["env"] = Update_environment(ids_ti, info_ti, environment_tim1, location_ti, status_tim1, interaction_with_patient_ti, 
        mu, mu_int, nu, nu_int, deltat, intervention, index_patient);
    global_environment[0] = environment_ti;
    // update status for time t = 1?
    
    //////////////////////////
    // SIMULATION (R's t>1) //
    //////////////////////////
    int sim_size = global_interaction.size();

    for (int t = 1; t < sim_size; t++){
        global_data_t = Get_t(global_data, t);
        ids_ti = global_data_t["id"];
        info_ti = global_data_t["info"];
        location_ti = global_data_t["location_ti"];
        interaction_with_patient_ti = global_data_t["interaction_with_patient"];

        interaction_ti = Get_t(global_interaction, t);
        status_tim1 = Get_status_t(global_status, ids_ti, t-1);
        status_ti = Get_status_t(global_status, ids_ti, t);

        ////////////////////////////
        // Update the environment //
        ////////////////////////////
        environment_tim1 = Get_t(global_environment, t-1);
        environment_ti = clone(environment_tim1);
        environment_ti["env"] = Update_environment(ids_ti, info_ti, environment_tim1, location_ti, status_tim1, interaction_with_patient_ti, 
            mu, mu_int, nu, nu_int, deltat, intervention, index_patient);
        global_environment[t] = environment_ti;

        ///////////////////
        // Update Lambda //
        ///////////////////
        Rcpp::DataFrame temp_global_data = clone(global_data_t);
        temp_global_data["lambda_e"] = Lambda_e(info_ti, location_ti, environment_ti, beta_e, B, env_model, deltat);
        temp_global_data["lambda_c"] = Lambda_c(ids_ti, interaction_ti, status_ti, beta_c, deltat, intervention, rel_trans_risk, index_patient);
        global_data[t] = temp_global_data;
        ///////////////////////
        // Update the status //
        ///////////////////////
        Rcpp::DataFrame temp = Update_status_bis(
          global_status, 
          status_tim1, 
          temp_global_data["lambda_c"], 
          temp_global_data["lambda_e"], 
          info_ti, 
          ids_ti, 
          interaction_ti, 
          location_ti, 
          params,
          t);
        global_status = clone(temp);
    }

    Rcpp::DataFrame res = global_status;
    return res;    
}
