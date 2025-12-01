#include "dev-model-nodscov2_fun.h"
#include <Rcpp.h>
#include <iostream>

using namespace std;

// [[Rcpp::plugins(cpp11)]]
const int z = 10;
const double mIncub = (1.63) * (24 * 60 * 2); // * z;
const double sdIncub = (0.5) * (24 * 60 * 2); // * z;

const double m_incub_g = (4.07) * (24 * 60 * 2); // * z;
const double sd_incub_g = (2.12) * (24 * 60 * 2); // * z;
const double shape_incub_g = pow(m_incub_g, 2) / pow(sd_incub_g, 2);
const double scale_incub_g = pow(sd_incub_g, 2) / m_incub_g;

const double mInf = (1.63) * (24 * 60 * 2); // * z;
const double sdInf = (0.5) * (24 * 60 * 2); // * z;

const double m_inf_g = (4.07) * (24 * 60 * 2); // * z;
const double sd_inf_g = (2.12) * (24 * 60 * 2); // * z;
const double shape_inf_g = pow(m_inf_g, 2) / pow(sd_inf_g, 2);
const double scale_inf_g = pow(sd_inf_g, 2) / m_inf_g;


// R UNIQUE(X) FUNCTION
Rcpp::Environment base("package:base");
Function do_unique = base["unique"];
Function do_sample = base["sample"];



// info 
    // 1 = HCW
    // 0 = Patient
// status 
    // 0 = Susceptible
    // 1 = Exposed non infectious
    // 2 = Exposed and infectious
    // 3 = Symptomatic and infectious
    // 4 = Recovered


// [[Rcpp::export]]
Rcpp::IntegerVector Get_status_t(
    const Rcpp::DataFrame& global_status,
    const Rcpp::StringVector& ids_ti,
    const int& t
) {
    Rcpp::IntegerVector t_inf = global_status["t_inf"];
    Rcpp::IntegerVector t_incub = global_status["t_incub"];
    Rcpp::IntegerVector t_infectious_start = global_status["t_infectious_start"];
    Rcpp::IntegerVector t_recover = global_status["t_recover"];
    Rcpp::CharacterVector ids_status = global_status["id"];
    
    int n_ind_ti = ids_ti.size();
    Rcpp::IntegerVector status_ti(n_ind_ti, -1);
    
    for(int j = 0; j < n_ind_ti; j++) {
        Rcpp::String id_j = ids_ti[j];
        int index_id = -1;
        // Find the index in global_status where the 'id' matches
        for (int k = 0; k < t_inf.size(); k++) {
            if (id_j == ids_status[k]) {
                index_id = k;
                break;
            }
        }
    
        if (index_id == -1) {
            // Identifier not found in global_status, handle as needed
            status_ti[j] = -1; // Example: Set status to -1 for not found
        } else {
            // Determine status based on t_inf, t_incub, t_recover
            if (t_inf[index_id] != -1 && (t+1) >= t_inf[index_id] && (t+1) <= t_incub[index_id] && t_infectious_start[index_id] > (t+1)) {
                status_ti[j] = 1; // individual j is EXPOSED and NON INFECTIOUS
            } else if (t_inf[index_id] != -1 && (t+1) >= t_inf[index_id] && (t+1) <= t_incub[index_id] && t_infectious_start[index_id] <= (t+1)) {
                status_ti[j] = 2; // individual j is EXPOSED and INFECTIOUS
            } else if (t_inf[index_id] != -1 && (t+1) > t_incub[index_id] && (t+1) <= t_recover[index_id]) {
                status_ti[j] = 3; // individual j is SYMPTOMATIC and INFECTIOUS
            } else if (t_inf[index_id] != -1 && (t+1) > t_recover[index_id]) {
                status_ti[j] = 4; // individual j is RECOVERED
            } else {
                status_ti[j] = 0; // individual j is SUSCEPTIBLE
            }
        }
    }

    return status_ti;
}

//////////////////////////////////////////////
// [[Rcpp::export]]
int Get_status_j(
    const Rcpp::String& id,
    const Rcpp::DataFrame& global_status,
    const int& t
) {
    int status_j = -1; // if returns -1 --> error (id not in global_status)
    Rcpp::CharacterVector ids = global_status["id"];
    Rcpp::IntegerVector t_inf = global_status["t_inf"];
    Rcpp::IntegerVector t_incub = global_status["t_incub"];
    Rcpp::IntegerVector t_infectious_start = global_status["t_infectious_start"];
    Rcpp::IntegerVector t_recover = global_status["t_recover"];
    int index_j = -1;
    for (int k = 0; k < global_status.nrows(); k++){
        if (id == ids[k]){
            index_j = k;
        }
    }
    // WE CHECK THE STATUS FOR ONLY INDIVIDUAL J (with index_j)
    if (index_j != -1){
        if (t_inf[index_j] != -1 && (t+1) >= t_inf[index_j] && (t+1) <= t_incub[index_j] && t_infectious_start[index_j] > (t+1)){
            status_j = 1; // individual j is EXPOSED and NON INFECTIOUS
        } else if (t_inf[index_j] != -1 && (t+1) >= t_inf[index_j] && (t+1) <= t_incub[index_j] && t_infectious_start[index_j] <= (t+1)) {
            status_j = 2; // individual j is EXPOSED and INFECTIOUS
        } else if (t_inf[index_j] != -1 && (t+1) > t_incub[index_j] && (t+1) <= t_recover[index_j]){ //cpp index begins at 0 & R's at 1, we chose to use R's index for time
            status_j = 3; // individual j is INFECTIOUS and SYMPTOMATIC
        } else if (t_inf[index_j] != -1 && (t+1) > t_recover[index_j]){
            status_j = 4; // individual j is RECOVERED
        } else{
            status_j = 0; // individual j is SUSCEPTIBLE
        }
    }
    
    return status_j;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::String Sample_inf(
    const Rcpp::String& id,
    const Rcpp::List& list_inf_encountered,
    const Rcpp::CharacterVector& ids_ti,
    const Rcpp::IntegerVector& location_ti,
    const double& lambda_e_j,
    const double& lambda_c_j
) {
    if (list_inf_encountered.size() == 0){ // NO CONTACTS -> ENV
        int room_j = Get_loc_j(id, location_ti, ids_ti);
        Rcpp::String res = "ENVIRONMENT-";
        res += std::to_string(room_j);
        return res;
        
    } if (lambda_e_j <= 1e-10 && list_inf_encountered.size() > 0){ // NO ENV --> CONTACT
        // Uniform weights are implicitly used in Rcpp::sample if not provided
        Rcpp::String sampled = Rcpp::sample(list_inf_encountered, 1, true)[0];
        Rcpp::String res = "CONTACT-";
        res += sampled;
        return res; 

    } else { // CONTACT AND ENV
        int n_ind = list_inf_encountered.size();
        double ind_weight; // lambda_c !=0 so n_ind!=0
        double env_weight;
        if (Rcpp::traits::is_infinite<REALSXP>(lambda_e_j) == 1) {
          ind_weight = 0;
          env_weight = 1;
        } else {
          ind_weight = lambda_c_j / (n_ind * (lambda_c_j + lambda_e_j)); // lambda_c !=0 so n_ind!=0
          env_weight = lambda_e_j / (lambda_c_j + lambda_e_j);
        }
        
        // WEIGHTS VECTOR
        Rcpp::NumericVector weights(n_ind + 1, ind_weight); // Initialize with ind_weight, CAUTION, R its rep(1,3) but cpp its v(3,1)
        weights[0] = env_weight;
        
        // ELEMENTS VECTOR
        Rcpp::CharacterVector elements(n_ind + 1);
        elements[0] = "ENVIRONMENT-";
        for (int i = 0; i < n_ind; ++i) {
            Rcpp::String id_patient_i = list_inf_encountered[i];
            elements[i + 1] = id_patient_i;
        }
        // SAMPLE
        Rcpp::String sampled = Rcpp::sample(elements, 1, true, weights)[0];
        
        // RETURN THE CAUSE OF INFECTION
        if (sampled == "ENVIRONMENT-") {
            int room_j = Get_loc_j(id, location_ti, ids_ti);
            Rcpp::String res = "ENVIRONMENT-";
            res += std::to_string(room_j);
            return res;
        } else {
            Rcpp::String res = "CONTACT-";
            res += sampled;
            return res;
        }
    }
};

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::DataFrame Update_status_bis(
    const Rcpp::DataFrame& global_status,
    const Rcpp::IntegerVector& status_tim1,
    const Rcpp::NumericVector& lambda_c_ti,
    const Rcpp::NumericVector& lambda_e_ti,
    const Rcpp::IntegerVector& info_ti,
    const Rcpp::CharacterVector& ids_ti,
    const Rcpp::DataFrame& interactions_ti,
    const Rcpp::IntegerVector& location_ti,
    const Rcpp::List& params,
    const int& t
) {
    Rcpp::CharacterVector ids_total = global_status["id"];
    Rcpp::IntegerVector t_inf_tim1 = global_status["t_inf"];
    Rcpp::IntegerVector t_incub_tim1 = global_status["t_incub"];
    Rcpp::IntegerVector t_infectious_start_tim1 = global_status["t_infectious_start"];
    Rcpp::LogicalVector presymp_trans_tim1 = global_status["presymp_trans"];
    Rcpp::IntegerVector t_recover_tim1 = global_status["t_recover"];
    Rcpp::CharacterVector inf_by_tim1 = global_status["inf_by"];
    Rcpp::IntegerVector inf_room_tim1 = global_status["inf_room"];
    
    Rcpp::DataFrame global_status_updated = clone(global_status);
    Rcpp::IntegerVector t_inf_ti = clone(t_inf_tim1);
    Rcpp::IntegerVector t_incub_ti = clone(t_incub_tim1);
    Rcpp::IntegerVector t_infectious_start_ti = clone(t_infectious_start_tim1);
    Rcpp::LogicalVector presymp_trans_ti = clone(presymp_trans_tim1);
    Rcpp::IntegerVector t_recover_ti = clone(t_recover_tim1);
    Rcpp::CharacterVector inf_by_ti = clone(inf_by_tim1);
    Rcpp::IntegerVector inf_room_ti = clone(inf_room_tim1);
    
    Rcpp::NumericVector FOI(lambda_c_ti.size(), 1);
    
    for (int k=0; k<lambda_c_ti.size(); k++) {
      if (Rcpp::traits::is_infinite<REALSXP>(lambda_e_ti[k]) == 0) {
        //Rcout << lambda_e_ti[k] << " " << FOI[k] << std::endl;
        FOI[k] += -exp(- (lambda_c_ti[k]  + lambda_e_ti[k]));
      }
    }

    for (int j=0; j < ids_ti.size(); j++){
        Rcpp::String id_j = ids_ti[j];
        int index_j = -1;
        for (int k = 0; k < global_status.nrows(); k++){
            if (id_j == ids_total[k]){
                index_j = k;
            }
        }

        if (status_tim1[j] == 0 && R::runif(0, 1) <= FOI[j]){
            t_inf_ti[index_j] = t+1; // C++ INDEX BEGINS AT 0 / R BEGINS AT 1
            t_incub_ti[index_j] = t_inf_ti[index_j] + Incub_period_gamma_discretized(params);
            
            // Presymptomatic acquisition 
            LogicalVector x={true, false};
            NumericVector probs={0.2, 0.8};
            presymp_trans_ti[index_j] = false + sample(x, 1, false, probs)[0];

            if (presymp_trans_ti[index_j]) {
                int min_start = t_incub_ti[index_j] - 3 * 24 * 60 * 2;
                t_infectious_start_ti[index_j] = runif_int(max({min_start, t_inf_ti[index_j]}), t_incub_ti[index_j]-1);
            } else {
                t_infectious_start_ti[index_j] = t_incub_ti[index_j];
            }

            // End of infectious period
            t_recover_ti[index_j] = t_infectious_start_ti[index_j] + runif_int(3, 7) * 24 * 60 * 2; // C++ INDEX BEGINS AT 0 / R BEGINS AT 1

            // ROOM WHERE j IS INFECTED //
            inf_room_ti[index_j] = location_ti[j];
            
            // CAUSE OF INFECTION
            Rcpp::List list_inf_encountered = List_inf_encountered(ids_ti[j], interactions_ti, global_status, t);
            double lambda_e_j = lambda_e_ti[j];
            double lambda_c_j = lambda_c_ti[j];
            Rcpp::String sampled_j = Sample_inf(ids_ti[j], list_inf_encountered, ids_ti, location_ti, lambda_e_j, lambda_c_j);
            inf_by_ti[index_j] = sampled_j;
        }
    }
    global_status_updated["t_inf"] = t_inf_ti;
    global_status_updated["t_incub"] = t_incub_ti;
    global_status_updated["presymp_trans"] = presymp_trans_ti;
    global_status_updated["t_infectious_start"] = t_infectious_start_ti;
    global_status_updated["t_recover"] = t_recover_ti;
    global_status_updated["inf_room"] = inf_room_ti;
    global_status_updated["inf_by"] = inf_by_ti;

    return global_status_updated;
};
            

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::DataFrame Get_t(
    const Rcpp::List& Global_list,
    const int& t
){
    // WARNING: c++ COUNTS STARTS AT 0
    Rcpp::DataFrame df(Global_list[t]);
    return df;
};


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::NumericVector Update_environment(
    const Rcpp::CharacterVector& ids_ti,
    const Rcpp::IntegerVector& info_ti,
    const Rcpp::DataFrame& environment_tim1,
    const Rcpp::IntegerVector& location_ti,
    const Rcpp::IntegerVector& status_tim1,
    const Rcpp::IntegerVector& interaction_with_patient_ti,
    const double& mu,
    const double& mu_int,
    const double& nu,
    const double& nu_int,
    const double& deltat,
    const std::string& intervention,
    const std::string& index_patient
) {
    // THERE IS MULTIPLE WAYS TO ACHIEVE THIS
    // A. (naive) for loop on room ( for loop on individuals (check if patient/HCW and if the location == room r then check if infected etc))
    // B. for loop on room (Patient assigned to this room is infected?  THEN for loop on location(is the location == room r ? and is the individual infected?))
    // C. FIRST loop on rooms/patients (is the patient infected? if so, update the env) THEN another for loop on location (update the env for each row if the individual is infected)
    // C. FIRST loop on rooms/patients (is the patient infected? if so, update the env) THEN another loop on infected patients etc.
    // THIS FUNCTION USES THE C. METHOD
    // WARNING: IF THERE IS DOUBLE ROOMS, DONT DO THE DOUBLE EXPONENTIAL INACTIVATION OF E(t-1)
    
    Rcpp::IntegerVector rooms_environment = environment_tim1["id_room"];

    double individual_nu = nu;

    // EXPONENTIAL INACTIVATION 
    Rcpp::NumericVector env_tim1 = environment_tim1["env"];
    Rcpp::CharacterVector env_loc = environment_tim1["room"];
    Rcpp::NumericVector env_ti(env_tim1.size());
    
    if (intervention == "Mixed1" || intervention == "Improved ventilation patients") {
      for (int j = 0; j < env_tim1.size(); j++) {
        if (env_loc[j] != "Medical Staff Room" && env_loc[j] != "Office" && env_loc[j] != "Paramedical Staff Room" && env_loc[j] != "Nursing station") {
          env_ti[j] = env_tim1[j] * exp(-mu_int * deltat);
        } else {
          env_ti[j] = env_tim1[j] * exp(-mu * deltat);
        }
      }
     
    } else if (intervention == "Mixed2" || intervention == "Improved ventilation hcws") {
        for (int j = 0; j < env_tim1.size(); j++) {
            if (env_loc[j] == "Medical Staff Room" || env_loc[j] == "Office" || env_loc[j] == "Paramedical Staff Room" || env_loc[j] == "Nursing station") {
                env_ti[j] = env_tim1[j] * exp(-mu_int * deltat);
            } else {
                env_ti[j] = env_tim1[j] * exp(-mu * deltat);
            }
        }

    } else {
      env_ti = as<NumericVector>(environment_tim1["env"]) * exp(-mu * deltat);
    }
    
    // INFECTING INDIVIDUALS SHEDDING
    for (int j = 0; j < ids_ti.size(); ++j) {

        if (status_tim1[j] == 2 || status_tim1[j] == 3){ // 2 = INFECTIOUS -> SHEDDING
          
          int room_j = location_ti[j];
            // Search for the index of the room in environment dataframe
            int index_room_j = -1;
            for (int k = 0; k < rooms_environment.size(); ++k) {
                if (rooms_environment[k] == room_j) {
                    index_room_j = k;
                    break;
                }
            }
            
            if (intervention == "Symptomatic masking" || intervention == "Mixed1" || intervention == "Mixed2") {
              if (info_ti[j] == 0 && status_tim1[j] == 3) { // INFECTIOUS PATIENT WITH SYMPTOMS
                individual_nu = nu_int;
              } 
              
              if (info_ti[j] == 1 && status_tim1[j] == 3 && (room_j < 50 || (room_j == 54 && interaction_with_patient_ti[j] == 1) ) ) { // INFECTIOUS HCW WITH SYMPTOMS in a patient room or in corridor in interaction with patients
                individual_nu = nu_int;
              }
            }

            if (intervention == "Universal masking" && !(ids_ti[j] == index_patient && status_tim1[j] == 2) ) {
                individual_nu = nu_int;
            }
            
            if(index_room_j != -1){
              env_ti[index_room_j] += individual_nu * deltat;
            }
          }
    }

  return env_ti;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List List_encountered(
    const Rcpp::String& id,
    const Rcpp::DataFrame& interaction_ti
) {
  Rcpp::List list_id;
  Rcpp::CharacterVector from = interaction_ti["from"];
  Rcpp::CharacterVector to = interaction_ti["to"];

  for (int j = 0; j < interaction_ti.nrows(); ++j) {
    if (from[j] == id) {
        Rcpp::String push = to[j];
        list_id.push_back(push);
    }
    if (to[j] == id) {
        Rcpp::String push = from[j];
        list_id.push_back(push);
    }
  }
  // Need to have unique list of ids
  list_id = do_unique(list_id);
  return list_id;
};


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List List_inf_encountered(
    const Rcpp::String& id,
    const Rcpp::DataFrame& interaction_ti,
    const Rcpp::DataFrame& global_status,
    const int& t
) {
    Rcpp::List list_encountered = List_encountered(id, interaction_ti);
    Rcpp::List list_inf_encountered;
    for (int j = 0; j < list_encountered.size(); j++){
        Rcpp::String id_encountered = list_encountered[j];
        // 0 = Susceptible, 1 = Exposed not infectious, 2 = Exposed infectious, 3 = Symptomatic infected, 4 = Recovered
        int status_j = Get_status_j(id_encountered, global_status, t);
        if (status_j == 2 || status_j == 3){ //INFECTIOUS
            list_inf_encountered.push_back(id_encountered);
        }
    }

    return list_inf_encountered;
};




//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::NumericVector Lambda_c(
    const Rcpp::CharacterVector& ids_ti, 
    const Rcpp::DataFrame& interaction_ti,
    const Rcpp::IntegerVector& status_ti,
    const double& beta_c,
    const double& deltat,
    const std::string& intervention,
    const double& rel_trans_risk,
    const std::string& index_patient
) {
    Rcpp::NumericVector lambda_c_ti (ids_ti.size(), 0);
    Rcpp::List list_ind_r;
    int nb_inf_r; 

    for (int j = 0; j < ids_ti.size(); ++j){
        std::string id_j = Rcpp::as<std::string>(ids_ti[j]);
        list_ind_r = List_encountered(id_j, interaction_ti);
        bool infectee_patient = id_j.find("PA") != std::string::npos;

        nb_inf_r = 0;
        
        if(list_ind_r.size() > 0){
            // Dont use List_inf_encountered because we dont need global_status (only check individuals encountered) for Lambda_c
            for (int i = 0; i < list_ind_r.size(); ++i){
                std::string id_r = list_ind_r[i];
                // Search for the index of individual encountered in ids vector 
                int index_r = -1;
                for (int k = 0; k < ids_ti.size(); ++k) {
                    if (ids_ti[k] == id_r) {
                        index_r = k;
                        break;
                    }
                }
                
                // if individual r is INFECTIOUS & we found its index (for safety)
                if (index_r != -1 && (status_ti[index_r] == 2 || status_ti[index_r] == 3)) { // 2 = INFECTIOUS not symptomatic
                                                                                             // 3 = INFECTIOUS symptomatic

                    bool infector_patient = id_r.find("PA") != std::string::npos;
                    bool hcw_to_patient = !infector_patient && infectee_patient;

                    if ( intervention == "Hand hygiene" && hcw_to_patient ) { // HCWs wash their hands when in interaction with patients (does not depend on symptom status)
                        nb_inf_r += rel_trans_risk; 

                    } else if ( (intervention == "Mixed1" || intervention == "Mixed2" || intervention == "Symptomatic masking") && hcw_to_patient) {  // HCWs wear a mask when in interaction with patients (does not depend on symptom status)
                        nb_inf_r += rel_trans_risk; 

                    } else if ( // Symptomatic patients wear a mask
                        infector_patient && status_ti[index_r] == 3 &&
                        (intervention == "Symptomatic masking" || intervention == "Mixed1" || intervention == "Mixed2") 
                        )  {

                        nb_inf_r += rel_trans_risk; 

                    } else if (intervention == "Universal masking" && !(id_r == index_patient && status_ti[index_r] == 2)) { // Universal masking except when the infector is the index patient before their symptoms onset 
                        nb_inf_r += rel_trans_risk; 

                    } else {
                        nb_inf_r += 1;
                    }
                }
            }
          
            lambda_c_ti[j] = beta_c * deltat * nb_inf_r;
        }
    }

    return lambda_c_ti;
};


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::NumericVector Lambda_e(
    const Rcpp::IntegerVector& info_ti,
    const Rcpp::IntegerVector& location_ti,
    const Rcpp::DataFrame& environment_ti,
    const double& beta_e,
    const double& B,
    const String& env_model,
    const double& deltat
) {
    Rcpp::NumericVector lambda_e_ti (location_ti.size(), 0);
    Rcpp::NumericVector environment = environment_ti["env"];
    Rcpp::IntegerVector rooms_environment = environment_ti["id_room"];
    Rcpp::IntegerVector rooms_volume = environment_ti["volume"];
    
    double individual_weight = 1;
    for (int j = 0; j < location_ti.size(); ++j){
        // // TWO CASES (PATIENTS AND HCWS)
        // if (info_ti[j] == 0){
        //     // CASE 1. IF INDIVIDUAL j IS A PATIENT --> 
        //     individual_weight = 1;
        // } else if (info_ti[j] == 1){
        //     // CASE 2. IF INDIVIDUAL j IS A HCW
        //     individual_weight = 1;
        // } else {
        //     individual_weight = 1;
        // }
        int room_j = location_ti[j];
            // Search for the index of patient's room
            int index_room = -1;
            for (int k = 0; k < rooms_environment.size(); ++k) {
                if (rooms_environment[k] == room_j) {
                    index_room = k;
                    break;
                }
            }
            // VIRAL LOAD threshold
            if (env_model == "linear") {
              lambda_e_ti[j] = beta_e * deltat * individual_weight * (B*deltat/rooms_volume[index_room]) * environment[index_room];
            }
            
            if (env_model == "exponential") {
              lambda_e_ti[j] = beta_e * deltat * individual_weight * (B*deltat/rooms_volume[index_room]) * expm1(environment[index_room]);
            }
            
            if (env_model == "log-linear") {
              lambda_e_ti[j] = beta_e * deltat * individual_weight * (B*deltat/rooms_volume[index_room]) * log10(environment[index_room]);
            }

        }  
        return lambda_e_ti;  
    };

//////////////////////////////////////////////
// [[Rcpp::export]]
int Incub_period_gamma() {
    // Incubation period -> Gamma distribution (shape,scale)
    double incubation_period_seconds = R::rgamma(shape_incub_g, scale_incub_g);
    int incubation_period_subdivisions = static_cast<int>(incubation_period_seconds / 30.0);
    
    return incubation_period_subdivisions;
};


//////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List param_gamma_discretized() {
    int max_time = 10 * 24 * 60 * 2; // Max incubation duration : 10 days
    IntegerVector to_sample_from;
    for (int k=0; k<max_time+1; k++) to_sample_from.push_back(k);
    
    NumericVector probs = Rcpp::pgamma(to_sample_from, shape_incub_g, scale_incub_g, true, false);
    probs.push_front(0.0);
    probs = diff(probs);
    double norm_cons = R::pgamma(max_time, shape_incub_g, scale_incub_g, true, false);
    for (int k=0; k < probs.size(); k++) probs[k] /= norm_cons;

    Rcpp::List params = List::create(Named("to_sample_from") = to_sample_from, _["probs"] = probs);

    return params;
}


//////////////////////////////////////////////
// [[Rcpp::export]]
int Incub_period_gamma_discretized(Rcpp::List params) {
    // Incubation period -> Gamma distribution (shape,scale)
    // double incubation_period_seconds = R::rgamma(shape_incub_g, scale_incub_g);
    //int incubation_period_subdivisions = static_cast<int>(incubation_period_seconds / 30.0);
    // return incubation_period_subdivisions;

    IntegerVector to_sample_from = params["to_sample_from"];
    NumericVector probs = params["probs"];
    int incub_period = sample(to_sample_from, 1, false, probs)[0];

    return incub_period;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
int Incub_period_lognormal() {
    // Incubation period -> Log-normal distribution (meanlog, sdlog)
    double incubation_period_seconds = R::rlnorm(mIncub, sd_incub_g);
    int incubation_period_subdivisions = static_cast<int>(incubation_period_seconds / 30.0);
    
    return incubation_period_subdivisions;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
int Inf_period_gamma() {
    // Infection period -> Gamma distribution (shape,scale)
    double infection_period_seconds = R::rgamma(shape_inf_g, scale_inf_g);
    int infection_period_subdivisions = static_cast<int>(infection_period_seconds / 30.0);
    
    return infection_period_subdivisions;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
int Inf_period_lognormal() {
    // Infection period -> Log-normal distribution (meanlog, sdlog)
    double infection_period_seconds = R::rlnorm(mInf, sd_inf_g);
    int infection_period_subdivisions = static_cast<int>(infection_period_seconds / 30.0);
    
    return infection_period_subdivisions;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
int Inf_period_uniform() {
    // Infection period -> Uniform distribution
    double infection_period_days = R::runif(3, 7);
    int infection_period_subdivisions = static_cast<int>( (infection_period_days * 3600*24) / 30.0);
    
    return infection_period_subdivisions;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
int Incub_period_uniform() {
    // Incub period -> Uniform distribution
    double incubation_period_days = R::runif(1, 3);
    int incubation_period_subdivisions = static_cast<int>( (incubation_period_days * 3600*24) / 30.0);
    
    return incubation_period_subdivisions;
};

//////////////////////////////////////////////
// [[Rcpp::export]]
int runif_int(int lower_value, int upper_value) {

    double p = 1 / double(upper_value-lower_value+1);
    IntegerVector to_sample_from;
    NumericVector probs;
    for (int i=lower_value; i < upper_value+1; i++) {
        to_sample_from.push_back(i);
        probs.push_back(p);
    }
    int out = sample(to_sample_from, 1, false, probs)[0];

    return out;
}

//////////////////////////////////////////////
// [[Rcpp::export]]
int Get_loc_j(
    const Rcpp::String& id,
    const Rcpp::IntegerVector& location_ti,
    const Rcpp::CharacterVector& ids_ti
) {
    // Search for individual's index
    int index_j = -1;
    for (int k = 0; k < ids_ti.size(); k++){
        if (id == ids_ti[k]){
            index_j = k;
            break;
        }
    }
    return location_ti[index_j];
}
