
## Create global interaction list
cppFunction('Rcpp::DataFrame get_global_interaction(int t, Rcpp::DataFrame interactions) {

  CharacterVector from = interactions["from"];
  CharacterVector from_status = interactions["from_status"];
  CharacterVector to = interactions["to"];
  CharacterVector to_status = interactions["to_status"];
  NumericVector date_posix = interactions["date_posix"];
  NumericVector length = interactions["length"];
  
  IntegerVector index;
  
  for (int i=0; i<date_posix.size(); i++) {
    if (date_posix[i] <= (t-1)*30 && date_posix[i] + length[i] >= (t)*30 ) {
      index.push_back(i);
    }
  }
  
  from = from[index];
  to = to[index];
  from_status = from_status[index];
  to_status = to_status[index];
  
  Rcpp::DataFrame interactions_t;
  if (index.size() == 0) {
      interactions_t = DataFrame::create(
      _["from"] = CharacterVector({}), 
      _["to"] = CharacterVector({}),
      _["time"] = NumericVector({}), 
      _["from_status"] = CharacterVector({}),
      _["to_status"] = CharacterVector({})
      );
      
  } else {
    interactions_t = DataFrame::create(
      _["from"] = from, 
      _["to"] = to, 
      _["time"] = (index.size(), t), 
      _["from_status"] = from_status,
      _["to_status"] = to_status 
      );
  }
  
  return interactions_t;
}')

# Create clusters list
cppFunction('std::list<std::vector<std::string>> get_clusters(Rcpp::DataFrame df) {

  std::vector<std::string> from = df["from"];
  std::vector<std::string> to = df["to"];
  std::vector<std::vector<std::string>> clusters_list;
  std::list<std::vector<std::string>> clusters_list_out;
  
  if (df.nrows() != 0) {
    for (int i=0; i<df.nrow(); i++) {
      if (i==0) {
        clusters_list.push_back({from[i], to[i]});
        
      } else {
        
        int absent = 0;
        for (int j=0; j<clusters_list.size(); j++) {
          std::vector<std::string> current_cluster = clusters_list[j];
          bool from_included = std::find(current_cluster.begin(), current_cluster.end(), from[i]) != current_cluster.end();
          bool to_included = std::find(current_cluster.begin(), current_cluster.end(), to[i]) != current_cluster.end();
          
          if (from_included && !to_included) clusters_list[j].push_back(to[i]);
          if (!from_included && to_included) clusters_list[j].push_back(from[i]);
          if (!from_included && !to_included) absent++;
        }
        
        if (absent == clusters_list.size()) clusters_list.push_back({from[i], to[i]});
      }
    }
    
    // Convert to list of vector
    for (int i=0;i<clusters_list.size();i++) clusters_list_out.push_back(clusters_list[i]);
  } 
  
  return clusters_list_out;
}')

# Create location matrix
cppFunction('CharacterMatrix get_global_location(
std::vector<std::vector<std::vector<std::string>>> clusters, 
Rcpp::DataFrame admission,
int n_subdivisions,
Rcpp::DataFrame rooms
) {
  
  // Initialize 3D location array
  size_t dimX = n_subdivisions;
  size_t dimY = admission.nrows();
  CharacterMatrix location_t(dimX,dimY);
  for (int x=0; x<dimX; x++) {
    for (int y=0; y<dimY; y++){
      location_t(x,y) = NA_STRING;
    }
  }
  
  // Names vectors to facilitate data accession
  IntegerVector ind_index ;
  for (int i=0; i<admission.nrows(); i++) ind_index.push_back(i);
  ind_index.names() = admission["id"];
  CharacterVector ind_status = admission["cat"] ;
  ind_status.names() = admission["id"];
  CharacterVector all_ids = admission["id"];
  CharacterVector id_rooms = rooms["id_room"];
  id_rooms.names() = rooms["id"];
  
  // Loop over clusters 
  for (int t=0; t<n_subdivisions; t++) {
  
    std::vector<std::vector<std::string>> clusters_t = clusters[t];
    
    for (auto &clusters_t_n: clusters_t) {
      
      // Get composition of the cluster
      int n_medical = 0;
      int n_patient = 0;
      CharacterVector patients_in_cluster;
      CharacterVector patient_rooms;
      
      for (auto &ind: clusters_t_n) {
        if (strcmp(ind_status[ind], "Medical") == 0) n_medical++;
        if (strcmp(ind_status[ind], "Patient") == 0) {
          n_patient++;
          patients_in_cluster.push_back(ind);
          patient_rooms.push_back(id_rooms[ind]);
        }
      }
      
      // Get location according to cluster composition 
      String cluster_location("");
      if (n_medical == clusters_t_n.size()) cluster_location.push_back(id_rooms["M-O"]); // Medical office 
      if (n_patient == 0 && n_medical != clusters_t_n.size()) cluster_location.push_back(id_rooms["NS"]); // Nursing station
      if (n_patient == 1) cluster_location+=patient_rooms[0]; // Patient room
      if (n_patient > 1) cluster_location.push_back(id_rooms["C"]); // Corridor
      
      // Assign cluster location
      for (auto &ind: clusters_t_n) {
        location_t(t, ind_index[ind]) = cluster_location;
      }
    }
  }
  
  colnames(location_t) = all_ids;
  return location_t;
}')