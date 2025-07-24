

#ifndef xxFI_H
#define xxFI_H


template <typename T1, typename T2>


void cpp_Find_Index(arma::uvec& iVec, T1 Vec, T2 val){
  
  int n = Vec.size();
  arma::uvec tmpVec(1);
  
  T2 tmp;
  
  for(int i=1;i<=n;i++){
    tmp = Vec[i-1];
    
    if(tmp==val){
      tmpVec[0] = i;
      iVec.insert_rows(iVec.n_rows, tmpVec);
      
    }
    
  }
  
}



arma::uvec numeric_Find_Index(arma::vec Vec, int val){
  
  arma::uvec uVec;
  cpp_Find_Index(uVec, Vec, val);
  
  return(uVec);
}



arma::uvec str_Find_Index(Rcpp::StringVector Vec, Rcpp::String val){
  
  arma::uvec uVec;
  cpp_Find_Index(uVec, Vec, val);
  
  return(uVec);
}




template <typename T1, typename T2>
void cpp_Find_Not_Index(arma::uvec& iVec, T1 Vec, T2 val){
  
  int n = Vec.size();
  arma::uvec tmpVec(1);
  
  T2 tmp;
  
  for(int i=1;i<=n;i++){
    tmp = Vec[i-1];
    
    if(tmp!=val){
      tmpVec[0] = i;
      iVec.insert_rows(iVec.n_rows, tmpVec);
      
    }
    
  }
  
}

arma::uvec numeric_Find_Not_Index(arma::vec Vec, int val){
  
  arma::uvec uVec;
  cpp_Find_Not_Index(uVec, Vec, val);
  return uVec;
  
}


arma::uvec str_Find_Not_Index(Rcpp::StringVector Vec, Rcpp::String val){
  
  arma::uvec uVec;
  cpp_Find_Not_Index(uVec, Vec, val);
  return uVec;
  
}




#endif
