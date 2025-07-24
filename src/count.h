
#ifndef xxCOUNT_H
#define xxCOUNT_H

namespace count{


arma::vec numeric_Count(arma::vec Vec, int val=0){
  
  int n = Vec.size();
  
  arma::vec ansVec;
  arma::vec tmpVec(1);
  
  int newval=val;
  int oldval=val;
  int chkval=val;
  
  int nInc=0;
  int nCount=0;
  
  oldval = Vec(0);
  tmpVec[0] = oldval;
  
  ansVec.insert_rows(ansVec.n_rows, tmpVec);
  nInc++;
  
  for(int i=1;i<=n;i++){
    
    newval = Vec(i-1);
    if(newval!=oldval){
      
      nCount=0;
      for(int j=1;j<=nInc;j++){
        chkval = ansVec(j-1);
        if(chkval==newval){
          break;
        }else{
          nCount++;
        }
      }
      
      if(nInc==nCount){
        
        tmpVec[0] = newval;
        ansVec.insert_rows(ansVec.n_rows, tmpVec);
        nInc++;
        oldval = newval;
        
      }
      
    }
    
  }
  
  return ansVec;
  
  
}


Rcpp::StringVector str_Count(Rcpp::StringVector Vec, Rcpp::String val="dummy"){
  
  int n = Vec.size();
  
  Rcpp::StringVector ansVec;
  
  Rcpp::String newval=val;
  Rcpp::String oldval=val;
  Rcpp::String chkval=val;
  
  int nInc=0;
  int nCount=0;
  
  oldval = Vec(0);
  ansVec.push_back(oldval);
  nInc++;
  
  for(int i=1;i<=n;i++){
    
    newval = Vec(i-1);
    if(newval!=oldval){
      
      nCount=0;
      for(int j=1;j<=nInc;j++){
        chkval = ansVec(j-1);
        if(chkval==newval){
          break;
        }else{
          nCount++;
        }
      }
      
      if(nInc==nCount){
        
        ansVec.push_back(newval);
        nInc++;
        oldval = newval;
        
      }
      
    }
    
  }
  
  return ansVec;
  
  
}

}

#endif
