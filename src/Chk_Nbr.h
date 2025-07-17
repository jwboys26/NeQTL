

#ifndef CHK_NBR_H
#define CHK_NBR_H


void Revised_Get_qNk_Idx_Chr(arma::mat DM, arma::vec& ivec, arma::vec& Prob_q0, arma::vec& Prob_q1);
void cpp_Conv_Str2Numeric(Rcpp::StringVector strVec);



void Revised_Get_q_or_lNk_Idx_Chr(arma::mat DM, 
                                  arma::vec& ivec,
                                  arma::vec& Prob_q0,
                                  arma::vec& Prob_q1, arma::uvec colidx){
  
  arma::vec Total_GeneVec = DM.col(1);
  arma::vec Total_EnhVec = DM.col(2);
  arma::vec Total_qkVec =  DM.col(3);
  arma::vec Total_nIDVec = DM.col(4);
  
  
  int gID;
  
  //arma::vec GeneVec = count::numeric_Count(Total_GeneVec);
  
  arma::vec GeneVec = unique(Total_GeneVec);
  
  int nGene = GeneVec.size();
  int nEnh=0;
  
  arma::mat PrimeDM;
  arma::mat NeighborDM;
  
  int q1_Total=0;
  double q1_Neighbor=0;
  int q1_Prime=0;
  
  int enh_k;
  int nNeighbor = 0;
  
  double Cond_Prob = 0;
  
  
  
  arma::vec tmpVec(1);
  
  for(int j=1; j<= nGene; j++){
    gID = GeneVec[j-1];
    
    arma::uvec gidxVec = find(Total_GeneVec== gID);
    
    arma::mat DMj = DM.submat(gidxVec, colidx);
    
    
    arma::vec sub_qk = DMj.col(3);
    q1_Total = sum(sub_qk);
    
    arma::vec sub_enh = DMj.col(2);
    //arma::vec EnhVec= count::numeric_Count(sub_enh, q1_Total);
    arma::vec EnhVec= unique(sub_enh);
    nEnh = EnhVec.size();
    
    arma::vec sub_nID = DMj.col(4);
    
    
    if(nEnh>1){
      
      for(int k=1;k<=nEnh;k++){
        
        enh_k = EnhVec[k-1];
        
        arma::uvec enh_idx = find(sub_enh== enh_k);
        arma::mat PrimeDM = DMj.submat(enh_idx, colidx);
        
        arma::uvec enh_idx2 = find(sub_enh != enh_k);
        arma::mat NeighborDM = DMj.submat(enh_idx2, colidx);;
        
        nNeighbor = enh_idx2.size();
        
        arma::vec sub_qk_Neighbor = NeighborDM.col(3);
        q1_Neighbor = sum(sub_qk_Neighbor);
        
        Cond_Prob = q1_Neighbor/nNeighbor;
        tmpVec[0] = Cond_Prob;
        
        
        q1_Prime = q1_Total - q1_Neighbor;
        
        
        if(q1_Neighbor!=0){
          arma::vec sub_idvec = PrimeDM.col(4);
          ivec.insert_rows(ivec.n_rows, sub_idvec);
          
        }
        
        if(q1_Prime==0){
          Prob_q0.insert_rows(Prob_q0.n_rows, tmpVec);
          
        }else{
          Prob_q1.insert_rows(Prob_q1.n_rows, tmpVec);
          
        }
        
      }
    }
    
    
  }
  
  
  
}


void cpp_Conv_Str2Numeric(Rcpp::StringVector strVec){
  
  Rcpp::StringVector titleVec = count::str_Count(strVec);
  Rcpp::String tmpstr;
  
  int nlen = titleVec.size();
  
  arma::uvec numVec(nlen);
  
  for(int i=1;i<=nlen;i++){
    tmpstr = titleVec[i-1];
    
  }
  
}





#endif