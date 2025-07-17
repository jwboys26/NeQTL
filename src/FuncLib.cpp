#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

#include "count.h"
#include "Chk_Nbr.h"
#include "fi.h"
#include "CreateDM.h"


int GetBernoulli(double prob){
  
  int out=0;
  
  double val = R::rbinom(1, prob);
  out = floor(val);
  return out;
  
}




// [[Rcpp::export]]
arma::mat Revised_MakeTable(unsigned int n, arma::uvec vec_Index_lk, arma::uvec vec_Index_qk){
  

  unsigned int nl_Index = vec_Index_lk.n_elem;
  unsigned int nq_Index = vec_Index_qk.n_elem;
  
  //printf("nl_Index is %d", nl_Index);
  
  arma::mat out(2,2);
  out.zeros();
  
  int a11=0;
  int a12=0;
  int a21=0;
  int a22=0;
  
  int lk=0, qk=0;
  unsigned int nbl=1, nbq=1;
  
  
  for(unsigned int i=1;i<=n;i++){
    
    if(nl_Index==0){
      lk=0;
    }else{
      
      lk=0;
      
      for(unsigned int j=nbl;j<=nl_Index;j++){
        
        if(i==vec_Index_lk[j-1]){
          lk = 1;
          nbl += 1;
          break;
        }else{
          lk = 0;
          break;
        }
      }
      
      
      
    }
    
    
    if(nq_Index==0){
      qk=0;
    }else{
      
      qk=0;
      for(unsigned int k=nbq;k<=nq_Index;k++){
        
        if(i== vec_Index_qk[k-1]){
          qk = 1;
          nbq += 1;
          break;
        }else{
          qk = 0;
          break;
        }
        
      }
      
      
    }
    
    
    if((lk==0) & (qk==0) ){
      a11 += 1;
    }else if((lk==0) & (qk==1) ){
      a12 += 1;
      //printf("lk, qk, i are %d, %d, %d\n", lk, qk, i);
    }else if((lk==1) & (qk==0) ){
      a21 += 1;
    }else{
      a22 += 1;
    }
    
  }
  
  out(0,0)=a11;
  out(0,1)=a12;
  out(1,0)=a21;
  out(1,1)=a22;
  
  return out;
  
}





//[[Rcpp::export]]
List Revised_Get_skql(arma::vec vec_skl, arma::vec vec_skq, 
                      arma::vec vec_Index_lk, arma::vec vec_Index_qk){
  
  int n = vec_skl.n_elem;
  int nl_Index = vec_Index_lk.n_elem;
  int nq_Index = vec_Index_qk.n_elem;
  
  List lst(8);
  
  double qk=0, lk=0;
  double skl, skq=0;
  
  
  int nq_qk00=1, nq_qk01=1, nq_qk10=1, nq_qk11 = 1;
  int nl_qk00=1, nl_qk01=1, nl_qk10=1, nl_qk11 = 1;
  
  arma::vec skq_qk0lk0(n);
  arma::vec skq_qk0lk1(n);
  arma::vec skq_qk1lk0(n);
  arma::vec skq_qk1lk1(n); 
  
  arma::vec skl_qk0lk0(n);
  arma::vec skl_qk0lk1(n);
  arma::vec skl_qk1lk0(n);
  arma::vec skl_qk1lk1(n); 
  
  
  int nbl = 1;
  int nbq = 1;
  
  for(int i=1;i<=n;i++){
    
    //vec_Index_lk
    
    if(nl_Index==0){
      lk=0;      
    }else{
      lk=0;
      for(int j=nbl;j<=nl_Index;j++){
        
        if(i==vec_Index_lk[j-1]){
          lk = 1;
          nbl += 1;
          break;
        }else{
          lk = 0;
          break;
        }
      }
    }
    
    if(nq_Index==0){
      qk=0;
    }else{
      
      qk=0;
      for(int k=nbq;k<=nq_Index;k++){
        
        if(i== vec_Index_qk[k-1]){
          qk = 1;
          nbq += 1;
          break;
        }else{
          qk = 0;
          break;
        }
        
      }
      
    }
    
    
    
    skl = vec_skl[i-1];
    skq = vec_skq[i-1];
    
    if(qk==0){
      if(lk==0){
        //skl = UD(i, nCol_skl);
        skl_qk0lk0[nl_qk00-1]=skl;
        nl_qk00 += 1;
        
        //skq = UD(i, nCol_skq);
        skq_qk0lk0[nq_qk00-1]=skq;
        nq_qk00 += 1;
        
        
      }else{    ////// qk=0, lk=1
        
        //skl = UD(i, nCol_skl);
        skl_qk0lk1[nl_qk01-1]=skl;
        nl_qk01 += 1;
        
        //skq = UD(i, nCol_skq);
        skq_qk0lk1[nq_qk01-1]=skq;
        nq_qk01 += 1;
        
      }
    }else{
      
      if(lk==0){    /// qk=1 lk=0
        //skl = UD(i, nCol_skl);
        skl_qk1lk0[nl_qk10-1]=skl;
        nl_qk10 += 1;
        
        //skq = UD(i, nCol_skq);
        skq_qk1lk0[nq_qk10-1]=skq;
        nq_qk10 += 1;
        
        
      }else{    ////// qk=1, lk=1
        
        //skl = UD(i, nCol_skl);
        skl_qk1lk1[nl_qk11-1]=skl;
        nl_qk11 += 1;
        
        //skq = UD(i, nCol_skq);
        skq_qk1lk1[nq_qk11-1]=skq;
        nq_qk11 += 1;
        
      }
      
    }
    
    
  }
  
  
  
  
  
  nq_qk00 = nq_qk00-1;
  nq_qk01 = nq_qk01-1;
  nq_qk10 = nq_qk10-1;
  nq_qk11 = nq_qk11-1;
  
  nl_qk00 = nl_qk00-1;
  nl_qk10 = nl_qk10-1;
  nl_qk01 = nl_qk01-1;
  nl_qk11 = nl_qk11-1;
  
  
  lst[0] = skq_qk0lk0.subvec(0, nq_qk00-1);
  lst[1] = skq_qk0lk1.subvec(0, nq_qk01-1);
  lst[2] = skq_qk1lk0.subvec(0, nq_qk10-1);
  lst[3] = skq_qk1lk1.subvec(0, nq_qk11-1);
  
  lst[4] = skl_qk0lk0.subvec(0, nl_qk00-1);
  lst[5] = skl_qk0lk1.subvec(0, nl_qk01-1);
  lst[6] = skl_qk1lk0.subvec(0, nl_qk10-1);
  lst[7] = skl_qk1lk1.subvec(0, nl_qk11-1);
  
  return lst;
  
}



//[[Rcpp::export]]
List Revised_Get_Cond_Prob(const Rcpp::DataFrame DM, bool b_qk=true,  bool bDisplay=false, int n_Elems=1e+5){
  
  arma::vec Total_ChrVec = DM["Chr_nID"];
  arma::vec Total_GeneVec = DM["Gene_nID"];
  arma::vec Total_EnhVec = DM["enh_start"];
  arma::vec Total_qkVec = DM["qk"];
  arma::vec Total_lkVec = DM["lk"];
  
  arma::vec Total_nIDVec = DM["nID"];
  
  int K = Total_nIDVec.size();
  
  arma::mat Num_Mat(K, 5);
  Num_Mat.col(0) = Total_ChrVec;
  Num_Mat.col(1) = Total_GeneVec;
  Num_Mat.col(2) = Total_EnhVec;
  
  if(b_qk){
    Num_Mat.col(3) = Total_qkVec;
  }else{
    Num_Mat.col(3) = Total_lkVec;
  }
  Num_Mat.col(4) = Total_nIDVec;
  
  
  int chr;
  
  arma::vec IndVec;
  
  arma::vec Prob_q_or_l1_Neighbor_q_or_l1_Primary; 
  arma::vec Prob_q_or_l1_Neighbor_q_or_l0_Primary;
  
  
  Rcpp::DataFrame DMc;
  arma::uvec chr_subidx;
  
  arma::uvec colidx(5);
  colidx[0]=0;
  colidx[1]=1;
  colidx[2]=2;
  colidx[3]=3;
  colidx[4]=4;
  
  
  //arma::vec ChrVec = count::numeric_Count(Total_ChrVec);
  arma::vec ChrVec = unique(Total_ChrVec);
  int nChr = ChrVec.size();
  
  int N_Size = 0;
  
  arma::uvec rand_index;
  arma::mat sub_Mat;
  
  for(int i=1;i<=nChr; i++){
    chr = ChrVec[i-1];
    chr_subidx = find(Total_ChrVec==chr);
    
    N_Size = chr_subidx.size();
    
    
    if(N_Size>n_Elems){
      rand_index = arma::conv_to<arma::uvec>::from(arma::randi(n_Elems, arma::distr_param(0, N_Size-1)));
      sub_Mat = Num_Mat.submat(chr_subidx.elem( sort(rand_index) ), colidx);
    }else{
      sub_Mat = Num_Mat.submat(chr_subidx, colidx);
      
    }
    
    arma::vec idxVec;
    arma::vec Prob_q_or_l1; 
    arma::vec Prob_q_or_l0;
    
    Revised_Get_q_or_lNk_Idx_Chr(sub_Mat, idxVec, Prob_q_or_l0, Prob_q_or_l1, colidx);
    
    IndVec.insert_rows(IndVec.n_rows, idxVec);
    Prob_q_or_l1_Neighbor_q_or_l0_Primary.insert_rows(Prob_q_or_l1_Neighbor_q_or_l0_Primary.n_rows, Prob_q_or_l0);
    Prob_q_or_l1_Neighbor_q_or_l1_Primary.insert_rows(Prob_q_or_l1_Neighbor_q_or_l1_Primary.n_rows, Prob_q_or_l1);
    
    if(bDisplay){
      Rcout<<"Finished chr"<<chr<<"\n";
    }
    
  }
  
  
  
  List lst(3);
  lst[0] = IndVec;
  lst[1] = Prob_q_or_l1_Neighbor_q_or_l1_Primary;
  lst[2] = Prob_q_or_l1_Neighbor_q_or_l0_Primary;
  
  return lst;
  
}






//[[Rcpp::export]] 
List Get_Final_Result_Neighbor(arma::vec lk_index, arma::vec qk_index, arma::vec qNk_index,
                                  arma::vec prob_skq_qk1lk0, arma::vec prob_skq_qk1lk1,
                                  arma::vec prob_skq_qk0lk0, arma::vec prob_skq_qk0lk1,
                                  arma::vec prob_skl_qk1lk0, arma::vec prob_skl_qk1lk1,
                                  arma::vec prob_skl_qk0lk0, arma::vec prob_skl_qk0lk1,
                                  arma::vec prob_qNk_qk_0, arma::vec prob_qNk_qk_1,
                                  double prob_lk0_qk0, double prob_lk0_qk1,
                                  double prob_lk1_qk0, double prob_lk1_qk1,
                                  double prob_qk0_lk0, double prob_qk0_lk1,
                                  double prob_qk1_lk0, double prob_qk1_lk1,
                                  double prob_qk0, double prob_qk1, double prob_lk0, double prob_lk1, int K){
  
  
  int nl_Index = lk_index.n_elem;
  int nq_Index = qk_index.n_elem;
  
  
  List lst(5);
  
  arma::vec prob_lk1_all(K);
  arma::vec prob_qk1_all(K);
  
  arma::vec new_lk_index(K);
  arma::vec new_qk_index(K);
  
  int nbl=1, nbq=1;
  int lk=0, qk=0;
  
  int nL=1, nQ=1;
  
  
  double qk1_num=0, qk1_dnom=0, skq_qk1_lk=0, lk_qk1=0, skq_qk0_lk=0, lk_qk0=0;
  double lk1_num=0, lk1_dnom=0, skl_lk1_qk=0, qk_lk1=0, skl_lk0_qk=0, qk_lk0=0;
  
  double skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
  
  double qNk_q1=0, qNk_q0=0;
  
  int BerVal = 0;
  long double tmp=0;
  
  for(int i=1;i<=K;i++){
    
    if(nl_Index==0){
      lk=0;
    }else{
      
      lk=0;
      
      for(int j=nbl;j<=nl_Index;j++){
        
        if(i==lk_index[j-1]){
          lk = 1;
          nbl += 1;
          break;
        }else{
          lk = 0;
          break;
        }
      }
      
    }
    
    
    if(lk==0){
      
      skq_qk1_lk = prob_skq_qk1lk0[i-1];
      skq_qk0_lk = prob_skq_qk0lk0[i-1];
      
      lk_qk1 = prob_lk0_qk1;
      lk_qk0 = prob_lk0_qk0;
      
      
      
      ///////////////////////   skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
      skl_lk_qk1 = prob_skl_qk1lk0[i-1];
      skl_lk_qk0 = prob_skl_qk0lk0[i-1];
      ////////////////////////
      
      
      
      
    }else{
      
      skq_qk1_lk = prob_skq_qk1lk1[i-1];
      skq_qk0_lk = prob_skq_qk0lk1[i-1];
      
      lk_qk1 = prob_lk1_qk1;
      lk_qk0 = prob_lk1_qk0;
      
      ///////////////////////   skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
      skl_lk_qk1 = prob_skl_qk1lk1[i-1];
      skl_lk_qk0 = prob_skl_qk0lk1[i-1];
      ////////////////////////
      
      
      
      
    }
    
    
    if(nq_Index==0){
      qk=0;
    }else{
      
      qk=0;
      for(int k=nbq;k<=nq_Index;k++){
        
        if(i== qk_index[k-1]){
          qk = 1;
          nbq += 1;
          break;
        }else{
          qk = 0;
          break;
        }
        
      }
      
      
    }
    
    
    
    
    
    if(qk==0){
      
      skl_lk1_qk = prob_skl_qk0lk1[i-1];
      skl_lk0_qk = prob_skl_qk0lk0[i-1];
      
      
      qk_lk1 = prob_qk0_lk1;
      qk_lk0 = prob_qk0_lk0;
      
      ///////////////////////   skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
      skq_qk_lk1 = prob_skq_qk0lk1[i-1];
      skq_qk_lk0 = prob_skq_qk0lk0[i-1];
      ////////////////////////
      
      qNk_q0 = prob_qNk_qk_0[i-1];
      
    }else{
      
      skl_lk1_qk = prob_skl_qk1lk1[i-1];
      skl_lk0_qk = prob_skl_qk1lk0[i-1];
      
      qk_lk1 = prob_qk1_lk1;
      qk_lk0 = prob_qk1_lk0;
      
      ///////////////////////   skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
      skq_qk_lk1 = prob_skq_qk1lk1[i-1];
      skq_qk_lk0 = prob_skq_qk1lk0[i-1];
      ////////////////////////
      
      qNk_q1 = prob_qNk_qk_1[i-1];
      
    }
    
    
    
    
    if(skq_qk1_lk!=0){
      
      qk1_num=skl_lk_qk1*skq_qk1_lk*lk_qk1*qNk_q1*prob_qk1;
      qk1_dnom = skl_lk_qk0 * skq_qk0_lk * lk_qk0 *qNk_q0* prob_qk0 + qk1_num;
      
      if(qk1_dnom==0){
        tmp=0;
      }else{
        tmp = qk1_num/qk1_dnom ;
      }
      
      prob_qk1_all[i-1] = tmp;
      
      BerVal = GetBernoulli(tmp);
      if(BerVal==1){
        new_qk_index[nQ-1] = i;
        nQ += 1;
      }
      
    }else{
      prob_qk1_all[i-1] = 0;
    }
    
    //printf("finish Q1\n");
    
    
    
    if(skl_lk1_qk!=0){
      
      lk1_num= skq_qk_lk1* skl_lk1_qk*qk_lk1*prob_lk1;
      lk1_dnom = skq_qk_lk0* skl_lk0_qk * qk_lk0 * prob_lk0 + lk1_num;
      
      
      if(lk1_dnom==0){
        tmp=0;
      }else{
        tmp = lk1_num/lk1_dnom;
      }
      
      prob_lk1_all[i-1] = tmp;
      
      //printf("finish , %.3f, %.3f\n", lk1_dnom, tmp);
      
      BerVal = GetBernoulli(tmp);
      if(BerVal==1){
        new_lk_index[nL-1] = i;
        nL += 1;
      }
      
    }else{
      prob_lk1_all[i-1] = 0;
    }
    
    
    
    
  }
  
  nQ -= 1;
  nL -= 1;
  
  
  
  lst[0]=prob_lk1_all;
  lst[1]=prob_qk1_all;
  
  lst[2]=new_lk_index.subvec(0, nL-1);
  lst[3]=new_qk_index.subvec(0, nQ-1);
  
  
  return lst;
  
}












//[[Rcpp::export]] 
List Gibbs_Sampler(arma::vec lk_index, arma::vec qk_index,
                       arma::vec prob_skq_qk1lk0, arma::vec prob_skq_qk1lk1,
                       arma::vec prob_skq_qk0lk0, arma::vec prob_skq_qk0lk1,
                       arma::vec prob_skl_qk1lk0, arma::vec prob_skl_qk1lk1,
                       arma::vec prob_skl_qk0lk0, arma::vec prob_skl_qk0lk1,
                       arma::vec prob_qNk_qk_0, arma::vec prob_qNk_qk_1,
                       arma::vec prob_lNk_lk_0, arma::vec prob_lNk_lk_1,
                       double prob_lk0_qk0, double prob_lk0_qk1,
                       double prob_lk1_qk0, double prob_lk1_qk1,
                       double prob_qk0_lk0, double prob_qk0_lk1,
                       double prob_qk1_lk0, double prob_qk1_lk1,
                       double prob_qk0, double prob_qk1, double prob_lk0, double prob_lk1, int K){
  
  
  int nl_Index = lk_index.n_elem;
  int nq_Index = qk_index.n_elem;
  
  
  List lst(5);
  
  arma::vec prob_lk1_all(K);
  arma::vec prob_qk1_all(K);
  
  arma::vec new_lk_index(K);
  arma::vec new_qk_index(K);
  
  int nbl=1, nbq=1;
  int lk=0, qk=0;
  
  int nL=1, nQ=1;
  
  
  double qk1_num=0, qk1_dnom=0, skq_qk1_lk=0, lk_qk1=0, skq_qk0_lk=0, lk_qk0=0;
  double lk1_num=0, lk1_dnom=0, skl_lk1_qk=0, qk_lk1=0, skl_lk0_qk=0, qk_lk0=0;
  
  double skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
  
  double qNk_q1=0, qNk_q0=0;
  double lNk_l1=0, lNk_l0=0;
  
  int BerVal = 0;
  long double tmp=0;
  
  for(int i=1;i<=K;i++){
    
    if(nl_Index==0){
      lk=0;
    }else{
      
      lk=0;
      
      for(int j=nbl;j<=nl_Index;j++){
        
        if(i==lk_index[j-1]){
          lk = 1;
          nbl += 1;
          break;
        }else{
          lk = 0;
          break;
        }
      }
      
    }
    
    
    if(lk==0){
      
      skq_qk1_lk = prob_skq_qk1lk0[i-1];
      skq_qk0_lk = prob_skq_qk0lk0[i-1];
      
      lk_qk1 = prob_lk0_qk1;
      lk_qk0 = prob_lk0_qk0;
      
      
      
      ///////////////////////   skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
      skl_lk_qk1 = prob_skl_qk1lk0[i-1];
      skl_lk_qk0 = prob_skl_qk0lk0[i-1];
      ////////////////////////
      
      lNk_l0 = prob_lNk_lk_0[i-1];
      
      
      
    }else{
      
      skq_qk1_lk = prob_skq_qk1lk1[i-1];
      skq_qk0_lk = prob_skq_qk0lk1[i-1];
      
      lk_qk1 = prob_lk1_qk1;
      lk_qk0 = prob_lk1_qk0;
      
      ///////////////////////   skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
      skl_lk_qk1 = prob_skl_qk1lk1[i-1];
      skl_lk_qk0 = prob_skl_qk0lk1[i-1];
      ////////////////////////
      
      
      lNk_l1 = prob_lNk_lk_1[i-1];
      
      
    }
    
    
    if(nq_Index==0){
      qk=0;
    }else{
      
      qk=0;
      for(int k=nbq;k<=nq_Index;k++){
        
        if(i== qk_index[k-1]){
          qk = 1;
          nbq += 1;
          break;
        }else{
          qk = 0;
          break;
        }
        
      }
      
      
    }
    
    
    
    
    
    if(qk==0){
      
      skl_lk1_qk = prob_skl_qk0lk1[i-1];
      skl_lk0_qk = prob_skl_qk0lk0[i-1];
      
      
      qk_lk1 = prob_qk0_lk1;
      qk_lk0 = prob_qk0_lk0;
      
      ///////////////////////   skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
      skq_qk_lk1 = prob_skq_qk0lk1[i-1];
      skq_qk_lk0 = prob_skq_qk0lk0[i-1];
      ////////////////////////
      
      qNk_q0 = prob_qNk_qk_0[i-1];
      
    }else{
      
      skl_lk1_qk = prob_skl_qk1lk1[i-1];
      skl_lk0_qk = prob_skl_qk1lk0[i-1];
      
      qk_lk1 = prob_qk1_lk1;
      qk_lk0 = prob_qk1_lk0;
      
      ///////////////////////   skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
      skq_qk_lk1 = prob_skq_qk1lk1[i-1];
      skq_qk_lk0 = prob_skq_qk1lk0[i-1];
      ////////////////////////
      
      qNk_q1 = prob_qNk_qk_1[i-1];
      
    }
    
    
    
    
    if(skq_qk1_lk!=0){
      
      qk1_num=skl_lk_qk1*skq_qk1_lk*lk_qk1*qNk_q1*prob_qk1;
      qk1_dnom = skl_lk_qk0 * skq_qk0_lk * lk_qk0 *qNk_q0* prob_qk0 + qk1_num;
      
      if(qk1_dnom==0){
        tmp=0;
      }else{
        tmp = qk1_num/qk1_dnom ;
      }
      
      prob_qk1_all[i-1] = tmp;
      
      BerVal = GetBernoulli(tmp);
      if(BerVal==1){
        new_qk_index[nQ-1] = i;
        nQ += 1;
      }
      
    }else{
      prob_qk1_all[i-1] = 0;
    }
    
    //printf("finish Q1\n");
    
    
    
    if(skl_lk1_qk!=0){
      
      lk1_num= skq_qk_lk1* skl_lk1_qk*qk_lk1*lNk_l1*prob_lk1;
      lk1_dnom = skq_qk_lk0* skl_lk0_qk * qk_lk0 *lNk_l0* prob_lk0 + lk1_num;
      
      
      if(lk1_dnom==0){
        tmp=0;
      }else{
        tmp = lk1_num/lk1_dnom;
      }
      
      prob_lk1_all[i-1] = tmp;
      
      //printf("finish , %.3f, %.3f\n", lk1_dnom, tmp);
      
      BerVal = GetBernoulli(tmp);
      if(BerVal==1){
        new_lk_index[nL-1] = i;
        nL += 1;
      }
      
    }else{
      prob_lk1_all[i-1] = 0;
    }
    
    
    
    
  }
  
  nQ -= 1;
  nL -= 1;
  
  
  
  lst[0]=prob_lk1_all;
  lst[1]=prob_qk1_all;
  
  lst[2]=new_lk_index.subvec(0, nL-1);
  lst[3]=new_qk_index.subvec(0, nQ-1);
  
  
  return lst;
  
}




//[[Rcpp::export]] 
List Gibbs_Sampler2(arma::vec lk_index, arma::vec qk_index, bool b_use_lNk,
                   arma::vec prob_skq_qk1lk0, arma::vec prob_skq_qk1lk1,
                   arma::vec prob_skq_qk0lk0, arma::vec prob_skq_qk0lk1,
                   arma::vec prob_skl_qk1lk0, arma::vec prob_skl_qk1lk1,
                   arma::vec prob_skl_qk0lk0, arma::vec prob_skl_qk0lk1,
                   arma::vec prob_qNk_qk_0, arma::vec prob_qNk_qk_1,
                   arma::vec prob_lNk_lk_0, arma::vec prob_lNk_lk_1,
                   double prob_lk0_qk0, double prob_lk0_qk1,
                   double prob_lk1_qk0, double prob_lk1_qk1,
                   double prob_qk0_lk0, double prob_qk0_lk1,
                   double prob_qk1_lk0, double prob_qk1_lk1,
                   double prob_qk0, double prob_qk1, double prob_lk0, double prob_lk1, int K){
  
  
  int nl_Index = lk_index.n_elem;
  int nq_Index = qk_index.n_elem;
  
  
  List lst(5);
  
  arma::vec prob_lk1_all(K);
  arma::vec prob_qk1_all(K);
  
  arma::vec new_lk_index(K);
  arma::vec new_qk_index(K);
  
  int nbl=1, nbq=1;
  int lk=0, qk=0;
  
  int nL=1, nQ=1;
  
  
  double qk1_num=0, qk1_dnom=0, skq_qk1_lk=0, lk_qk1=0, skq_qk0_lk=0, lk_qk0=0;
  double lk1_num=0, lk1_dnom=0, skl_lk1_qk=0, qk_lk1=0, skl_lk0_qk=0, qk_lk0=0;
  
  double skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
  
  double qNk_q1=0, qNk_q0=0;
  double lNk_l1=0, lNk_l0=0;
  
  int BerVal = 0;
  long double tmp=0;
  
  for(int i=1;i<=K;i++){
    
    if(nl_Index==0){
      lk=0;
    }else{
      
      lk=0;
      
      for(int j=nbl;j<=nl_Index;j++){
        
        if(i==lk_index[j-1]){
          lk = 1;
          nbl += 1;
          break;
        }else{
          lk = 0;
          break;
        }
      }
      
    }
    
    
    if(lk==0){
      
      skq_qk1_lk = prob_skq_qk1lk0[i-1];
      skq_qk0_lk = prob_skq_qk0lk0[i-1];
      
      lk_qk1 = prob_lk0_qk1;
      lk_qk0 = prob_lk0_qk0;
      
      
      
      ///////////////////////   skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
      skl_lk_qk1 = prob_skl_qk1lk0[i-1];
      skl_lk_qk0 = prob_skl_qk0lk0[i-1];
      ////////////////////////
      
      lNk_l0 = prob_lNk_lk_0[i-1];
      
      
      
    }else{
      
      skq_qk1_lk = prob_skq_qk1lk1[i-1];
      skq_qk0_lk = prob_skq_qk0lk1[i-1];
      
      lk_qk1 = prob_lk1_qk1;
      lk_qk0 = prob_lk1_qk0;
      
      ///////////////////////   skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
      skl_lk_qk1 = prob_skl_qk1lk1[i-1];
      skl_lk_qk0 = prob_skl_qk0lk1[i-1];
      ////////////////////////
      
      
      lNk_l1 = prob_lNk_lk_1[i-1];
      
      
    }
    
    
    if(nq_Index==0){
      qk=0;
    }else{
      
      qk=0;
      for(int k=nbq;k<=nq_Index;k++){
        
        if(i== qk_index[k-1]){
          qk = 1;
          nbq += 1;
          break;
        }else{
          qk = 0;
          break;
        }
        
      }
      
      
    }
    
    
    
    
    
    if(qk==0){
      
      skl_lk1_qk = prob_skl_qk0lk1[i-1];
      skl_lk0_qk = prob_skl_qk0lk0[i-1];
      
      
      qk_lk1 = prob_qk0_lk1;
      qk_lk0 = prob_qk0_lk0;
      
      ///////////////////////   skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
      skq_qk_lk1 = prob_skq_qk0lk1[i-1];
      skq_qk_lk0 = prob_skq_qk0lk0[i-1];
      ////////////////////////
      
      qNk_q0 = prob_qNk_qk_0[i-1];
      
    }else{
      
      skl_lk1_qk = prob_skl_qk1lk1[i-1];
      skl_lk0_qk = prob_skl_qk1lk0[i-1];
      
      qk_lk1 = prob_qk1_lk1;
      qk_lk0 = prob_qk1_lk0;
      
      ///////////////////////   skl_lk_qk0=0, skl_lk_qk1=0, skq_qk_lk0=0, skq_qk_lk1=0;
      skq_qk_lk1 = prob_skq_qk1lk1[i-1];
      skq_qk_lk0 = prob_skq_qk1lk0[i-1];
      ////////////////////////
      
      qNk_q1 = prob_qNk_qk_1[i-1];
      
    }
    
    
    
    
    if(skq_qk1_lk!=0){
      
      qk1_num=skl_lk_qk1*skq_qk1_lk*lk_qk1*qNk_q1*prob_qk1;
      qk1_dnom = skl_lk_qk0 * skq_qk0_lk * lk_qk0 *qNk_q0* prob_qk0 + qk1_num;
      
      if(qk1_dnom==0){
        tmp=0;
      }else{
        tmp = qk1_num/qk1_dnom ;
      }
      
      prob_qk1_all[i-1] = tmp;
      
      BerVal = GetBernoulli(tmp);
      if(BerVal==1){
        new_qk_index[nQ-1] = i;
        nQ += 1;
      }
      
    }else{
      prob_qk1_all[i-1] = 0;
    }
    
    //printf("finish Q1\n");
    
    
    
    if(skl_lk1_qk!=0){
      
      if(b_use_lNk == TRUE){
        lk1_num= skq_qk_lk1* skl_lk1_qk*qk_lk1*lNk_l1*prob_lk1;
        lk1_dnom = skq_qk_lk0* skl_lk0_qk * qk_lk0 *lNk_l0* prob_lk0 + lk1_num;
        
      }else{
        lk1_num= skq_qk_lk1* skl_lk1_qk*qk_lk1*prob_lk1;
        lk1_dnom = skq_qk_lk0* skl_lk0_qk * qk_lk0 * prob_lk0 + lk1_num;
        
      }
      
      
      if(lk1_dnom==0){
        tmp=0;
      }else{
        tmp = lk1_num/lk1_dnom;
      }
      
      prob_lk1_all[i-1] = tmp;
      
      //printf("finish , %.3f, %.3f\n", lk1_dnom, tmp);
      
      BerVal = GetBernoulli(tmp);
      if(BerVal==1){
        new_lk_index[nL-1] = i;
        nL += 1;
      }
      
    }else{
      prob_lk1_all[i-1] = 0;
    }
    
    
    
    
  }
  
  nQ -= 1;
  nL -= 1;
  
  
  
  lst[0]=prob_lk1_all;
  lst[1]=prob_qk1_all;
  
  lst[2]=new_lk_index.subvec(0, nL-1);
  lst[3]=new_qk_index.subvec(0, nQ-1);
  
  
  return lst;
  
}











