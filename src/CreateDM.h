
#ifndef CREATEDM_H
#define CREATEDM_H


Rcpp::StringVector Get_sub_strVec(Rcpp::StringVector strVec, arma::uvec idxVec){

  int nID = idxVec.size();

  arma::uword idx=0;

  //Rcpp::StringVector sub_strVec(nID);
  Rcpp::StringVector sub_strVec;
  Rcpp::String tmpStr;

  for(int i=1;i<=nID;i++){
    idx = idxVec[i-1];
    //sub_strVec[i-1] = strVec[idx-1];
    sub_strVec.push_back(strVec[idx-1]);

  }

  return sub_strVec;
}


Rcpp::DataFrame Create_subDM(Rcpp::DataFrame DM, arma::uvec idxVec){


  arma::uvec g_Vec = DM["Gene_nID"];
  arma::uvec chr_Vec = DM["Chr_nID"];

  arma::uvec id_Vec = DM["nID"];
  arma::uvec qk_Vec = DM["qk"];
  arma::uvec enh_Vec = DM["enh_start"];

  arma::uvec idxVec2 = idxVec-1;

  arma::uvec sub_g_Vec = g_Vec.elem(idxVec2);
  arma::uvec sub_chr_Vec = chr_Vec.elem(idxVec2);

  arma::uvec sub_id_Vec = id_Vec.elem(idxVec2);
  arma::uvec sub_qk_Vec = qk_Vec.elem(idxVec2);
  arma::uvec sub_enh_Vec = enh_Vec.elem(idxVec2);


  Rcpp::DataFrame subDM =

    Rcpp::DataFrame::create(Rcpp::Named("Gene_nID") = sub_g_Vec,
                            Rcpp::Named("Chr_nID") = sub_chr_Vec,
                            Rcpp::Named("nID") = sub_id_Vec,
                            Rcpp::Named("qk") = sub_qk_Vec,
                            Rcpp::Named("enh_start") = sub_enh_Vec);



  return subDM;

}


Rcpp::DataFrame Create_subDM2(Rcpp::DataFrame DM, arma::uvec idxVec){



  arma::uvec id_Vec = DM["nID"];
  arma::uvec qk_Vec = DM["qk"];
  arma::uvec enh_Vec = DM["enh_start"];

  arma::uvec idxVec2 = idxVec-1;


  arma::uvec sub_id_Vec = id_Vec.elem(idxVec2);
  arma::uvec sub_qk_Vec = qk_Vec.elem(idxVec2);
  arma::uvec sub_enh_Vec = enh_Vec.elem(idxVec2);


  Rcpp::DataFrame subDM =

    Rcpp::DataFrame::create(Rcpp::Named("nID") = sub_id_Vec,
                            Rcpp::Named("qk") = sub_qk_Vec,
                            Rcpp::Named("enh_start") = sub_enh_Vec);



  return subDM;

}




#endif
