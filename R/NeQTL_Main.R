



#'
#'##Performs Bayesian inference for eQTL analysis
#'@param Input_Data input data which information of gene, chromosome, enhancer, etc.
#'@param total_iter number of iterations. Default is set at 100.
#'@param qR threshold value for eQTL. Default is 0.1.
#'@param lR threshold value for chromatin interaction. Default is 0.1.
#'@param b_use_lNk option for incorporating a neighboring relationship for chromatin interaction. Default is FALSE.
#'@param bDisp  option for displaying the iteration. Default is FALSE.
#'@return list of the following values:
#'\describe{
#'\item{lk_index}{Indice of chromatin interactions.}
#'\item{qk_index}{Indice of eQTL SNPs.}
#'\item{qNk_index}{Indice of neighboring eQTLs.}
#'\item{prob_qk1_all}{Probability of a SNP being eQTL.}
#'\item{prob_lk1_all}{Probablity of a chromatin interaction.}
#'\item{prob_qNk_qk_0}{Conditional probability of neighboring eQTL SNPs given that it is not an eQTL SNP.}
#'\item{prob_qNk_qk_1}{Conditional probability of neighboring eQTL SNPs given that it is an eQTL SNP.}
#'\item{skq_qk0lk0}{eQTL score given no eQTL SNP (qk=0) and no chromatin interaction (lk=0).}
#'\item{skq_qk0lk1}{eQTL score given no eQTL SNP (qk=0) and chromatin interaction (lk=1).}
#'\item{skq_qk1lk0}{eQTL score given eQTL SNP (qk=1) and no chromatin interaction (lk=0).}
#'\item{skq_qk1lk1}{eQTL score given eQTL SNP (qk=1) and chromatin interaction (lk=1).}
#'\item{skl_qk0lk0}{chromatin interaction score given no eQTL SNP (qk=0) and no chromatin interaction (lk=0).}
#'\item{skl_qk0lk1}{chromatin interaction score given no eQTL SNP (qk=0) and chromatin interaction (lk=1).}
#'\item{skl_qk1lk0}{chromatin interaction score given eQTL SNP (qk=1) and no chromatin interaction (lk=0).}
#'\item{skl_qk1lk1}{chromatin interaction score given eQTL SNP (qk=1) and chromatin interaction (lk=1).}
#'
#'
#'}
#'
#'
#'
#'
#'@examples
#'####################
#'
#' ####Import or upload an input data. Here, an example data (Valid_Input) will be used.
#' data(List_Example_Data)
#' Valid_Input = List_Example_Data$Valid_Input
#' head(Valid_Input)
#' dim(Valid_Input)
#' Final_lst = Run_NeQTL2(Valid_Input, qR=0.1, lR=0.1, b_use_lNk=FALSE, total_iter=10, bDisp=FALSE)
#'#####################
#'
#'
#'@export
#'@importFrom Rcpp evalCpp
#'@useDynLib NeQTL

Run_NeQTL2 = function(Input_Data, qR = 0.1, lR = 0.1, b_use_lNk=FALSE, total_iter=100, bDisp=FALSE){

  K = nrow(Input_Data)

  ############################### Load coefficient of Logistic

  nCol=5
  Design_X = matrix(1, K, nCol)
  Design_X[,2] = abs(Input_Data$distance)
  Design_X[,3] = Input_Data$Gene_Activity
  Design_X[,4] = Input_Data$Enh_Activity
  Design_X[,5] = Input_Data$Correlation


  fitVal = exp(Design_X %*% logis_coef)/(1+exp(Design_X%*% logis_coef))
  Input_Data$skl = as.vector(fitVal)
  ###############################


  lInitial = Get_Quantile(Input_Data$skl, lR)
  whl = which(Input_Data$skl>=lInitial)
  Input_Data$lk = rep(0,times=K)
  Input_Data$lk[whl] = 1


  qInitial = Get_Quantile(Input_Data$skq, qR)
  whq = which(Input_Data$skq>=qInitial)
  Input_Data$qk = rep(0,times=K)
  Input_Data$qk[whq] = 1



  seqID = 1:K
  Input_Data = cbind(Input_Data, "nID"= seqID)
  Input_Data = cbind(Input_Data, "nEnh"= rep(0, times=K))
  Input_Data = cbind(Input_Data, "nGene"= rep(0, times=K))

  ##########
  ########## Changing Gene and Chr to numerical
  Vec = Input_Data$gene_id
  Vec2 = as.factor(Vec)
  Gene_nID = as.numeric(Vec2)

  Vec = Input_Data$chr
  Vec2 = as.factor(Vec)
  Chr_nID = as.numeric(Vec2)


  Input_Data = cbind(Input_Data, "Gene_nID"= Gene_nID)
  Input_Data = cbind(Input_Data, "Chr_nID"= Chr_nID)



  ###################### Make gene_id as numeric

  Enh_Categ = as.factor(Input_Data$enh_start)
  tbl = table(Enh_Categ)

  EnhVec = names(tbl)
  nNum_Enh = length(EnhVec)

  Gene_Categ = as.factor(Input_Data$gene_id)
  tbl = table(Gene_Categ)

  GeneVec = names(tbl)
  nNum_Gene = length(GeneVec)

  #####



  ################ Chromosome Vector
  Chr_Categ = as.factor(Input_Data$chr)
  tbl = table(Chr_Categ)

  ChrVec = names(tbl)
  nNum_Chr = length(ChrVec)

  ################


  lk_index = which(Input_Data$lk == 1)
  qk_index = which(Input_Data$qk == 1)

  ##################


  ### we have qNk_index below
  DM = Input_Data[which(Input_Data$lk==1), c("Chr_nID", "Gene_nID", "chr", "gene_id", "enh_start", "qk", "lk", "nID")]
  lst = Revised_Get_Cond_Prob(DM, b_qk=TRUE, bDisplay=FALSE)

  qNk_index=lst[[1]]
  prob_q1_Neighbor_q1_Primary=lst[[2]]
  prob_q1_Neighbor_q0_Primary=lst[[3]]

  rm(DM)
  rm(lst)


  if(b_use_lNk==TRUE){

    DM = Input_Data[, c("Chr_nID", "Gene_nID", "chr", "gene_id", "enh_start", "qk", "lk", "nID")]
    lst = Revised_Get_Cond_Prob(DM, b_qk=FALSE, bDisplay=FALSE)

    lNk_index=lst[[1]]
    prob_l1_Neighbor_l1_Primary=lst[[2]]
    prob_l1_Neighbor_l0_Primary=lst[[3]]

    rm(DM)
    rm(lst)

  }



  ### Skq socre for eQTL
  ### Skl score for chromatin interaction

  sq <- min(Input_Data[ , "skq" ])
  eq <- max(Input_Data[ , "skq" ])
  sl <- min(Input_Data[ , "skl" ])
  el <- max(Input_Data[ , "skl" ])



  ##############################
  for(iter in 1:total_iter) {

    print( paste0("start iteration: ", iter) )


    freq_table = Revised_MakeTable(K, lk_index, qk_index)



    # diagonal, row-sum, col-sum of freq table not zero
    if( ( sum(diag(freq_table)) != 0 ) | ( sum(c(apply(freq_table,1,sum), apply(freq_table,2,sum)) == 0 )==0)  ){

      lst = Revised_Get_skql(Input_Data$skl, Input_Data$skq, lk_index, qk_index)

      skq_qk0lk0 = lst[[1]]
      skq_qk0lk1 = lst[[2]]
      skq_qk1lk0 = lst[[3]]
      skq_qk1lk1 = lst[[4]]

      skl_qk0lk0 = lst[[5]]
      skl_qk0lk1 = lst[[6]]
      skl_qk1lk0 = lst[[7]]
      skl_qk1lk1 = lst[[8]]
      #####################
      rm(lst)

      # Prob(skq | qk=0, lk=0)
      kden             =  density(skq_qk0lk0, from=sq, to=eq, n=10000)
      prob_skq_qk0lk0  =  kden_prob(fkden=kden, fxval=Input_Data[ ,"skq"])


      kden             =  density(skq_qk0lk1, from=sq, to=eq, n=10000)
      prob_skq_qk0lk1  =  kden_prob(fkden=kden, fxval=Input_Data[ ,"skq"])


      kden             =  density(skq_qk1lk0, from=sq, to=eq, n=10000)
      prob_skq_qk1lk0  =  kden_prob(fkden=kden, fxval=Input_Data[ ,"skq"])



      kden             =  density(skq_qk1lk1, from=sq, to=eq, n=10000)
      prob_skq_qk1lk1  =  kden_prob(fkden=kden, fxval=Input_Data[ ,"skq"])


      # Prob(skl | qk=0, lk=0)
      kden             =  density(skl_qk0lk0, from=sl, to=el, n=10000)
      prob_skl_qk0lk0  =  kden_prob(fkden=kden, fxval=Input_Data[ ,"skl"])


      kden             =  density(skl_qk1lk0, from=sl, to=el, n=10000)
      prob_skl_qk1lk0  =  kden_prob(fkden=kden, fxval=Input_Data[ ,"skl"])



      kden             =  density(skl_qk0lk1, from=sl, to=el, n=10000)
      prob_skl_qk0lk1  =  kden_prob(fkden=kden, fxval=Input_Data[ ,"skl"])


      kden             =  density(skl_qk1lk1, from=sl, to=el, n=10000)
      prob_skl_qk1lk1  =  kden_prob(fkden=kden, fxval=Input_Data[ ,"skl"])



      # Prob(lk=0 | qk=0), scaler
      prob_lk0_qk0 = freq_table[1, 1] / sum(freq_table[, 1])
      prob_lk1_qk0 = freq_table[2, 1] / sum(freq_table[, 1])
      prob_lk0_qk1 = freq_table[1, 2] / sum(freq_table[, 2])
      prob_lk1_qk1 = freq_table[2, 2] / sum(freq_table[, 2])



      # Prob(qk=0 | lk=0), scaler
      prob_qk0_lk0 = freq_table[1, 1] / sum(freq_table[1, ])
      prob_qk1_lk0 = freq_table[1, 2] / sum(freq_table[1, ])
      prob_qk0_lk1 = freq_table[2, 1] / sum(freq_table[2, ])
      prob_qk1_lk1 = freq_table[2, 2] / sum(freq_table[2, ])


      tmp      = apply( freq_table, 2, sum)
      prob_qk0 = tmp[1] / K
      prob_qk1 = tmp[2] / K


      # Prob(lk=0)
      tmp      = apply( freq_table, 1, sum)
      prob_lk0 = tmp[1] / K
      prob_lk1 = tmp[2] / K


      #median(prob_q1_Neighbor_q1_Primary)

      prob_q0_Neighbor_q0_Primary=1-prob_q1_Neighbor_q0_Primary
      prob_q0_Neighbor_q1_Primary=1-prob_q1_Neighbor_q1_Primary



      prob_qNk_qk_1 = RFInvUVec(prob_q0_Neighbor_q1_Primary, K)  # rep(p0, times=K)   #
      prob_qNk_qk_1[qNk_index] = RFInvUVec(prob_q1_Neighbor_q1_Primary, length(qNk_index))   # p1 #


      prob_qNk_qk_0 = RFInvUVec(prob_q0_Neighbor_q0_Primary, K)   # rep(p0, times=K)   #
      prob_qNk_qk_0[qNk_index] = RFInvUVec(prob_q1_Neighbor_q0_Primary, length(qNk_index)) #p1


      ##########

      if(b_use_lNk==TRUE){

        prob_l0_Neighbor_l0_Primary=1-prob_l1_Neighbor_l0_Primary
        prob_l0_Neighbor_l1_Primary=1-prob_l1_Neighbor_l1_Primary


        p0 = mean(prob_l0_Neighbor_l1_Primary)
        prob_lNk_lk_1 = rep(p0, times=K) # RFInvUVec(prob_l0_Neighbor_l1_Primary, K)
        p1 = mean(prob_l1_Neighbor_l1_Primary)
        prob_lNk_lk_1[lNk_index] = p1 # RFInvUVec(prob_l1_Neighbor_l1_Primary, length(lNk_index))   # p1 #


        p0 = mean(prob_l0_Neighbor_l0_Primary)
        prob_lNk_lk_0 = rep(p0, times=K)   # RFInvUVec(prob_l0_Neighbor_l0_Primary, K)
        p1 = mean(prob_l1_Neighbor_l0_Primary)
        prob_lNk_lk_0[lNk_index] = p1 # RFInvUVec(prob_l1_Neighbor_l0_Primary, length(lNk_index)) #p1


      }else{
        prob_lNk_lk_1 = rep(0, times=K)
        prob_lNk_lk_0 = rep(0, times=K)
      }


      #####################################
      #####################################
      #sourceCpp(cppLibrary)

      lst2 = Gibbs_Sampler2(lk_index, qk_index, b_use_lNk,
                            prob_skq_qk1lk0, prob_skq_qk1lk1,
                            prob_skq_qk0lk0, prob_skq_qk0lk1,
                            prob_skl_qk1lk0, prob_skl_qk1lk1,
                            prob_skl_qk0lk0, prob_skl_qk0lk1,
                            prob_qNk_qk_0, prob_qNk_qk_1,
                            prob_lNk_lk_0, prob_lNk_lk_1,
                            prob_lk0_qk0, prob_lk0_qk1,
                            prob_lk1_qk0, prob_lk1_qk1,
                            prob_qk0_lk0, prob_qk0_lk1,
                            prob_qk1_lk0, prob_qk1_lk1,
                            prob_qk0, prob_qk1, prob_lk0, prob_lk1, K)


      prob_lk1_all = lst2[[1]]
      prob_qk1_all = lst2[[2]]
      lk_index = lst2[[3]]
      qk_index = lst2[[4]]

      rm(lst2)


      ############################## Update lk and qk

      #nlk = length(lk_index)
      #whn = sort(sample.int(nlk, floor(0.995*nlk)))
      #lk_index = lk_index[whn]


      Input_Data$lk = rep(0, times=K)
      Input_Data$lk[lk_index] = 1

      Input_Data$qk = rep(0, times=K)
      Input_Data$qk[qk_index] = 1


      ## Prob(qNk|qk=1 in primary enhancer)
      #####################################
      #########  For Conditional Probability

      DM = Input_Data[lk_index, c("Chr_nID", "Gene_nID", "chr", "gene_id", "enh_start", "qk", "lk", "nID")]
      Condlst = Revised_Get_Cond_Prob(DM, b_qk=TRUE, bDisplay=FALSE)


      qNk_index=Condlst[[1]]
      prob_q1_Neighbor_q1_Primary=Condlst[[2]]
      prob_q1_Neighbor_q0_Primary=Condlst[[3]]

      rm(DM)
      rm(Condlst)


      ##############################


      ## Prob(lNk|lk=1 in primary enhancer)
      #####################################
      #########  For Conditional Probability

      if(b_use_lNk==TRUE){

        DM = Input_Data[ , c("Chr_nID", "Gene_nID", "chr", "gene_id", "enh_start", "qk", "lk", "nID")]
        Condlst = Revised_Get_Cond_Prob(DM, b_qk=FALSE, bDisplay=FALSE)

        lNk_index=Condlst[[1]]
        prob_l1_Neighbor_l1_Primary=Condlst[[2]]
        prob_l1_Neighbor_l0_Primary=Condlst[[3]]

        rm(DM)
        rm(Condlst)


      }else{
        lNk_index=c(0)
      }



      ##############################

      if(bDisp==FALSE){
        tmpstr = paste(iter," iterations are finished.\n", sep="")
        print(tmpstr)
      }




    }else{
      break
    } ## diagonal, row-sum, col-sum of freq table not zero
  }



  Final_lst = list(
    lk_index=lk_index,
    qk_index=qk_index,
    qNk_index=qNk_index,
    lNk_index=lNk_index,

    prob_qk1_all=prob_qk1_all,
    prob_lk1_all=prob_lk1_all,
    prob_qNk_qk_0=prob_qNk_qk_0,
    prob_qNk_qk_1=prob_qNk_qk_1,
    prob_q1_Neighbor_q1_Primary=prob_q1_Neighbor_q1_Primary,
    prob_q1_Neighbor_q0_Primary=prob_q1_Neighbor_q0_Primary,
    skq_qk0lk0=skq_qk0lk0,
    skq_qk0lk1=skq_qk0lk1,
    skq_qk1lk0=skq_qk1lk0,
    skq_qk1lk1=skq_qk1lk1,
    skl_qk0lk0=skl_qk0lk0,
    skl_qk0lk1=skl_qk0lk1,
    skl_qk1lk0=skl_qk1lk0,
    skl_qk1lk1=skl_qk1lk1 )

  return(Final_lst)


}



