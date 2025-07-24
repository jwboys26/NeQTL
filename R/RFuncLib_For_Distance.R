

#'## Get two histograms of distances between genes and eQTLs obtained from two methods.
#'@param Dist_a a vector of distance between genes and eQTL obtained from the method a (denoted by blue in the histogram)
#'@param Dist_b a vector of distance between genes and eQTL obtained from the method b (denoted by pink in the histogram)
#'@param mainstr a title of the histogram
#'@param legendstr a vector of legend titles for the method a and method b. Default is c("Blue: Method A", "Pink: Method B") 
#'
#'@return a plot of two histogram. 
#'@seealso Dist_SNP_Gene()
#'@examples
#'####################
#'
#' 
#'data(List_Naive_Bayes)
#'Prob_Naive_qk = List_Naive_Bayes$Prob_Naive_qk
#'
#'#####Compare it with NeQTL.
#'
#'data(List_NeQTL)
#'Prob_NeQTL_qk = List_NeQTL$Prob_NeQTL_qk
#'#####################
#'
#'
#'
#'data(DistVec)   ### Use DistVec for this example only. 
#'percent = 0.05
#'lst = Dist_SNP_Gene(Prob_NeQTL_qk, DistVec, percent)
#'Dist_SNP_Gene_NeQTL = lst$Dist_eQTL_Gene
#'Threshold_NeQTL = lst$ProbVal
#'lst = Dist_SNP_Gene(Prob_Naive_qk, DistVec, percent)
#'Dist_SNP_Gene_Naive = lst$Dist_eQTL_Gene
#'Threshold_Naive = lst$ProbVal
#'
#'
#'mainstr="Distances between eQTL SNP's and genes"
#'legendstr=c("Blue: NeQLT", "Pink: Naive")
#'GetHist2(Dist_SNP_Gene_NeQTL, Dist_SNP_Gene_Naive, mainstr=mainstr, legendstr=legendstr)
#'
#'@export

GetHist2 = function(Dist_a, Dist_b, mainstr="Plot of two histograms", legendstr=c("Blue: Method A", "Pink: Method B")){
  
  xlabstr="Distance"
  
  minVal = min(c(Dist_a, Dist_b))
  maxVal = max(c(Dist_a, Dist_b))
  
  ax = pretty(minVal:maxVal, n=100)
  
  hista = hist(Dist_a, breaks=ax, plot=FALSE)
  histb = hist(Dist_b, breaks=ax, plot=FALSE)
  
  
  c1 <- rgb(173,216,230, maxColorValue = 255, alpha = 80, names = "lt.blue")
  c2 <- rgb(255,192,203, maxColorValue = 255, alpha = 80, names = "lt.pink")
  
  
  plot(hista, col=c1, main=mainstr, xlab=xlabstr)
  plot(histb, col=c2, add=TRUE)
  
  legend("topright", col = c(c1, c2),
         legend = legendstr)
  
}



#'# Get a distance between putative eQTL SNPs and associated genes. 
#'@param prob_qk a vector of posterior probabilities of SNP to be an eQTL
#'@param DistVec a vector of distance between genes and SNPs.
#'@param percent a percentage usded for thresholding.
#'
#'@return a list of
#'\describe{
#'\item{Dist_eQTL_Gene}{a vector of distance between putative eQTL SNPs and genes}
#'\item{ProbVal}{a probability used to determine eQTL, which is a (1-percent) percentile of prob_qk}
#'}
#'
#'@seealso GetHist2()
#'@examples
#'####################
#'
#' 
#'#'
#'data(List_NeQTL)
#'Prob_NeQTL_qk = List_NeQTL$Prob_NeQTL_qk
#'#####################
#'
#'
#'data(DistVec)   ### Use DistVec for this example only. 
#'percent = 0.05
#'lst = Dist_SNP_Gene(Prob_NeQTL_qk, DistVec, percent)
#'Dist_SNP_Gene_NeQTL = lst$Dist_eQTL_Gene
#'Threshold_NeQTL = lst$ProbVal
#'
#'hist(Dist_SNP_Gene_NeQTL, main="Distribution of distances between eQTLs and gene")
#'
#'
#'@export


Dist_SNP_Gene = function(prob_qk, DistVec, percent){
  
  
  ProbVec = prob_qk
  ProbQuan = quantile(ProbVec, (1-percent))
  ProbVal = ProbQuan[[1]][1]
  
  whProb = which(ProbVec>=ProbVal)
  Dist_eQTL_Gene = DistVec[whProb]
  
  
  lst= list(Dist_eQTL_Gene=Dist_eQTL_Gene, ProbVal=ProbVal)
  
  return(lst)
}


#'### Adjusting a pair of distances (between eQTL SNPs and genes) obtained by method 1 and method 2. 
#'@param Original_Dist a vector of ditances of all SNPs and genes
#'@param Dist1 the distance vector obtained from the method 1 
#'@param Dist2 the distance vector obtained from the method 2 
#'@param prob_qk1 the vector of posterior probabilities obtained from the method 1
#'@param prob_qk2 the vector of posterior probabilities obtained from the method 2
#'@param ProbVal1 the threshold posterior probability obtained from the method 1
#'@param ProbVal2 the threshold posterior probability obtained from the method 2
#'@param gold_label a gold label dataset
#'
#'@return a list of
#'\describe{
#'\item{Adj_Dist1}{Adjusted distance for the method 1}
#'\item{Adj_Dist2}{Adjusted distance for the method 2}
#'\item{Adj_prob_qk1}{Adjusted posterior probobailities for the method 1}
#'\item{Adj_prob_qk2}{Adjusted posterior probobailities for the method 2}
#'}
#'
#'@seealso GetHist2(), Dist_SNP_Gene()
#'@examples
#' 
#'
#'data(List_NeQTL)
#'Prob_NeQTL_qk = List_NeQTL$Prob_NeQTL_qk
#'
#'
#'data(DistVec)   ### Use DistVec for this example only. 
#'percent = 0.05
#'lst = Dist_SNP_Gene(Prob_NeQTL_qk, DistVec, percent)
#'Dist_SNP_Gene_NeQTL = lst$Dist_eQTL_Gene
#'Threshold_NeQTL = lst$ProbVal
#'
#'data(List_Naive_Bayes)
#'Prob_Naive_qk = List_Naive_Bayes$Prob_Naive_qk
#'
#'lst = Dist_SNP_Gene(Prob_Naive_qk, DistVec, percent)
#'Dist_SNP_Gene_Naive = lst$Dist_eQTL_Gene
#'Threshold_Naive = lst$ProbVal
#'
#'data(List_Gold_Label)
#'gold_label = List_Gold_Label$muther_eqtl     ### Use Muther eQTL data
#'
#'lst = Balancing_Dist(DistVec, Dist_SNP_Gene_NeQTL, Dist_SNP_Gene_Naive, 
#'Prob_NeQTL_qk, Prob_Naive_qk, Threshold_NeQTL, Threshold_Naive, gold_label)
#'
#'Adj_Dist_NeQTL = lst$Adj_Dist1
#'Adj_Dist_Naive = lst$Adj_Dist2
#' 
#'
#'mainstr="Adjusted distances between eQTL SNP's and genes"
#'legendstr=c("Blue: NeQLT", "Pink: Naive")
#'
#'GetHist2(Adj_Dist_NeQTL, Adj_Dist_Naive, mainstr=mainstr, legendstr=legendstr)
#'
#'
#'
#'@export

Balancing_Dist = function(Original_Dist, Dist1, Dist2, prob_qk1, prob_qk2, ProbVal1, ProbVal2, gold_label){
  
  lst_ax = Get_ax_dVec_cVec(Dist1, Dist2)
  
  ax   = lst_ax$ax
  dVec = lst_ax$dVec
  cVec = lst_ax$cVec
  
  nc = length(cVec)
  
  
  Adj_Dist1 = Original_Dist
  Adj_Dist2 = Original_Dist
  
  Adj_gold_label1 = gold_label
  Adj_gold_label2 = gold_label
  
  
  for(i in 1:nc){
    
    if(i==1){
      starting_dist = 0
    }else{
      starting_dist = dVec[i-1]
    }
    ending_dist = dVec[i]
    
    
    qk_index_a_i = which((prob_qk1>=ProbVal1 ) & (Adj_Dist1>starting_dist) &(Adj_Dist1<=ending_dist))
    qk_index_b_i = which((prob_qk2>=ProbVal2 ) & (Adj_Dist2>starting_dist) &(Adj_Dist2<=ending_dist))
    
    nq_a_i = length(qk_index_a_i)
    nq_b_i = length(qk_index_b_i)
    
    counta = nq_a_i
    countb = nq_b_i
    
    maxcount = max(counta, countb)
    
    if(counta>countb){
      
      nNew_q = countb
      if(nq_a_i>(counta-countb) ){
        removed_index = sort( sample.int(nq_a_i, (counta-countb) ) )
        removed_qk_index_a_i = qk_index_a_i[removed_index]
      }else{
        removed_qk_index_a_i = qk_index_a_i
      }
      
      Adj_Dist1 = Adj_Dist1[-removed_qk_index_a_i]
      prob_qk1 = prob_qk1[-removed_qk_index_a_i]
      Adj_gold_label1 = Adj_gold_label1[-removed_qk_index_a_i]
      
    }else if(counta<countb){
      
      if( nq_b_i>(countb-counta) ){
        removed_index = sort( sample.int(nq_b_i, (countb-counta) ) )
        removed_qk_index_b_i = qk_index_b_i[removed_index]
        
      }else{
        removed_qk_index_b_i = qk_index_b_i
      }
      
      
      Adj_Dist2 = Adj_Dist2[-removed_qk_index_b_i]
      prob_qk2 = prob_qk2[-removed_qk_index_b_i]
      Adj_gold_label2 = Adj_gold_label2[-removed_qk_index_b_i]
      
      
    }
    
  }
  
  
  
  whProb1 = which(prob_qk1>=ProbVal1)
  AdjDist_SNP_Gene1 = Adj_Dist1[whProb1]
  
  whProb2 = which(prob_qk2>=ProbVal2)
  AdjDist_SNP_Gene2 = Adj_Dist2[whProb2]
  
  
  lst = list(Adj_Dist1=AdjDist_SNP_Gene1, Adj_Dist2=AdjDist_SNP_Gene2, Adj_prob_qk1 = prob_qk1, Adj_prob_qk2 = prob_qk2, Adj_gold_label1=Adj_gold_label1, Adj_gold_label2=Adj_gold_label2)
  return(lst)
  
}










GetTrueFalse = function(GroundVec, PredVec){

  nLen = length(GroundVec)

  TPIndex = which((GroundVec==1)&(PredVec==1))
  nTP = length(TPIndex)

  FPIndex = which((GroundVec==0)&(PredVec==1))
  nFP = length(FPIndex)

  FNIndex= which((GroundVec==1)&(PredVec==0))
  nFN = length(FNIndex)

  TNIndex= which((GroundVec==0)&(PredVec==0))
  nTN = length(TNIndex)

  ans = c(nTP, nFP, nFN, nTN)
  return(ans)

}











Get_ax_dVec_cVec = function(Dist_a, Dist_b){

  minVal = min(c(Dist_a, Dist_b))
  maxVal = max(c(Dist_a, Dist_b))

  ax = pretty(minVal:maxVal, n=100)

  hista = hist(Dist_a, breaks=ax, plot=FALSE)
  histb = hist(Dist_b, breaks=ax, plot=FALSE)

  cVec = hista$count
  nc = length(cVec)

  d0 = hista$breaks

  dVec = rep(0, times=length(d0)-1)
  dVec = d0[2: length(d0)]

  lst=list(ax=ax, dVec=dVec, cVec=cVec)

  return(lst)
}











######## for  adjusted distance


#'## Get two accurcy measures using the adjusted distances between genes and eQTLs obtained from two methods.
#'@param DistVec a vector of the original distance between genes and eQTL 
#'@param Post_Prob1 a vector of posterior probabilities of eQTL obtained from method 1
#'@param Post_Prob2 a vector of posterior probabilities of eQTL obtained from method 2
#'@param gold_label a vector of gold labels
#'@param bDisp an option to display the intermediate result.
#'
#'@return a list of two 100 x 4 matrices, each of which contains true positive, false positive, false negative, true negative rates after using the adjusted distance obtained by two methods. 
#'@seealso \code{List_Acc_Example}, \code{Draw_Adjusted_ROC()}
#'@examples
#'####################
#'
#'
#'data(List_Gold_Label)
#'gold_label = List_Modified_Gold_Label$muther_eqtl     ### Use Muther eQTL data
#'
#'data(DistVec)   ### Use DistVec for this example only. 
#'
#'data(List_NeQTL)
#'Prob_NeQTL_qk = List_NeQTL$Prob_NeQTL_qk
#'
#'data(List_Naive_Bayes)
#'Prob_Naive_qk = List_Naive_Bayes$Prob_Naive_qk
#'
#'List_Acc = Compute_Accuracy(DistVec, Prob_NeQTL_qk, Prob_Naive_qk, gold_label, bDisp=FALSE)
#'
#'Acc_NeQTL = List_Acc$Accuracy_Mat1
#'Acc_Naive = List_Acc$Accuracy_Mat2
#'
#'head(Acc_NeQTL)
#'head(Acc_Naive)
#'
#'
#'##### Draw the adjusted ROC
#'Draw_Adjusted_ROC(Acc_NeQTL, Acc_Naive, "NeQTL", "Naive")
#'
#'@export



Compute_Accuracy = function(DistVec, Post_Prob1, Post_Prob2, gold_label, bDisp=FALSE){


  pRVec = seq(0.01, 1, 0.01)
  nR = length(pRVec)

  Accuracy_Mat1 = matrix(0, nR, 4)
  Accuracy_Mat2 = matrix(0, nR, 4)




  for(i in 1:nR){

    percent = pRVec[i]

    lst = Dist_SNP_Gene(Post_Prob1, DistVec, percent)
    Dist_SNP_Gene1 = lst$Dist_eQTL_Gene
    Threshold1 = lst$ProbVal


    lst2 = Dist_SNP_Gene(Post_Prob2, DistVec, percent)
    Dist_SNP_Gene2 = lst2$Dist_eQTL_Gene
    Threshold2 = lst2$ProbVal

    ###################################

    lst = Balancing_Dist(DistVec, Dist_SNP_Gene1, Dist_SNP_Gene2, Post_Prob1, Post_Prob2, Threshold1, Threshold2, gold_label)
    
    
    #
    Adj_Post_Prob1 = lst$Adj_prob_qk1
    Adj_Post_Prob2 = lst$Adj_prob_qk2

    Adj_gold_label1 = lst$Adj_gold_label1
    Adj_gold_label2 = lst$Adj_gold_label2


    AdjDist_SNP_Gene1 = lst$Adj_Dist1

    AdjDist_SNP_Gene2 = lst$Adj_Dist2

    ##################################


    Gra1 = Adj_gold_label1
    Pra1 = rep(0, times=length(Adj_Post_Prob1))
    Pra1[which(Adj_Post_Prob1>=Threshold1)] = 1
    lst1 = GetTrueFalse(Gra1, Pra1)

    nTP1=lst1[1]
    nFP1=lst1[2]
    nFN1=lst1[3]
    nTN1=lst1[4]

    Accuracy_Mat1[i, 1] = nTP1
    Accuracy_Mat1[i, 2] = nFP1
    Accuracy_Mat1[i, 3] = nFN1
    Accuracy_Mat1[i, 4] = nTN1



    Gra2 = Adj_gold_label2
    Pra2 = rep(0, times=length(Adj_Post_Prob2))
    Pra2[which(Adj_Post_Prob2>=Threshold2)] = 1
    lst2 = GetTrueFalse(Gra2, Pra2)

    nTP2=lst2[1]
    nFP2=lst2[2]
    nFN2=lst2[3]
    nTN2=lst2[4]

    Accuracy_Mat2[i, 1] = nTP2
    Accuracy_Mat2[i, 2] = nFP2
    Accuracy_Mat2[i, 3] = nFN2
    Accuracy_Mat2[i, 4] = nTN2


    if(bDisp==TRUE){
      tmpstr = paste(i, "% is completed.", sep=""  )
      print(tmpstr)
    }



  }

  lst = list(Accuracy_Mat1=Accuracy_Mat1, Accuracy_Mat2=Accuracy_Mat2)

  return(lst)
}




#'## Draw the ROC curves of two methods using the adjusted distances.
#'@param Accuracy_Mat1 a 100 x 4 matrix of the adjusted accuracy measures of method 1.
#'@param Accuracy_Mat2 a 100 x 4 matrix of the adjusted accuracy measures of method 2.
#'@param method_name1 the name of method 1. This will be expressed in the plot.
#'@param method_name2 the name of method 2. This will be expressed in the plot.
#'@param bDisp an option to display the intermediate result.
#'
#'@return two ROC curves using the adjusted distance obtained by two methods. 
#'@seealso \code{List_Acc_Example}, \code{Compute_Accuracy}
#'@examples
#'####################
#'
#'
#'## Here we will use already-adjusted accuracy measures we found in our study. Thus, we skip the process of adjusting measures.
#'data(List_Acc_Example)
#'
#'Acc_NeQTL = List_Acc_Example$Accuracy_Mat1
#'Acc_Naive = List_Acc_Example$Accuracy_Mat2
#'
#'head(Acc_NeQTL)
#'head(Acc_Naive)
#'
#'#'Draw_Adjusted_ROC(Acc_NeQTL, Acc_Naive, "NeQTL", "Naive")
#'
#'
#'@export



Draw_Adjusted_ROC = function(Accuracy_Mat1, Accuracy_Mat2,
                      method_name1, method_name2){
  
  Mat1 = Accuracy_Mat1
  Mat2 = Accuracy_Mat2
  nLen = nrow(Mat1)
  
  FPVec1 = rep(0, times=nLen)
  TPVec1 = rep(0, times=nLen)
  
  FPVec2 = rep(0, times=nLen)
  TPVec2 = rep(0, times=nLen)
  
  
  AUC1=0
  AUC2=0
  
  for(i in 1:nLen){
    
    nTP=Mat1[i,1]
    nFP=Mat1[i,2]
    nFN=Mat1[i,3]
    nTN=Mat1[i,4]
    
    TPVec1[i] = nTP/(nTP+nFN)
    FPVec1[i] = nFP/(nFP+nTN)
    
    
    nTP=Mat2[i,1]
    nFP=Mat2[i,2]
    nFN=Mat2[i,3]
    nTN=Mat2[i,4]
    
    TPVec2[i] = nTP/(nTP+nFN)
    FPVec2[i] = nFP/(nFP+nTN)
    
    
  }
  
  
  TPVec1 = c(0, TPVec1)
  FPVec1 = c(0, FPVec1)
  
  TPVec2 = c(0, TPVec2)
  FPVec2 = c(0, FPVec2)
  
  AUC1 = round(Compute_AUC(Mat1), digits=3)
  AUC2 = round(Compute_AUC(Mat2), digits=3)
  
  
  pchvec = c(3,4,8)
  ltyvec = c(1,2,4)
  colvec = c("lightsalmon", "dodgerblue3")
  lwdval=1.5
  lwdval2=1.5
  cexval=1.2
  
  
  
  mainstr=paste("ROC curve")
  
  
  plot(FPVec1,TPVec1, xlim=c(0,1), ylim=c(0,1), lty=ltyvec[1], col=colvec[1],
       xlab = "FPR", ylab="TPR",
       main= mainstr , lwd = 2, type="l")
  par(new=TRUE)
  plot(FPVec2, TPVec2,  xlim=c(0,1), ylim=c(0,1), lty=ltyvec[2],  col=colvec[2],
       xlab = "FPR", ylab="TPR", lwd = 2, type="l")
  
  tmp1 = paste(method_name1, " (AUC=", AUC1, ")", sep="")
  tmp2 = paste(method_name2, " (AUC=", AUC2, ")", sep="")
  
  legend("bottomright", col = colvec,
         legend = c(tmp1, tmp2 ),
         lty=c(ltyvec[1], ltyvec[2]), 
         cex=cexval, lwd=lwdval, box.lty=2, box.lwd=2)
  abline(0,1)
  
  
}






