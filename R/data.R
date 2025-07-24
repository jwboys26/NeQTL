
#'Dataset of original gold labels
#'
#'#####################
#'@format A list containing original references of gold labels (data frames that include gene id, chromosome, and gv position):
#' \describe{
#'   \item{battle_eqtl}{gold label data of eQTL from Battle (2014).}
#'   \item{muther_eqtl}{gold label data of eQTL from Multiple Tissue Human Expression Resource (Muther) project.}
#'   \item{geuvadis_eqtl}{gold label data of eQTL from Genetic European Variation in Disease (Geuvadis) project.}
#' }
#'
#'@docType data
#'@keywords datasets
#'@usage data(List_Gold_Label)
#'@seealso List_Modified_Gold_Label, Mod_Gold_Label()
#'@references
#'[1] Battle, A. et al. (2014) Characterizing the genetic basis of transcriptome diversity through RNA-sequencing of 922 individuals. Genome. Res. 24: 14-24
#'
#'@examples
#'
#'data(List_Gold_Label)
#'battle = List_Gold_Label$battle_eqtl
#'muther = List_Gold_Label$muther_eqtl
#'geuvadis = List_Gold_Label$geuvadis_eqtl
#'
"List_Gold_Label"





#'Example of valid input
#'
#'#####################
#'@format A list containing examples of valid and invalid input data:
#' \describe{
#'   \item{Valid_Input}{Example of a valid input data which contains all required varialbe.}
#'   \item{Invalid_Input1}{Example of an invalid input data.}
#'   \item{Invalid_Input2}{Another example of an invalid input data.}
#' }
#'
#'@docType data
#'@keywords datasets
#'@usage data(List_Example_Data)
#'@examples
#'
#'data(List_Example_Data)
#'Input_Data = List_Example_Data$Valid_Input
"List_Example_Data"



#'Result of accuracy measures obtained from NeQTL and Naive Bayesian
#'
#'#####################
#'@format A list of 100 (= number of percentiles) pairs of 4 accuracy measures (true positive, false positive, false negative, true negative) after adjusting the distance obtained by NeQTL and naive Bayesian:
#' \describe{
#'   \item{Accuracy_Mat1}{100x4 matrix of 4 accuracy measures of NeQTL over 100 iterations using the adjusted distance}
#'   \item{Accuracy_Mat2}{100x4 matrix of 4 accuracy measures of naive Bayesian over 100 iterations using the adjusted distance}
#'   
#' }
#'
#'@docType data
#'@keywords datasets
#'@usage data(List_Acc_Example)
#'@seealso \code{Compute_Accuracy}, \code{Draw_Adjusted_ROC}
#'
#'@examples
#'
#'data(List_Acc_Example)
#'Acc_NeQTL = List_Acc_Example$Accuracy_Mat1
#'Acc_Naive = List_Acc_Example$Accuracy_Mat2
#'
#'Draw_Adjusted_ROC(Acc_NeQTL, Acc_Naive, "NeQTL", "Naive")
#'
#'
"List_Acc_Example"



#'Dataset of modified gold labels for sample examples only. For your research, use List_Gold_Label. 
#'
#'#####################
#'@format A list containing modified gold labels (vectors of true or false):
#' \describe{
#'   \item{battle_eqtl}{gold label data of eQTL from Battle (2014).}
#'   \item{muther_eqtl}{gold label data of eQTL from Multiple Tissue Human Expression Resource (Muther) project.}
#'   \item{geuvadis_eqtl}{gold label data of eQTL from Genetic European Variation in Disease (Geuvadis) project.}
#' }
#'
#'@docType data
#'@keywords datasets
#'@usage data(List_Modified_Gold_Label)
#'@seealso List_Gold_Label, Mod_Gold_Label()
#'@references
#'[1] Battle, A. et al. (2014) Characterizing the genetic basis of transcriptome diversity through RNA-sequencing of 922 individuals. Genome. Res. 24: 14-24
#'@examples
#'
#'data(List_Modified_Gold_Label)
#'battle = List_Modified_Gold_Label$battle_eqtl
#'muther = List_Modified_Gold_Label$muther_eqtl
#'geuvadis = List_Modified_Gold_Label$geuvadis_eqtl
#'
#'
#'
"List_Modified_Gold_Label"



#'Sample result obtained from NeQTL
#'
#'#####################
#'@format A list of
#' \describe{
#'   \item{Prob_NeQTL_qk}{posterior probabilities for qk=1}
#'}
#'
#'@docType data
#'@keywords datasets
#'@usage data(List_NeQTL)
#'@seealso List_Naive_Bayes
#'@examples
#'
#'data(List_NeQTL)
#'Prob_NeQTL_qk = List_NeQTL$Prob_NeQTL_qk
#'
#'#### example of computing the auc using battle gold label
#'data(List_Modified_Gold_Label)
#'battle = List_Modified_Gold_Label$battle_eqtl
#'
#'df = cbind("gold_label"=battle, "NeQTL"=Prob_NeQTL_qk)
#'df = data.frame(df)
#'
#'library(ROCR)
#'#### NeQTL's auc
#'pred  = prediction(df$NeQTL, df$gold_label)
#'perf  = performance(pred, "tpr", "fpr")
#'auc_b = performance(pred, 'auc')
#'auc_b = unlist(slot(auc_b,"y.values"))
#'auc_b =  round(auc_b, digits=3)
#'plot(perf, xlim=c(0,1), ylim=c(0,1), col='red', main=paste("ROC curve of eQTL prediction"), lwd = 2)
#'print(auc_b)
#'
#'
#'
"List_NeQTL"


#'Sample result obtained from the application of naive Bayesian method, which doesn't use the neighboring relationship. 
#'
#'#####################
#'@format A list of
#' \describe{
#'   \item{Prob_Naive_qk}{posterior probabilities for qk=1}
#'}
#'
#'@docType data
#'@keywords datasets
#'@usage data(List_Naive_Bayes)
#'@seealso List_NeQTL
#'@examples
#'
#'data(List_Naive_Bayes)
#'Prob_Naive_qk = List_Naive_Bayes$Prob_Naive_qk
#'
#'#####Compare it with NeQTL.
#'
#'data(List_NeQTL)
#'Prob_NeQTL_qk = List_NeQTL$Prob_NeQTL_qk
#'
#'#### example of computing the auc using muther gold label
#'data(List_Modified_Gold_Label)
#'muther = List_Modified_Gold_Label$muther_eqtl
#'
#'df = cbind("gold_label"=muther, "NeQTL"=Prob_NeQTL_qk, "Naive"=Prob_Naive_qk)
#'df = data.frame(df)
#'
#'library(ROCR)
#'#### NeQTL's auc
#'pred  = prediction(df$NeQTL, df$gold_label)
#'perf  = performance(pred, "tpr", "fpr")
#'auc_b = performance(pred, 'auc')
#'auc_b = unlist(slot(auc_b,"y.values"))
#'auc_b =  round(auc_b, digits=3)
#'plot(perf, xlim=c(0,1), ylim=c(0,1), col='red', main=paste("ROC curve of eQTL prediction"), lwd = 2)
#'print(auc_b)
#'
#'
#'pred  = prediction(df$Naive, df$gold_label)
#'perf  = performance(pred, "tpr", "fpr")
#'auc_t = performance(pred, 'auc')
#'auc_t = unlist(slot(auc_t, "y.values"))
#'auc_t =  round(auc_t, digits=3)
#'par(new=TRUE)
#'plot(perf, xlim=c(0,1), ylim=c(0,1), col='blue', lwd = 2)
#'legend("bottomright", lwd= c(2,2), col = c("red", "blue"),
#'       legend = c(paste0("NeQTL, auc=", auc_b), paste0("Naive Bayesian, auc=", auc_t)))
#'
#'print(auc_t)
"List_Naive_Bayes"



#'
#'
#'#####################
#'@format A vector of distance between SNP and genes from real data used in this study.
#'@docType data
#'@keywords datasets
#'@usage data(DistVec)
#'@examples
#'
#'data(DistVec)
#'summary(DistVec)
#'
"DistVec"



