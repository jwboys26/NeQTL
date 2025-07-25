% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RFuncLib_For_Distance.R
\name{Balancing_Dist}
\alias{Balancing_Dist}
\title{### Adjusting a pair of distances (between eQTL SNPs and genes) obtained by method 1 and method 2.}
\usage{
Balancing_Dist(
  Original_Dist,
  Dist1,
  Dist2,
  prob_qk1,
  prob_qk2,
  ProbVal1,
  ProbVal2,
  gold_label
)
}
\arguments{
\item{Original_Dist}{a vector of ditances of all SNPs and genes}

\item{Dist1}{the distance vector obtained from the method 1}

\item{Dist2}{the distance vector obtained from the method 2}

\item{prob_qk1}{the vector of posterior probabilities obtained from the method 1}

\item{prob_qk2}{the vector of posterior probabilities obtained from the method 2}

\item{ProbVal1}{the threshold posterior probability obtained from the method 1}

\item{ProbVal2}{the threshold posterior probability obtained from the method 2}

\item{gold_label}{a gold label dataset}
}
\value{
a list of
\describe{
\item{Adj_Dist1}{Adjusted distance for the method 1}
\item{Adj_Dist2}{Adjusted distance for the method 2}
\item{Adj_prob_qk1}{Adjusted posterior probobailities for the method 1}
\item{Adj_prob_qk2}{Adjusted posterior probobailities for the method 2}
}
}
\description{
### Adjusting a pair of distances (between eQTL SNPs and genes) obtained by method 1 and method 2.
}
\examples{


data(List_NeQTL)
Prob_NeQTL_qk = List_NeQTL$Prob_NeQTL_qk


data(DistVec)   ### Use DistVec for this example only. 
percent = 0.05
lst = Dist_SNP_Gene(Prob_NeQTL_qk, DistVec, percent)
Dist_SNP_Gene_NeQTL = lst$Dist_eQTL_Gene
Threshold_NeQTL = lst$ProbVal

data(List_Naive_Bayes)
Prob_Naive_qk = List_Naive_Bayes$Prob_Naive_qk

lst = Dist_SNP_Gene(Prob_Naive_qk, DistVec, percent)
Dist_SNP_Gene_Naive = lst$Dist_eQTL_Gene
Threshold_Naive = lst$ProbVal

data(List_Gold_Label)
gold_label = List_Gold_Label$muther_eqtl     ### Use Muther eQTL data

lst = Balancing_Dist(DistVec, Dist_SNP_Gene_NeQTL, Dist_SNP_Gene_Naive, 
Prob_NeQTL_qk, Prob_Naive_qk, Threshold_NeQTL, Threshold_Naive, gold_label)

Adj_Dist_NeQTL = lst$Adj_Dist1
Adj_Dist_Naive = lst$Adj_Dist2


mainstr="Adjusted distances between eQTL SNP's and genes"
legendstr=c("Blue: NeQLT", "Pink: Naive")

GetHist2(Adj_Dist_NeQTL, Adj_Dist_Naive, mainstr=mainstr, legendstr=legendstr)



}
\seealso{
GetHist2(), Dist_SNP_Gene()
}
