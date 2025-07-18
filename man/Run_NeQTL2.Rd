% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NeQTL_Main.R
\name{Run_NeQTL2}
\alias{Run_NeQTL2}
\title{##Performs Bayesian inference for eQTL analysis}
\usage{
Run_NeQTL2(
  Input_Data,
  qR = 0.1,
  lR = 0.1,
  b_use_lNk = FALSE,
  total_iter = 100,
  bDisp = FALSE
)
}
\arguments{
\item{Input_Data}{input data which information of gene, chromosome, enhancer, etc.}

\item{qR}{threshold value for eQTL. Default is 0.1.}

\item{lR}{threshold value for chromatin interaction. Default is 0.1.}

\item{b_use_lNk}{option for incorporating a neighboring relationship for chromatin interaction. Default is FALSE.}

\item{total_iter}{number of iterations. Default is set at 100.}

\item{bDisp}{option for displaying the iteration. Default is FALSE.}
}
\value{
list of the following values:
\describe{
\item{lk_index}{Indice of chromatin interactions.}
\item{qk_index}{Indice of eQTL SNPs.}
\item{qNk_index}{Indice of neighboring eQTLs.}
\item{prob_qk1_all}{Probability of a SNP being eQTL.}
\item{prob_lk1_all}{Probablity of a chromatin interaction.}
\item{prob_qNk_qk_0}{Conditional probability of neighboring eQTL SNPs given that it is not an eQTL SNP.}
\item{prob_qNk_qk_1}{Conditional probability of neighboring eQTL SNPs given that it is an eQTL SNP.}
\item{skq_qk0lk0}{eQTL score given no eQTL SNP (qk=0) and no chromatin interaction (lk=0).}
\item{skq_qk0lk1}{eQTL score given no eQTL SNP (qk=0) and chromatin interaction (lk=1).}
\item{skq_qk1lk0}{eQTL score given eQTL SNP (qk=1) and no chromatin interaction (lk=0).}
\item{skq_qk1lk1}{eQTL score given eQTL SNP (qk=1) and chromatin interaction (lk=1).}
\item{skl_qk0lk0}{chromatin interaction score given no eQTL SNP (qk=0) and no chromatin interaction (lk=0).}
\item{skl_qk0lk1}{chromatin interaction score given no eQTL SNP (qk=0) and chromatin interaction (lk=1).}
\item{skl_qk1lk0}{chromatin interaction score given eQTL SNP (qk=1) and no chromatin interaction (lk=0).}
\item{skl_qk1lk1}{chromatin interaction score given eQTL SNP (qk=1) and chromatin interaction (lk=1).}


}
}
\description{
##Performs Bayesian inference for eQTL analysis
}
\examples{
####################

####Import or upload an input data. Here, an example data (Valid_Input) will be used.
data(List_Example_Data)
Valid_Input = List_Example_Data$Valid_Input
head(Valid_Input)
dim(Valid_Input)
Final_lst = Run_NeQTL2(Valid_Input, qR=0.1, lR=0.1, b_use_lNk=FALSE, total_iter=10, bDisp=FALSE)
#####################


}
