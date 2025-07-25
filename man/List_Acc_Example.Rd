% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{List_Acc_Example}
\alias{List_Acc_Example}
\title{Result of accuracy measures obtained from NeQTL and Naive Bayesian}
\format{
A list of 100 (= number of percentiles) pairs of 4 accuracy measures (true positive, false positive, false negative, true negative) after adjusting the distance obtained by NeQTL and naive Bayesian:
\describe{
  \item{Accuracy_Mat1}{100x4 matrix of 4 accuracy measures of NeQTL over 100 iterations using the adjusted distance}
  \item{Accuracy_Mat2}{100x4 matrix of 4 accuracy measures of naive Bayesian over 100 iterations using the adjusted distance}
  
}
}
\usage{
data(List_Acc_Example)
}
\description{
#####################
}
\examples{

data(List_Acc_Example)
Acc_NeQTL = List_Acc_Example$Accuracy_Mat1
Acc_Naive = List_Acc_Example$Accuracy_Mat2

Draw_Adjusted_ROC(Acc_NeQTL, Acc_Naive, "NeQTL", "Naive")


}
\seealso{
\code{Compute_Accuracy}, \code{Draw_Adjusted_ROC}
}
\keyword{datasets}
