% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NeQTL_Main.R
\name{Generate_skl}
\alias{Generate_skl}
\title{##Generate the chromatin interaction score skl.}
\usage{
Generate_skl(df)
}
\arguments{
\item{df}{a data frame that contains
\describe{
\item{distance}{the distance between a SNP and an associated gene.}
\item{Gene_Activity}{an activity score of a gene associated with a SNP.}
\item{Enh_Activity}{the activity score of an enhancer associated with a SNP.}
\item{Correlation}{Correlation between a SNP and a gene.}
}}
}
\value{
a vector which contains skl's and whose length is equal to the row dimension of df.
}
\description{
##Generate the chromatin interaction score skl.
}
\examples{
####################

####Import or upload an input data. Here, an example data (Valid_Input) will be used.
data(List_Example_Data)
Valid_Input = List_Example_Data$Valid_Input
head(Valid_Input)
Vec_skl = Generate_skl(Valid_Input)
summary(Vec_skl)
#####################


}
