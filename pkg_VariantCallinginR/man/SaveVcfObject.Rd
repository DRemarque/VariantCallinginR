% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SaveVcfObject.R
\name{SaveVcfObject}
\alias{SaveVcfObject}
\title{SaveVcfObject()}
\usage{
SaveVcfObject(vcf_list = "No Default", Output_Base_Name = "VcfObject")
}
\arguments{
\item{vcf_list}{The R object into which the vcf data was stored by readVcf.}

\item{Output_Base_Name}{The desired output base name for the created excel file. \link{default="VcfObject"}}
}
\description{
This function allows the user to save a vcf R object created by readVcf() into an excel file where sample specific parameters are stored in separate sheets. Note that this saving method is only functional for small datasets and not whole genome data.
}
\examples{
# Simple Save
SaveVcfObject(vcf_list=vcf_object,Output_Base_Name="My_Variants")
}
\keyword{excel}
\keyword{save}
\keyword{vcf}
\keyword{xlsx}
