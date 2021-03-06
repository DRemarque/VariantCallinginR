% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VariantPlot.R
\name{VariantPlot}
\alias{VariantPlot}
\title{VariantPlot()}
\usage{
VariantPlot(
  vcf_list = "No Default",
  plot = "No Default",
  save_name = "No Default",
  server = F
)
}
\arguments{
\item{vcf_list}{R object into which vcf data was stored with ReadVcf()}

\item{plot}{If only a specific parameter has to be plotted, insert the name with plot='parameter name'. This only works for the parameters plotted when plot='No Default'. Other parameters give an empty plot.}

\item{save_name}{The generated quality plot can be be saved as a pdf file by setting save to the desired base save name (without .pdf extension)}

\item{server}{To save the quality plots on a server without gpu access, set this variable to TRUE to use an alternative save method.}
}
\description{
Based on ggplot2, this function is used to plot relevant quality parameters of variants stored in an R object with ReadVcf()
}
\examples{
# Plot only the quality to depth ratio
VariantPlot(vcf_list=vcf_object,plot="QD")
# Plot all graphs and save the plot on a server
VariantPlot(vcf_list=vcf_object,save_name="Graph",server=T)
}
\keyword{filter}
\keyword{vcf}
\keyword{vcfR}
