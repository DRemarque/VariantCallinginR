% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MergeAndGenotype.R
\name{MergeAndGenotype}
\alias{MergeAndGenotype}
\title{MergeAndGenotype()}
\usage{
MergeAndGenotype(
  Reference_Genome = "No Default",
  GVCF_Files = "No Default",
  Output_Base_Name = "Merged_GVCF",
  pathtoGATK = "./gatk",
  pathtobash = "C:/Program Files/Git/bin/bash.exe"
)
}
\arguments{
\item{Reference_Genome}{The reference genome unto which the sequence data has been mapped (including .fasta or .fna file extension).}

\item{GVCF_Files}{The g.vcf file names to merge into one. Input the names into a character string like so: c('one.g.vcf','two.g.vcf') or in a single string like so: "one.g.vcf two.g.vcf"}

\item{Output_Base_Name}{The desired output base name for the variant call format (VCF) file. \link{default="Gatk_Call"}}

\item{pathtoGATK}{Only for running GATK on linux: the file path to the GATK software when a custom installation path has been used. The default is ./gatk, which assumes gatk is added to your $PATH.}

\item{pathtobash}{Only for running GATK within a docker (Mainly for windows users): the file path to git bash (bash.exe) when installed in a custom file location. Linux users can alter this directory path to make use of GATK within the docker environment as well.}
}
\description{
This function allows the user to send commands to the GATK, running either natively on linux and Mac or within a docker for Windows, from within R. This call merges g.vcf files into a final .vcf file.
}
\examples{
# Basic example
MergeAndGenotype(Reference_Genome="ref.fasta",GVCF_Files=c("one.g.vcf","two.g.vcf"),Output_Base_Name="Final_Variants")
}
\keyword{GATK}
\keyword{call}
\keyword{variant}
