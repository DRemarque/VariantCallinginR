% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FastqQual.R
\name{FastqQual}
\alias{FastqQual}
\title{FastqQual()}
\usage{
FastqQual(
  Forward_Fastq = "No Default",
  Reverse_Fastq = "No Default",
  Paired_End_Data = TRUE,
  Adapter_Check = TRUE,
  Create_Graph = TRUE,
  Output_Base_Name = "FastqQual",
  server = F
)
}
\arguments{
\item{Forward_Fastq}{The single-read or forward paired-end fastq file name (including .fastq file extension)}

\item{Reverse_Fastq}{When Paired_End_Data=TRUE, the reverse paired-end fastq file name (including .fastq file extension). This argument is ignored when Paired_End_Data=FALSE.}

\item{Paired_End_Data}{The input data type, where FALSE corresponds with single-end reads. \link{default=TRUE}}

\item{Adapter_Check}{Whether to check for the presence of adapters when Paired-End_Data=TRUE. \link{default=TRUE}}

\item{Create_Graph}{Whether to create a quality graph when Paired-End_Data=TRUE. \link{default=TRUE}}

\item{Output_Base_Name}{The base name of the generated quality plots. No file extension is required. \link{default="FastqQual"}}

\item{server}{To save the quality plots on a server without gpu access, set this variable to TRUE to use an alternative save method.}
}
\description{
Based on ShortRead, SystempipeR and Rbowtie2, this function allows you to create quality plots and potential adapter sequences of both single and paired-end sequence data.
}
\examples{
To do a general check of paired-end sequence data:
FastqQual(Forward_Fastq="example_fw.fastq", Reverse_Fastq="example_rv.fastq")
To do a general check of single-end sequence data:
FastqQual(Forward_Fastq="example.fastq", Paired_End_Data=FALSE)
To search only for paired-end sequence adapters:
FastqQual(Forward_Fastq="example_fw.fastq", Reverse_Fastq="example_rv.fastq",Create_Graph=FALSE)
To create only a quality graph of paired-end sequence data without searching for adapters:
FastqQual(Forward_Fastq="example_fw.fastq", Reverse_Fastq="example_rv.fastq",Adapter_Check=FALSE)
}
\keyword{fastq}
\keyword{quality}
