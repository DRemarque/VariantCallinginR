% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AlignQual.R
\name{AlignQual}
\alias{AlignQual}
\title{AlignQual()}
\usage{
AlignQual(
  Alignment_File_Name = "No Default",
  Reference_Length = "No Default",
  Read_Length = "No Default",
  server = F
)
}
\arguments{
\item{Alignment_File_Name}{The file name of the alignment to be assessed (including .bam file extension)}

\item{Reference_Length}{The length of the utilised reference genome in basepairs. This is used to calculate the coverage based on the Lander/Waterman equation (Coverage=Length*Number of mapped reads/ Genome Length)}

\item{Read_Length}{The general length of the sequenced reads in basepairs. This generally between 100 and 150 bp for Illumina.}

\item{server}{To save the quality plots on a server without gpu access, set this variable to TRUE to use an alternative save method.}
}
\description{
Based on RSamtools and ggplot2, this function is used collect quality data from an alignment BAM file, which is then reported in a mapQ graph, a table read numbers (un)mapped  and ,optionally, coverage and alignment rates are outputted as well in the console. .
}
\examples{
# Plot the quality of a file
AlignQual("myalignment.bam")
# Calculate coverage as well
AlignQual("myalignment.bam",Reference_Length=1000,Read_Length=150)
}
\keyword{alignment}
\keyword{bowtie2}
\keyword{fastq}
