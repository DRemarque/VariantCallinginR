#' AlignQual()
#'
#' Based on RSamtools and ggplot2, this function is used collect quality data from an alignment BAM file, which is then reported in a mapQ graph, a table read numbers (un)mapped  and ,optionally, coverage and alignment rates are outputted as well in the console. .
#' @param Alignment_File_Name The file name of the alignment to be assessed (including .bam file extension)
#' @param Reference_Length The length of the utilised reference genome in basepairs. This is used to calculate the coverage based on the Lander/Waterman equation (Coverage=Length*Number of mapped reads/ Genome Length)
#' @param Read_Length The general length of the sequenced reads in basepairs. This generally between 100 and 150 bp for Illumina.
#' @param server To save the quality plots on a server without gpu access, set this variable to TRUE to use an alternative save method.
#' @keywords fastq alignment bowtie2
#' @export
#' @examples
#' # Plot the quality of a file
#' AlignQual("myalignment.bam")
#' # Calculate coverage as well
#' AlignQual("myalignment.bam",Reference_Length=1000,Read_Length=150)

AlignQual<-function(Alignment_File_Name="No Default",Reference_Length="No Default",Read_Length="No Default",server=F){

# Check Input
if(Alignment_File_Name=="No Default"){return(message("No alignment file has been submitted"))
}else if(Reference_Length=="No Default" & Read_Length=="No Default"){
    message("Coverage and alignment rates will not be calculated due to missing Reference_Length")}

# Fetch Mapping Qualities
message("Collecting data...")
param<-Rsamtools::ScanBamParam(what=c("mapq"))
BAM_Frame<-as.data.frame(Rsamtools::scanBam(file=Alignment_File_Name ,param=param))

if(server==F){
# Plot Mapping Qualities
  message(" ")
  message("Plotting quality...")
Bam_Qual_Plot<-paste(gsub(".bam","", Alignment_File_Name),"_Align_Qual.pdf",sep="")
pdf(file=file.path(getwd(),Bam_Qual_Plot), width=11.7, height=5.7)
ggplot2::ggplot(data=BAM_Frame)+ggplot2::geom_density(ggplot2::aes(x=mapq),stat="bin")
dev.off()
  message(" ")
  message("The Mapping Quality Graph has been saved in ",getwd()," as ",Bam_Qual_Plot)
}else if(server==T){
  message(" ")
  message("Plotting quality...")
Bam_Qual_Plot<-paste(gsub(".bam","", Alignment_File_Name),"_Align_Qual.pdf",sep="")
bitmap(file=file.path(getwd(),Bam_Qual_Plot), width=11.7, height=5.7, type="pdfwrite", res=100)
ggplot2::ggplot(data=BAM_Frame)+ggplot2::geom_density(ggplot2::aes(x=mapq),stat="bin")
dev.off()
  message(" ")
  message("The Mapping Quality Graph has been saved in ",getwd()," as ",Bam_Qual_Plot)
  }

# View Mapped/Unmapped ratios
  message(" ")
  message("Collecting BAM Statistics")
BAM_Summary<-capture.output(Rsamtools::quickBamFlagSummary(Alignment_File_Name))
Statistics<-paste(gsub(".bam","", Alignment_File_Name),"_Statistics.txt",sep="")
cat(BAM_Summary, file=Statistics,sep="\n")
  message(" ")
  message("The General_Alignment_Statistics have been saved as ", Statistics," in ", getwd())

# Calculate coverage and alignment rates and if possible
if(Reference_Length!="No Default" & Read_Length!="No Default"){
  stats<-read.delim(Statistics,stringsAsFactor=F)
  # For paired end data
  if(is.na(stats[25])){
    Coverage<-round((
      as.numeric(sub("\\D*(\\d+).*", "\\1", stats[3,]))-
        as.numeric(sub("\\D*(\\d+)\\D*(\\d+).*", "\\2", stats[20,]))-
        as.numeric(sub("\\D*(\\d+)\\D*(\\d+).*", "\\2", stats[25,])))*Read_Length
      /Reference_Length)
    Coverage<-paste(Coverage,"x",sep="")
    Rate<-round((as.numeric(sub("\\D*(\\d+).*", "\\1", stats[3,]))-
                   as.numeric(sub("\\D*(\\d+)\\D*(\\d+).*", "\\2", stats[20,]))-
                   as.numeric(sub("\\D*(\\d+)\\D*(\\d+).*", "\\2", stats[25,])))/
                  as.numeric(sub("\\D*(\\d+).*", "\\1", stats[3,]))*100,2)
    Rate<-paste(Rate,"%",sep="")
  }else{
  # For single end data
    Coverage<-round(
      (as.numeric(sub("\\D*(\\d+).*", "\\1", stats[3,]))-
        as.numeric(sub("\\D*(\\d+)\\D*(\\d+).*", "\\2", stats[20,]))-
        as.numeric(sub("\\D*(\\d+)\\D*(\\d+).*", "\\2", stats[25,])))*Read_Length/Reference_Length)
    Coverage<-paste(Coverage,"x",sep="")
    Rate<-round((as.numeric(sub("\\D*(\\d+).*", "\\1", stats[3,]))-
                 as.numeric(sub("\\D*(\\d+)\\D*(\\d+).*", "\\2", stats[19,]))
                 /as.numeric(sub("\\D*(\\d+).*", "\\1", stats[3,]))*100),2)
    Rate<-paste(Rate,"%",sep="")
  }
  message("The coverage has been estimated to be ",Coverage, " with an alignment rate of ",Rate)
  message("")
}
}
