#' Bowtie2Align()
#'
#' Based on Rbowtie2 and RSamtools, this function is used to map both single and paired-end sequence data to a reference genome, followed by converting the SAM output into a sorted and indexed BAM file.
#' @param Reference_Genome The reference genome unto which the sequence data will be mapped (including .fasta or .fna file extension).
#' @param Forward_Fastq The single-read or forward paired-end fastq file name (including .fastq file extension)
#' @param Reverse_Fastq When Paired_End_Data=TRUE, the reverse paired-end fastq file name (including .fastq file extension). This argument is ignored when Paired_End_Data=FALSE.
#' @param Paired_End_Data The input data type, where FALSE corresponds with single-end reads. [default=TRUE]
#' @param Index_File_Name The base name for the Bowtie2 reference index. If no index was made before, a new one will be constructed. Otherwise, the index with this basename will be used.
#' @param Output_Base_Name The base name for the alignment files and reports. No file extension is required. [default=Bowtie2_Alignment]
#' @param flag Alignment data is supplied with data information such as read group identifiers (--rg-id), Sample name (--rg SM:), library (rg LB:), Predicted median insert size (--rg PI:) and platform (--rg PL:). [default="--rg-id 1 --rg SM:Sample --rg LB:1 --rg PI:400 --rg PL:ILLUMINA"]
#' @param Threads The number of threads to use in alignment calculations. [default=2]
#' @keywords fastq alignment bowtie2
#' @export
#' @examples
#' Basic Alignment of paired-end example data on a local machine (Low number of threads):
#' Bowtie2Align(Reference_Genome = "example.fasta",Forward_Fastq = "example_fw.fastq",Reverse_Fastq = "example_rv.fastq",Paired_End_Data = T, Index_File_Name = "example", Output_Base_Name = "example", Threads=1)
#' Basic Alignment of single-end example data on a server (High number of threads):
#' Bowtie2Align(Reference_Genome = "example.fasta",Forward_Fastq = "example_fw.fastq",Paired_End_Data = F, Index_File_Name = "example", Output_Base_Name = "example", Threads=24)

Bowtie2Align<-function(Reference_Genome="No Default", Forward_Fastq="No Default", Reverse_Fastq="No Default", Paired_End_Data=TRUE, Index_File_Name="index", Output_Base_Name="Bowtie2_Alignment", Threads=4, flag="--rg-id 1 --rg SM:Sample --rg LB:1 --rg PI:400 --rg PL:ILLUMINA"){
  # Check Input
  if(Reference_Genome=="No Default"){return(message("No Reference genome has been submitted"))
  }else if(file.exists(Forward_Fastq)==FALSE){return(message("No existing Fastq files have been submitted"))
  }else if(file.exists(Reverse_Fastq)==FALSE & Paired_End_Data==TRUE){return(message("No existing reverse fastq file has been submitted, but paired_end_data=TRUE"))
  }else if(file.exists(paste(Output_Base_Name,".bam",sep=""))){return(message("The submitted output files already exist in the working directory. Please choose a different Output_Base_Name."))}

# Setup General Parameter Names
threads<-paste("--threads ",Threads)
threads2<-paste(flag, " --time --threads ",Threads)

# Alignment for Windows
if(.Platform$OS.type=="windows"){

# Build Index if required
if(!file.exists(paste(Index_File_Name,".1.bt2",sep=""))){
    message("Building Reference Genome index")
    message("")
  command<-paste(sep="",'"',.libPaths()[1],'/Rbowtie2/bowtie2-build-s.exe" -c ',file.path(getwd(), Reference_Genome)," ",file.path(getwd(),Index_File_Name)," ",threads)
  system(command, intern=TRUE)
}else{message("Previously made Reference Genome index located");message("")}

# Alignment for single end data
if(Paired_End_Data==FALSE){
  message("Aligning Reads. This will take a while")
  message("")
command<-paste(sep="",'"',.libPaths()[1],'/Rbowtie2/bowtie2-align-s.exe" -x ',file.path(getwd(), Index_File_Name)," -r ",file.path(getwd(),Forward_Fastq)," -S ",file.path(getwd(),Output_Base_Name)," ",threads2)
Log_Table_Alignment<-system(command, intern=TRUE)

}else{
#  Alignment for paired-end data
  message("Aligning Reads. This will take a while")
  message("")
command<-paste(sep="",'"',.libPaths()[1],'/Rbowtie2/bowtie2-align-s.exe" -x ',file.path(getwd(), Index_File_Name)," -1 ",file.path(getwd(),Forward_Fastq)," -2 ",file.path(getwd(),Reverse_Fastq) ," -S ",Output_Base_Name," ",threads2)
Log_Table_Alignment<-system(command, intern=TRUE)
}

}else{
#Alignment for Unix and Mac

# Build Index if required
if(!file.exists(paste(Index_File_Name,".1.bt2",sep=""))){
  file.create(Index_File_Name)
    message("Building Reference Genome index")
    message("")
  index<-Rbowtie2::bowtie2_build(references=Reference_Genome,
                                 bt2Index=file.path(getwd(),Index_File_Name),
                                 overwrite=TRUE,
                                 threads)
  }else{message("Previously made Reference Genome index located");message("")}

# Alignment for single end data
if(Paired_End_Data==FALSE){
  message("Aligning Reads. This will take a while")
  message("")
  Log_Table_Alignment<-Rbowtie2::bowtie2(bt2Index = file.path(getwd(), Index_File_Name),
                                         samOutput = file.path(getwd(), Output_Base_Name),
                                         seq1=Forward_Fastq,
                                         overwrite=TRUE,
                                         threads2)
}else{
# Alignment for paired-end data
  message("Aligning Reads. This will take a while")
  message("")
Log_Table_Alignment<-Rbowtie2::bowtie2(bt2Index = file.path(getwd(), Index_File_Name),
                                       samOutput = file.path(getwd(), Output_Base_Name),
                                       seq1=Forward_Fastq,
                                       seq2=Reverse_Fastq,
                                       overwrite=TRUE,
                                       threads2)
}}
  message("The alignment has finished and has been saved as ", Output_Base_Name)
  message(" ")
# Store alignment rates in a special .txt file
  message("Storing alignment rates")
  message(" ")
Alignment_File_NameX<-paste(gsub(".sam","", Output_Base_Name),"_Alignment_Rates.txt",sep="")
cat(Log_Table_Alignment, file=Alignment_File_NameX,sep="\n")
  message("Aligment rates have been stored as ",Alignment_File_NameX)
  message("")

  message("Converting and processing output...")
# Convert the SAM File to a BAM File
  message("Converting Bam to Sam")
  message("")
BAM_File_Alignment<-Rsamtools::asBam(Output_Base_Name)

# Sort the BAM File
  message("Sorting Bam by coordinate")
  message("")
sortname<-gsub(".bam","", Output_Base_Name)
Rsamtools::sortBam(BAM_File_Alignment,sortname)

# Index the BAM File
  message("Indexing BAM file")
Rsamtools::indexBam(BAM_File_Alignment)

# Final Messages
  message("The entire alignment has been performed succesfully!")
  message("The Alignment Rates Are:")
print(Log_Table_Alignment)
  message(" ")
  message("The output file ",BAM_File_Alignment," and the corresponding .BAM index file can be found at ",getwd())
}
