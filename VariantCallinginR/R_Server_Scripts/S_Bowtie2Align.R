#
message("This is the Bowtie2Align() R Server Script")
message("Starting up")
# Index Input
option_list = list(
  optparse::make_option(c("--Reference_Genome"), type="character", default="No Default",
                        help="The reference genome unto which the sequence data will be mapped (including .fasta or .fna file extension).",metavar = "character"),
  optparse::make_option(c("--Forward_Fastq"), type="character", default="No Default",
                        help="The single-read or forward paired-end fastq file name (including .fastq file extension)",metavar = "character"),
  optparse::make_option(c("--Reverse_Fastq"), type="character", default="No Default",
                        help="When Paired_End_Data=TRUE, the reverse paired-end fastq file name (including .fastq file extension). This argument is ignored when Paired_End_Data=FALSE.",metavar = "character"),
  optparse::make_option(c("--Paired_End_Data"), type="logical", default=FALSE,
                        help="The input data type, where FALSE corresponds with single-end reads. [default=FALSE]",metavar = "logical"),
  optparse::make_option(c("--Index_File_Name"), type="character", default="Index",
                        help="The base name for the Bowtie2 reference index. If no index was made before, a new one will be constructed. Otherwise, the index with this basename will be used.",metavar = "character"),
  optparse::make_option(c("--Output_Base_Name"), type="character", default="Bowtie2Align_file",
                        help="The base name for the alignment files and reports. No file extension is required. [default=Bowtie2Align_file]",metavar = "character"),
  optparse::make_option(c("--flag"), type="character", default="--rg-id 1 --rg SM:Sample --rg LB:1 --rg PI:400 --rg PL:ILLUMINA",
                        help="Alignment data is supplied with data information such as read group identifiers (--rg-id), Sample name (--rg SM:), library (rg LB:), Predicted median insert size (--rg PI:) and platform (--rg PL:). [default='--rg-id 1 --rg SM:Sample --rg LB:1 --rg PI:400 --rg PL:ILLUMINA']",metavar = "character"),
  optparse::make_option(c("--Threads"), type="integer", default=1,
                        help="The number of threads to use in alignment calculations. [default=1]",metavar="integer")
)
opt_parser = optparse::OptionParser(option_list=option_list)
opt= optparse::parse_args(opt_parser)

Bowtie2Align<-function(Reference_Genome=opt$Reference_Genome, Forward_Fastq=opt$Forward_Fastq, Reverse_Fastq=opt$Reverse_Fastq, Paired_End_Data=opt$Paired_End_Data, Index_File_Name=opt$Index_File_Name, Output_Base_Name=opt$Output_Base_Name, Threads=opt$Threads, flag=opt$flag){
  # Check Input
  if(Reference_Genome=="No Default"){return(message("No Reference genome has been submitted"))
  }else if(file.exists(Forward_Fastq)==FALSE){return(message("No existing Fastq files have been submitted"))
  }else if(file.exists(Reverse_Fastq)==FALSE & Paired_End_Data==TRUE){return(message("No existing reverse fastq file has been submitted, but paired_end_data=TRUE"))
  }else if(file.exists(paste(Alignment_File_Name,".bam",sep=""))){return(message("The submitted output files already exist in the working directory. Please choose a different Output_Base_Name."))}

# Setup General Parameter Names
threads<-paste("--threads ",Threads)
threads2<-paste(flags, " --time --threads ",Threads)

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
command<-paste(sep="",'"',.libPaths()[1],'/Rbowtie2/bowtie2-align-s.exe" -x ',file.path(getwd(), Index_File_Name)," -r ",file.path(getwd(),Forward_Fastq)," -S ",file.path(getwd(),Alignment_File_Name)," ",threads2)
Log_Table_Alignment<-system(command, intern=TRUE)

}else{
#  Alignment for paired-end data
  message("Aligning Reads. This will take a while")
  message("")
command<-paste(sep="",'"',.libPaths()[1],'/Rbowtie2/bowtie2-align-s.exe" -x ',file.path(getwd(), Index_File_Name)," -1 ",file.path(getwd(),Forward_Fastq)," -2 ",file.path(getwd(),Reverse_Fastq) ," -S ",Alignment_File_Name," ",threads2)
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
                                         samOutput = file.path(getwd(), Alignment_File_Name),
                                         seq1=Forward_Fastq,
                                         overwrite=TRUE,
                                         threads2)
}else{
# Alignment for paired-end data
  message("Aligning Reads. This will take a while")
  message("")
Log_Table_Alignment<-Rbowtie2::bowtie2(bt2Index = file.path(getwd(), Index_File_Name),
                                       samOutput = file.path(getwd(), Alignment_File_Name),
                                       seq1=Forward_Fastq,
                                       seq2=Reverse_Fastq,
                                       overwrite=TRUE,
                                       threads2)
}}
  message("The alignment has finished and has been saved as ", Alignment_File_Name)
  message(" ")
# Store alignment rates in a special .txt file
  message("Storing alignment rates")
  message(" ")
Alignment_File_NameX<-paste(gsub(".sam","", Alignment_File_Name),"_Alignment_Rates.txt",sep="")
cat(Log_Table_Alignment, file=Alignment_File_NameX,sep="\n")
  message("Aligment rates have been stored as ",Alignment_File_NameX)
  message("")

  message("Converting and processing output...")
# Convert the SAM File to a BAM File
  message("Converting Bam to Sam")
  message("")
BAM_File_Alignment<-Rsamtools::asBam(Alignment_File_Name)

# Sort the BAM File
  message("Sorting Bam by coordinate")
  message("")
sortname<-gsub(".bam","", BAM_File_Alignment)
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

Bowtie2Align()
