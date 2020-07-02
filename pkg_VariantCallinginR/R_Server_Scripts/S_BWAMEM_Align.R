#' BWAMEM_Align()
#'
#' This function allows the user to send commands to BWA-MEM, running either natively on linux and Mac or within a docker for Windows, from within R.
#' #' @param Reference_Genome
#' #' @param Forward_Fastq
#' #' @param Reverse_Fastq
#' #' @param Paired_End_Data
#' #' @param Output_Base_Name
#' #' @param flag
#' #' @param pathtobash
#' #' @keywords variant call GATK
#' @export
#' @examples
#'

#
message("This is the BWAMEM_Align() R Server Script")
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
                        help="The input data type, where FALSE corresponds with single-end reads. [default=FALSE]",metavar="logical"),
  optparse::make_option(c("--Output_Base_Name"), type="character", default="BWAMEM_Align_File",
                        help="The base name for the alignment files and reports. No file extension is required. [default=BWAMEM_Align_File]",metavar="character"),
  optparse::make_option(c("--pathtobash"), type="character", default="C:/Program Files/Git/bin/bash.exe",
                        help="Only for running BWA within a docker (Mainly for windows users): the file path to git bash (bash.exe) when installed in a custom file location. Linux users can alter this directory path to make use of GATK within the docker environment as well.",metavar="character"),
  optparse::make_option(c("--pathtoBWA"), type="character", default="bwa",
                        help="Only for running BWA on linux: the file path to the BWA software when a custom installation path has been used. The default is bwa, which assumes BWA is added to your $PATH.",metavar = "character"),
  optparse::make_option(c("--pathtoGATK"), type="character", default="./gatk",
                        help="Only for running GATK on linux: the file path to the GATK software when a custom installation path has been used. The default is ./gatk, which assumes gatk is added to your $PATH.",metavar = "character"),
  optparse::make_option(c("--flag"), type="character", default="--rg-id 1 --rg SM:Sample --rg LB:1 --rg PI:400 --rg PL:ILLUMINA",
                        help="Alignment data is supplied with data information such as read group identifiers (--rg-id), Sample name (--rg SM:), library (rg LB:), Predicted median insert size (--rg PI:) and platform (--rg PL:). [default='--rg-id 1 --rg SM:Sample --rg LB:1 --rg PI:400 --rg PL:ILLUMINA']",metavar="character")
)
opt_parser = optparse::OptionParser(option_list=option_list)
opt= optparse::parse_args(opt_parser)

# Build function
BWAMEM_Align<-function(Reference_Genome=opt$Reference_Genome, Forward_Fastq=opt$Forward_Fastq, Reverse_Fastq=opt$Reverse_Fastq, Paired_End_Data=opt$Paired_End_Data, Output_Base_Name=opt$Output_Base_Name, pathtobash=opt$pathtobash,pathtoBWA=opt$pathtoBWA,pathtoGATK=opt$pathtoGATK ,flag=opt$flag){
  # Check Input
  if(Reference_Genome=="No Default"){return(message("No Reference genome has been submitted"))
  }else if(file.exists(Forward_Fastq)==FALSE){return(message("No existing Fastq files have been submitted"))
  }else if(file.exists(Reverse_Fastq)==FALSE & Paired_End_Data==TRUE){return(message("No existing reverse fastq file has been submitted, but paired_end_data=TRUE"))
  }else if(file.exists(paste(Output_Base_Name,".bam",sep=""))){return(message("The submitted output files already exist in the working directory. Please choose a different Output_Base_Name."))
  }else if(.Platform$OS.type=="windows" & file.exists("C:/Program Files/Git/bin/bash.exe")==FALSE){return(message("Git bash could not be found at the defined location. Change pathtobash to the correct location or install git bash and the docker from https://docs.docker.com/toolbox/toolbox_install_windows/ or https://docs.docker.com/engine/install/#server"))
  }

# Variant Calling Windows version
if(.Platform$OS.type=="windows" | pathtobash!="C:/Program Files/Git/bin/bash.exe"){
  message("Using a docker environment")

# Adapt windows file path to Unix commands
if(.Platform$OS.type=="windows"){
  wd<-paste("//",gsub(":/","/",getwd()),sep="")
  }else{wd<-getwd()}

# Add the file paths in front of input names
Reference_GenomeP<-paste(wd,"/",Reference_Genome,sep="")
Output_Base_NameS<-paste(Output_Base_Name,".sam",sep="")
Output_Base_NameP<-paste(wd,"/",Output_Base_NameS,sep="")
BAM_File_Alignment<-paste(Output_Base_Name,".bam",sep="")

# Start docker
  message("Calling upon docker")
shell(cmd = "'docker-machine restart'", shell=pathtobash, intern=T, flag = "-c")

# Check software installation
  message("Check for BWA and GATK software")
shell(cmd = "'docker pull biocontainers/bwa:v0.7.17-3-deb_cv1
        docker pull broadinstitute/gatk:4.1.3.0'", shell=pathtobash, intern=T, flag = "-c")

# Prepare the reference genome if this has not been done
if(file.exists(paste(gsub("\\..*","",Reference_Genome),".fai",sep=""))==FALSE){
  # Index the fasta file
    message("Indexing reference fasta file")
  command<-paste("'docker run -v ",wd,"://c/Users/Public biocontainers/bwa:v0.7.17-3-deb_cv1 bwa index -a bwtsw -p ",Reference_GenomeP," -b 100000000'",sep="")
  shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")

  if(file.exists(paste(gsub("\\..*","",Reference_Genome),".fai",sep=""))){
    message("The reference index has been made successfully")
  }else{message("No index has been made. Something went wrong, please view the samtools output above for details")}}

# BWA Alignment
if(Paired_End_Data==FALSE){
  message("Aligning Reads. This will take a while")
  message("")
command<-paste("'docker run -v ",wd,"://c/Users/Public biocontainers/bwa:v0.7.17-3-deb_cv1 bwa mem -M -R '@RG\tID:MK\tLB:MK\tPL:ILLUMINA\tSM:MK' -t 4 ", Reference_GenomeP," ", Filtered_File_Name," -o //c/Users/Public/", Output_Base_NameS,"'",sep="")
shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")
}else{
  message("Aligning Reads. This will take a while")
  message("")
command<-paste("'docker run -v ",wd,"://c/Users/Public biocontainers/bwa:v0.7.17-3-deb_cv1 bwa mem -M -R '@RG\tID:MK\tLB:MK\tPL:ILLUMINA\tSM:MK' -t 4 ", Reference_GenomeP, Forward_Fastq," ", Reverse_Fastq," -o //c/Users/Public/", Output_Base_NameS,"'",sep="")
shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")
}

# Sortsam
  message("Sortting SAM by coordinate")
  message("")
command<-paste("'docker run -v ",wd,"://c/Users/Public broadinstitute/gatk:4.1.3.0 ./gatk SortSam -I ",Output_Base_NameP," -O //c/Users/Public/",Output_Base_NameS," -SO coordinate'",sep="")
shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")

# Mark Duplicates
  message("Marking duplicates")
  message("")
command<-paste("'docker run -v ",wd,"://c/Users/Public broadinstitute/gatk:4.1.3.0 ./gatk MarkDuplicates -I ",Output_Base_NameP," -O //c/Users/Public/",BAM_File_Alignment," -SO coordinate'",sep="")
shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")

# Shut down docker
  message("The alignment and sorting steps have finished. Shutting down...")
shell(cmd = "'docker-machine stop'", shell=pathtobash, intern=T, flag = "-c")

}else{
# Variant Calling Unix and Mac version
  message("Using native software")
# Add the file paths in front of input names
wd<-paste("//",gsub(":/","/",getwd()),sep="")
Reference_GenomeP<-paste(wd,"/",Reference_Genome,sep="")
Output_Base_NameS<-paste(Output_Base_Name,".sam",sep="")
Output_Base_NameP<-paste(wd,"/",Output_Base_NameS,sep="")
BAM_File_Alignment<-paste(Output_Base_Name,".bam",sep="")

# Prepare the reference genome if this has not yet been done
if(file.exists(paste(gsub("\\..*","",Reference_Genome),".fai",sep=""))==FALSE){
  # Index the fasta file
    message("Indexing reference fasta file")
  command<-paste(pathtoBWA," index -a bwtsw -p ",Reference_GenomeP," -b 100000000'",sep="")
  system(command, intern=TRUE)

  if(file.exists(paste(gsub("\\..*","",Reference_Genome),".fai",sep=""))){
    message("The reference index has been made successfully")
  }else{message("No index has been made. Something went wrong, please view the samtools output above for details")}}

# BWA Alignment
if(Paired_End_Data==FALSE){
    message("Aligning Reads. This will take a while")
    message("")
  command<-paste(pathtoBWA," mem -M -R '@RG\tID:MK\tLB:MK\tPL:ILLUMINA\tSM:MK' -t 4 ", Reference_GenomeP," ", Forward_Fastq," ", getwd(), "/", Output_Base_NameS,sep="")
  system(command, intern=TRUE)
}else{
    message("Aligning Reads. This will take a while")
    message("")
  command<-paste(pathtoBWA," mem -M -R '@RG\tID:MK\tLB:MK\tPL:ILLUMINA\tSM:MK' -t 4 ", Reference_GenomeP, Forward_Fastq," ", Reverse_Fastq," ", getwd(), "/", Output_Base_NameS,sep="")
  system(command, intern=TRUE)
}

# Sortsam
  message("Sorting SAM by coordinate")
  message("")
command<-paste(pathtoGATK," SortSam -I ",Output_Base_NameP," ", getwd(), "/", Output_Base_NameS," -SO coordinate",sep="")
system(command, intern=TRUE)

# Mark Duplicates
  message("Marking duplicates")
  message("")
command<-paste(pathtoGATK," MarkDuplicates -I ",Output_Base_NameP," ", getwd(), "/", BAM_File_Alignment," -SO coordinate",sep="")
system(command, intern=TRUE)
 message("The alignment has finished. The files can be found in ",getwd())
}
}
# Execute alignment
BWAMEM_Align()
