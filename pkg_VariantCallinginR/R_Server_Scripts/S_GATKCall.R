#
message("This is the GATKCall() R Server Script")
message("Starting up")
# Index Input
option_list = list(
  optparse::make_option(c("--Reference_Genome"), type="character", default="No Default",
                        help="The reference genome unto which the sequence data has been mapped (including .fasta or .fna file extension).",metavar = "character"),
  optparse::make_option(c("--Alignment_File_Name"), type="character", default="No Default",
                        help="The file name of the alignment data from which variants are to be called (including .bam )",metavar = "character"),
  optparse::make_option(c("--Output_Base_Name"), type="character", default="Gatk_Call",
                        help="The desired output base name for the variant call format (VCF) file. [default='Gatk_Call]'",metavar = "character"),
  optparse::make_option(c("--gvcf.mode"), type="logical", default=FALSE,
                        help="Whether to store the variant call data in a genomic file (g.vcf). Use this when calling variants from multiple sample alignments that are to be merged later into one large vcf. [default=FALSE]",metavar="logical"),
  optparse::make_option(c("--BAM_Validation"), type="logical", default=FALSE,
                        help="Whether to check for anomalies in the alignment file before performing variant calling. [default=FALSE]",metavar = "logical"),
  optparse::make_option(c("--BAM_Repair"), type="logical", default=FALSE,
                        help="Whether to repair or replace the alignment file header. Enable this if the .BAM file is malformed or if the sample name is incorrect by inserting the new header as a character. An example of such a header is '-SM Example -LB 1 -PL Illumina -PU 10578252' [default=FALSE] ",metavar = "logical"),
  optparse::make_option(c("--pathtoGATK"), type="character", default="./samtools",
                        help="Only for running GATK on linux when BamValidation=T or BamRepair=T: the file path to the samtools software when a custom installation path has been used. The default is ./samtools, which assumes samtools is added to your $PATH.",metavar = "character"),
  optparse::make_option(c("--pathtoGATK"), type="character", default="./gatk",
                        help="Only for running GATK on linux: the file path to the GATK software when a custom installation path has been used. The default is ./gatk, which assumes gatk is added to your $PATH.",metavar = "character"),
  optparse::make_option(c("--pathtobash"), type="character", default="C:/Program Files/Git/bin/bash.exe",
                        help="Only for running GATK within a docker (Mainly for windows users): the file path to git bash (bash.exe) when installed in a custom file location. Linux users can alter this directory path to make use of GATK within the docker environment as well.",metavar = "character")
)
opt_parser = optparse::OptionParser(option_list=option_list)
opt= optparse::parse_args(opt_parser)

# Build function
GATKCall<-function(Reference_Genome=opt$Reference_Genome, Alignment_File_Name=opt$Alignment_File_Name, Output_Base_Name=opt$Output_Base_Name, gvcf.mode=opt$gvcf.mode, BAM_Validation=opt$BAM_Validation, BAM_Repair=opt$BAM_Repair, pathtobash=opt$pathtobash, pathtoGATK=opt$pathtoGATK,pathtoSamtools=opt$pathtoSamtools){
  # Check Input
  if(file.exists(Reference_Genome)==FALSE){return(message("No existing Reference genome has been submitted"))
  }else if(file.exists(Alignment_File_Name)==FALSE){return(message("No existing alignment file has been submitted"))
  }else if(BAM_Repair==TRUE){return(message("BAM file repair has been turned on, but no replacement header has been supplied. Please do so by inserting the new header as a character behind BAM_Repair= "))
  }else if(.Platform$OS.type=="windows" & file.exists("C:/Program Files/Git/bin/bash.exe")==FALSE){return(message("Git bash could not be found at the defined location. Change pathtobash to the correct location or install git bash and the docker from https://docs.docker.com/toolbox/toolbox_install_windows/ or https://docs.docker.com/engine/install/#server"))}

  # Variant Calling Windows version
  if(.Platform$OS.type=="windows" | pathtobash!="C:/Program Files/Git/bin/bash.exe"){
    message("Using docker environment")

    # Adapt windows file path to Unix commands
    if(.Platform$OS.type=="windows"){
      wd<-paste("//",gsub(":/","/",getwd()),sep="")
    }else{wd<-getwd()}

    # Add the file paths in front of input names
    Reference_GenomeP<-paste(wd,"/",Reference_Genome,sep="")
    Alignment_File_NameP<-paste(wd,"/",Alignment_File_Name,sep="")
    pathtobash<-paste("'",pathtobash,"'",sep="")

    # Start docker
    message("Calling upon docker")
    shell(cmd = "'docker-machine restart'", shell=pathtobash, intern=T, flag = "-c")

    # Check software installation
    message("Check for GATK and Samtools software")
    shell(cmd = "'docker pull broadinstitute/gatk:4.1.3.0
          docker pull biocontainers/samtools:v1.9-4-deb_cv1'", shell=pathtobash, intern=T, flag = "-c")

    # Prepare the reference genome if this has not yet been done
    if(file.exists(paste(gsub("\\..*","",Reference_Genome),".dict",sep=""))==FALSE |
       file.exists(paste(gsub("\\..*","",Reference_Genome),".fai",sep=""))==FALSE){
      # Index the fasta file
      message("No previous index found. Indexing reference fasta file")
      command<-paste("'docker run -v ",wd,"://c/Users/Public biocontainers/samtools:v1.9-4-deb_cv1 samtools faidx ",Reference_GenomeP,"'",sep="")
      shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")

      if(file.exists(paste(gsub("\\..*","",Reference_Genome),".fai",sep=""))){
        message("The reference index has been made successfully")
      }else{message("Something went wrong and no index has been made. Please view the samtools output above for details")}

      # Create Fasta Dictionary
      message("Creating reference fasta dictrionary")
      command<-paste("'docker run -v ",wd,"://c/Users/Public broadinstitute/gatk:4.1.3.0 ./gatk CreateSequenceDictionary -R ",Reference_GenomeP," -O //c/Users/Public/",gsub("\\..*","",Reference_Genome),"'",sep="")
      shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")

      if(file.exists(paste(gsub("\\..*","",Reference_Genome),".dict",sep=""))){
        message("The reference dictionary has been made successfully")
      }else{message("No dictionary has been made. Something went wrong, please view the GATK output above for details")}}

    if(BAM_Validation==T){
      # Validate Bam File
      message("Validating Alignment File")
      command<-paste("'docker run -v ",wd,"://c/Users/Public broadinstitute/gatk:4.1.3.0 ./gatk ValidateSamFile -I ",Alignment_File_NameP," -MODE SUMMARY'",sep="")
      shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")}

    if(BAM_Repair!=F){
      message("Repairing Bam File")
      command<-paste("'docker run -v ",wd,"://c/Users/Public broadinstitute/gatk:4.1.3.0 ./gatk AddOrReplaceReadGroups -I ",Alignment_File_NameP," -O //c/Users/Public/",Alignment_File_Name," ",BAM_Repair,"'",sep="")
      shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")}

    # Call Variants
    if(gvcf.mode==FALSE){
      message("Building BAM index")
      command<-paste("'docker run -v ",wd,"://c/Users/Public broadinstitute/gatk:4.1.3.0 ./gatk BuildBamIndex -I ",Alignment_File_NameP,"'",sep="")
      shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")
      message(" ")
      message("Performing single sample variant calling")
      Output_Base_Name<-paste(Output_Base_Name,".vcf",sep="")
      command<-paste("'docker run -v ",wd,"://c/Users/Public broadinstitute/gatk:4.1.3.0 ./gatk HaplotypeCaller -R ",Reference_GenomeP," -I ",Alignment_File_NameP," -O //c/Users/Public/",Output_Base_Name,"'",sep="")
      shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")
      message("The variant call has been saved in ",getwd()," as ",Output_Base_Name)
    }else{
      message("Building BAM index")
      command<-paste("'docker run -v ",wd,"://c/Users/Public broadinstitute/gatk:4.1.3.0 ./gatk BuildBamIndex -I ",Alignment_File_NameP,"'",sep="")
      shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")
      message(" ")
      message("Performing gVCF variant calling")
      Output_Base_Name<-paste(Output_Base_Name,".g.vcf",sep="")
      command<-paste("'docker run -v ",wd,"://c/Users/Public broadinstitute/gatk:4.1.3.0 ./gatk HaplotypeCaller -R ",Reference_GenomeP," -I ",Alignment_File_NameP," -O //c/Users/Public/",Output_Base_Name," -ERC GVCF'",sep="")
      shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")
      message("The variant call has been saved in ",getwd()," as ",Output_Base_Name)
    }

    # Shut down docker
    message("Variant calling has been performed. Shutting down...")
    shell(cmd = "'docker-machine stop'", shell=pathtobash, intern=T, flag = "-c")

  }else{
    # Variant Calling Unix and Mac version
    message("Using native software")
    # Add the file paths in front of input names
    wd<-paste("//",gsub(":/","/",getwd()),sep="")
    Reference_GenomeP<-paste(wd,"/",Reference_Genome,sep="")
    Alignment_File_NameP<-paste(wd,"/",Alignment_File_Name,sep="")

    # Prepare the reference genome if this has not yet been done
    if(file.exists(paste(gsub("\\..*","",Reference_Genome),".dict",sep=""))==FALSE |
       file.exists(paste(gsub("\\..*","",Reference_Genome),".fai",sep=""))==FALSE){
      # Index the fasta file
      message("No previous index found. Indexing reference fasta file")
      command<-paste(pathtoSamtools," faidx ",Reference_GenomeP,sep="")
      system(command, intern=TRUE)

      # Create Fasta Dictionary
      message("Creating reference fasta dictrionary")
      command<-paste(pathtoGATK," CreateSequenceDictionary -R ",Reference_GenomeP," ", getwd(), "/",gsub("\\..*","",Reference_Genome),sep="")
      system(command, intern=TRUE)

      if(file.exists(paste(gsub("\\..*","",Reference_Genome),".fai",sep=""))){
        message("The reference index has been made successfully")
      }else{message("No index has been made. Something went wrong, please view the samtools output above for details")
      }
    }

    if(BAM_Validation==T){
      # Validate Bam File
      message("Validating Alignment File")
      command<-paste(pathtoGATK," ValidateSamFile -I ",Alignment_File_NameP," -MODE SUMMARY",sep="")
      system(command, intern=TRUE)

      if(file.exists(paste(gsub("\\..*","",Reference_Genome),".dict",sep=""))){
        message("The reference dictionary has been made successfully")
      }else{message("No dictionary has been made. Something went wrong, please view the GATK output above for details")}
    }

    if(BAM_Repair!=FALSE){
      message("Repairing Bam File")
      command<-paste(pathtoGATK," AddOrReplaceReadGroups -I ",Alignment_File_NameP," ", getwd(), "/",Alignment_File_Name," ",BAM_Repair,sep="")
      system(command, intern=TRUE)
    }

    # Call Variants
    if(gvcf.mode==FALSE){
      message("Buidling BAM Index")
      command<-paste(pathtoGATK," BuildBamIndex -I ",Alignment_File_NameP,sep="")
      system(command, intern=TRUE)
      message(" ")
      message("Performing single sample variant calling")
      Output_Base_Name<-paste(Output_Base_Name,".vcf",sep="")
      command<-paste(pathtoGATK," HaplotypeCaller -R ",Reference_GenomeP," -I ",Alignment_File_NameP," ", getwd(), "/",Output_Base_Name,sep="")
      system(command, intern=TRUE)
      message("The variant call has been saved in ",getwd()," as ",Output_Base_Name)
    }else{
      message("Buidling BAM Index")
      command<-paste(pathtoGATK," BuildBamIndex -I ",Alignment_File_NameP,sep="")
      system(command, intern=TRUE)
      message(" ")
      message("Performing gVCF variant calling")
      Output_Base_Name<-paste(Output_Base_Name,".g.vcf",sep="")
      command<-paste(pathtoGATK," HaplotypeCaller -R ",Reference_GenomeP," -I ",Alignment_File_NameP," ", getwd(), "/",Output_Base_Name," -ERC GVCF'",sep="")
      system(command, intern=TRUE)
      message("The variant call has been saved in ",getwd()," as ",Output_Base_Name)
    }
  }
}

# Execute variant call
GATKCall()
