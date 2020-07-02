#
message("This is the MergeAndGenotype() R Server Script")
message("Starting up")
# Index Input
option_list = list(
  optparse::make_option(c("--Reference_Genome"), type="character", default="No Default",
                        help="The reference genome unto which the sequence data has been mapped (including .fasta or .fna file extension).",metavar = "character"),
  optparse::make_option(c("--GVCF_Files"), type="logical", default="No Default",
                        help="The g.vcf file names to merge into one. Input the names into a character string like so: c('one.g.vcf','two.g.vcf') or as a single string 'one.g.vcf two.g.vcf'",metavar = "logical"),
  optparse::make_option(c("--Output_Base_Name"), type="logical", default="Merged_GVCF",
                        help="The desired output base name for the variant call format (VCF) file. [default='Gatk_Call'']",metavar = "logical"),
  optparse::make_option(c("--pathtoGATK"), type="character", default="./gatk",
                        help="Only for running GATK on linux: the file path to the GATK software when a custom installation path has been used. The default is ./gatk, which assumes gatk is added to your $PATH.",metavar = "character"),
  optparse::make_option(c("--pathtobash"), type="character", default="C:/Program Files/Git/bin/bash.exe",
                        help="Only for running GATK within a docker (Mainly for windows users): the file path to git bash (bash.exe) when installed in a custom file location. Linux users can alter this directory path to make use of GATK within the docker environment as well.",metavar = "character")
)
opt_parser = optparse::OptionParser(option_list=option_list)
opt= optparse::parse_args(opt_parser)

# Build function
MergeAndGenotype<-function(Reference_Genome=opt$Reference_Genome, GVCF_Files=opt$GVCF_Files, Output_Base_Name=opt$Output_Base_Name,pathtoGATK=opt$pathtoGATK, pathtobash=opt$pathtobash){
# Check Input
if(file.exists(Reference_Genome)==FALSE){return(message("No existing Reference genome has been submitted"))
}else if(length(GVCF_Files)<=1&grepl(" ",GVCF_Files)[1]==T){return(message("Submit at least two g.vcf files to merge."))}

# Index input
if(length(GVCF_Files)==1 & grepl(" ",GVCF_Files)[1]==T){
GVCF_Files<-unlist(strsplit(GVCF_Files," "))}
if(length(GVCF_Files)==1 & grepl(" , ",GVCF_Files)[1]==T){
  GVCF_Files<-unlist(strsplit(GVCF_Files," , "))}

GVCF_Files<-paste(paste("--variant ",GVCF_Files,sep=""),collapse=" ")
Output_Base_Name1<-paste(Output_Base_Name,"_temp.g.vcf",sep = "")
Output_Base_Name2<-paste(Output_Base_Name,".vcf",sep = "")

## Merge on Windows
if(.Platform$OS.type=="windows" | pathtobash!="C:/Program Files/Git/bin/bash.exe"){
  message("Using docker environment")

# Adapt windows file path to Unix commands
if(.Platform$OS.type=="windows"){
  wd<-paste("//",gsub(":/","/",getwd()),sep="")
}else{wd<-getwd()}

# Add the quotes around bash
pathtobash<-paste("'",pathtobash,"'",sep="")

# Start docker
message("Calling upon docker")
shell(cmd = "'docker-machine restart'", shell=pathtobash, intern=T, flag = "-c")

# Check software installation
message("Check for GATK and Samtools software")
shell(cmd = "'docker pull broadinstitute/gatk:4.1.3.0'", shell=pathtobash, intern=T, flag = "-c")
# Combine g.vcfs
  message("Combining g.vcf files...")
  message("")
command<-paste("'docker run -v ",wd,"://c/Users/Public broadinstitute/gatk:4.1.3.0 ./gatk CombineGVCFs -R ",Reference_Genome," ",GVCF_Files," -O ",Output_Base_Name1,"'",sep="")
shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")

# Genotype combined file
  message("Genotyping results and converting to vcf...")
  message("")
command<-paste(pathtoGATK," GenotypeGVCFs -R ",Reference_Genome," -V ",Output_Base_Name1," -O ",Output_Base_Name2)
shell(cmd = command, shell=pathtobash, intern=T, flag = "-c")
  message("Merging has been completed. The final file as been saved as ",Output_Base_Name2," in ",getwd())

# Shut down docker
  message("Shutting down...")
shell(cmd = "'docker-machine stop'", shell=pathtobash, intern=T, flag = "-c")

}else{
## Merge on linux and Mac
# Combine g.vcfs
  message("Combining g.vcf files...")
  message("")
command<-paste(pathtoGATK," CombineGVCFs -R ",Reference_Genome," ",GVCF_Files," -O ",Output_Base_Name1,sep="")
system(command, intern=TRUE)
# Genotype combined file
  message("Genotyping results and converting to vcf...")
  message("")
command<-paste(pathtoGATK," GenotypeGVCFs -R ",Reference_Genome," -V ",Output_Base_Name1," -O ",Output_Base_Name2)
  message("Merging has been completed. The final file as been saved as ",Output_Base_Name2," in ",getwd())
}
}

# Exectute merging
MergeAndGenotype()
