#' MergeAndGenotype()
#'
#' This function allows the user to send commands to the GATK, running either natively on linux and Mac or within a docker for Windows, from within R. This call merges g.vcf files into a final .vcf file.
#' @param Reference_Genome The reference genome unto which the sequence data has been mapped (including .fasta or .fna file extension).
#' @param GVCF_Files The g.vcf file names to merge into one. Input the names into a character string like so: c('one.g.vcf','two.g.vcf') or in a single string like so: "one.g.vcf two.g.vcf"
#' @param Output_Base_Name The desired output base name for the variant call format (VCF) file. [default="Gatk_Call"]
#' @param pathtoGATK Only for running GATK on linux: the file path to the GATK software when a custom installation path has been used. The default is ./gatk, which assumes gatk is added to your $PATH.
#' @param pathtobash Only for running GATK within a docker (Mainly for windows users): the file path to git bash (bash.exe) when installed in a custom file location. Linux users can alter this directory path to make use of GATK within the docker environment as well.
#' @keywords variant call GATK
#' @export
#' @examples
#' # Basic example
#' MergeAndGenotype(Reference_Genome="ref.fasta",GVCF_Files=c("one.g.vcf","two.g.vcf"),Output_Base_Name="Final_Variants")

MergeAndGenotype<-function(Reference_Genome="No Default", GVCF_Files="No Default", Output_Base_Name="Merged_GVCF",pathtoGATK='./gatk', pathtobash="C:/Program Files/Git/bin/bash.exe"){
# Check Input
if(file.exists(Reference_Genome)==FALSE){return(message("No existing Reference genome has been submitted"))
}else if(length(GVCF_Files)<=1){return(message("Submit at least two g.vcf files to merge."))}

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

