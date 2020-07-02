#' FetchFlanks()
#'
#' This function is used to retrieve the sequences surrounding the variants for primer construction. 
#' @param vcf_list R object into which vcf data was stored with ReadVcf()
#' @param Reference_Genome The reference genome unto which the sequence data has been mapped (including .fasta or .fna file extension).
#' @param validation Whether or not the submitted reference genome is compared to data from the vcf_list to verify that the correct reference has been submitted. This prevents the function from potentially generating flanking regions for non-existing alleles.
#' @param Flank_Length How many base pairs left and right of the variant are to be extracted. The default is 75 bp which produces a sequence of 150 bp (75 left of the variant and 75 right of the variant).
#' @param Output_Base_Name To save the flanking regions directly to an excel file, insert the desired base file name here. The file extension is automatically added to this name.
#' @keywords vcf vcfR flanks flanking region primer construction
#' @export
#' @examples
#' # Retrieve a short region of 10 bp and save it to an excel file:
#' FetchFlanks(vcf_list=vcf_object,Reference_Genome="MyRef.fasta",Flank_Length=10,Output_Base_Name="Regions"))


# Build function
FetchFlanks<-function(vcf_list="No Default", Reference_Genome="No Default",validation=TRUE, Flank_Length=75,Output_Base_Name="No Default"){
# Check input
if(file.exists(Reference_Genome)==FALSE){return(message("No existing Reference genome has been submitted"))
}
  
#Extract the sequence data
  message("Reading reference...")
seq<-seqinr::read.fasta(Reference_Genome)
  
if(validation==T){ 
  # As a control, the ref alleles are called from the chromosome to see if these are the same as the ones stated in the vcf file. This should reveal whether the used fasta file is correct
    message("Verifying reference genome input")
  # Walk through all chromosomes
  i<-1
  while(i<=length(levels(as.factor(vcf_list[[1]]$CHROM)))){
    ref<-seq[[which(names(seq)==levels(as.factor(vcf_list[[1]]$CHROM))[i])]][vcf_list[[1]]$POS[which(vcf_list[[1]]$CHROM==levels(as.factor(vcf_list[[1]]$CHROM))[i])]]
    i<-i+1}
    
if(length(ref[is.element(tolower(ref),tolower(vcf_list[[1]]$REF))])==length(vcf_list[[1]]$REF)){
   message("Reference verified: The reference genome corresponds to the vcf data.")
} else{return(message("Information of the reference and vcf file does not seem to correspond. Please check that the correct file input was given."))}
}

# Calculate the start and ends of the flanking regions.
L_start<-vcf_list[[1]]$POS-Flank_Length
L_end<-vcf_list[[1]]$POS-1
R_start<-vcf_list[[1]]$POS+1
R_end<-vcf_list[[1]]$POS+Flank_Length

# Merge the limits for parsing
L<-paste(L_start,L_end,sep=":")
R<-paste(R_start,R_end,sep=":")

# Create the command
z<-1:length(vcf_list$General_Parameters$POS)
commandL<-paste("Left_Flank[[",z,"]]<-toupper(paste(seq[[which(names(seq)=='",vcf_list[[1]]$CHROM,"')]][",L,"], collapse=''))",sep="")
commandR<-paste("Right_Flank[[",z,"]]<-toupper(paste(seq[[which(names(seq)=='",vcf_list[[1]]$CHROM,"')]][",R,"], collapse=''))",sep="")

# Create empty vectors that will be filled with the flanks
Left_Flank<-list()
Right_Flank<-list()
  
# Parse the flanking calls
eval(parse(text = commandL))
eval(parse(text = commandR))  
  
# Combine left and right flanks to receive the entire region
Sequence<-paste(unlist(Left_Flank),"[",vcf_list[[1]]$REF,"/",vcf_list[[1]]$ALT,"]",unlist(Right_Flank),sep="")
  
# Add Variant Name
Region_Name<-paste(vcf_list[[1]]$CHROM,vcf_list[[1]]$POS,sep="_")
Flanks<-cbind(Region_Name,Sequence)
  
  
# Optional: Save to excel
if(Output_Base_Name!="No Default"){
Output_Base_Name<-paste(Output_Base_Name,".xlsx",sep="")
openxlsx::write.xlsx2(Flanks,file=Output_Base_Name,sheetName="Flanking Regions",
                      col.names=T,row.names=F,append=F)}
  
# End of function
return(Flanks)
}
