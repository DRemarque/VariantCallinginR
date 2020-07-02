#' VariantFilter()
#'
#' Based on vcfR, this function is used to read variant data from a variant call format (VCF) file into a comprehensive list of dataframes. This list is split into general parameters and sample specific parameters.
#' @param vcf_list R object into which vcf data was stored with ReadVcf()
#' @param General_filter The filter criteria for general parameters, defined as a character string per criterium. View examples for details
#' @param Sample_Specific_Filter The filter criteria for sample parameters, defined as a character string per criterium. View examples for details
#' @param min_dist The minimum variant free space around any given variant in basepairs. Use this when the presence of nearby variants is undesirable, for example during variant genotyping primer construction. This filtering criterium takes all variants of the vcf into account. If only high quality variants should be taken into account, perform the filtering step twice.
#' @param min_AD The minimum summed allelic depth allowed. This is a measure of coverage used to discard variant with little evidende.
#' @param max_AD The maximum summed allelic depth allowed. This is a measure of oversequencing often observed for repeats and overabundant transposons.
#' @keywords vcf vcfR filter
#' @export
#' @examples
#' # Select only SNPs (note the double == sign)
#' VariantFilter(vcf_list=vcf_object,General_Filter="Type=='SNP'"))
#' # Select homozygotes with MQ higher than 2
#' VariantFilter(vcf_list=vcf_object, General_Filter="MQ>2",Sample_Specific_Filter="GT=='1/1'")


VariantFilter<-function(vcf_list="No Default", General_filter="No Default",Sample_Specific_Filter="No Default",min_dist="No Default",min_AD=0,max_AD=Inf){

# Check input
if(exists("vcf_list")==F){return(message("No vcf data has been submitted"))
}else if(General_filter[1]=="No Default" & Sample_Specific_Filter[1]=="No Default"&min_AD==0 &max_AD==Inf&min_dist=="No Default"){return(message("No filtering criteria were given"))}

# Setup output
i<-1
f<-list()
filt<-list()
filt1<-list()

# Filter based on general and sample specific parameters
if(General_filter[1]!="No Default"&Sample_Specific_Filter[1]!="No Default"){
# Index specific filters
z<-2
while(z<=length(vcf_list)){ #While for sample
    f<-paste("vcf_list[[",z,"]]$",Sample_Specific_Filter,sep="")
    f<-paste("filt1[[",z-1,"]]<-which(",paste(f,collapse = " & "),")",sep = "")
    # Look for passing values
    eval(parse(text = f))
z<-z+1}

# Index general filters
f<-paste("vcf_list[[1]]$",General_filter,sep="")
f<-paste("filt[[",z,"]]<-which(",paste(f,collapse = " & "),")",sep = "")
# Look for passing values
eval(parse(text = f))
# Index passed variant row numbers
filt<-append(unlist(filt1),unlist(filt))
filt<-filt[duplicated(filt)] #Keep only duplicates
filt<-filt[!duplicated(filt)] #Remove extra duplicates

# Filter based on sample specific parameters only
}else if(General_filter[1]=="No Default"&Sample_Specific_Filter[1]!="No Default"){
  # Index specific filters
  z<-2
  while(z<=length(vcf_list)){
    f<-paste("vcf_list[[",z,"]]$",Sample_Specific_Filter,sep="")
  f<-paste("filt[[",z-1,"]]<-which(",paste(f,collapse = " & "),")",sep = "")
  # Look for passing values
  eval(parse(text = f))
  z<-z+1}
  # Index passed variant row numbers
  filt<-unlist(filt)
# Filter based on general parameters only
}else if(General_filter[1]!="No Default"&Sample_Specific_Filter[1]=="No Default"){
  # Index general filters
  f<-paste("vcf_list[[1]]$",General_filter,sep="")
  f<-paste("filt<-which(",paste(f,collapse = " & "),")",sep = "")
  # Look for passing values
  eval(parse(text = f))
  # Index passed variant row numbers
  filt<-unlist(filt)
# End of filter
}else{filt<-unlist(filt)}

# Check for nearby variants (for primer development)
if(min_dist!="No Default"){
  # Calculate distance between variants
  Dist_Raw<-diff(vcf_list$General_Parameters$POS)
  Dist_Raw1<-c(Dist_Raw,1000)
  Dist_Raw2<-c(1000,Dist_Raw)
  # Remove variants that are too close to each other
  bp_filt<-which(Dist_Raw1>=min_dist &
                 Dist_Raw2>=min_dist)
  # Add the passing variants to the filtering list (filt)
  bp_filt<-unlist(bp_filt)
  bp_filt<-bp_filt[duplicated(bp_filt)]
  bp_filt<-bp_filt[!duplicated(bp_filt)]
  filt<-append(filt,bp_filt)
  filt<-filt[duplicated(filt)] #Keep only duplicates
  filt<-filt[!duplicated(filt)] #Remove extra duplicates
}

# Filter on allelic depth
if(min_AD!=0 | max_AD!=Inf){
  f<-list()
  z<-2
  while (z<=length(vcf_list)) {
    AD_Sum<-0
    tX<-paste("AD_Sum<-c(AD_Sum,",gsub(",","+",vcf_list[[z]]$AD),")",sep="")
    eval(parse(text=tX))
    AD_Sum<-AD_Sum[-1]
    f[[z-1]]<-which(AD_Sum>=min_AD & AD_Sum<=max_AD)
    z<-z+1
  }
  f<-unlist(f)
  f<-f[duplicated(f)]
  f<-f[!duplicated(f)]
  filt<-append(filt,f)
  filt<-filt[duplicated(filt)] #Keep only duplicates
  filt<-filt[!duplicated(filt)] #Remove extra duplicates
}

# Store only variants that passed all filters
i<-1
pass <- list()
while(i<=length(vcf_list)){
pass[[i]]<-vcf_list[[i]][filt,]
i<-i+1
names(pass)<-names(vcf_list)
}
  # End of function
  return(pass)
}
