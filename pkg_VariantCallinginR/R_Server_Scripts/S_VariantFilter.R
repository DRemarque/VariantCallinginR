#
message("This is the AlignQual() R Server Script")
message("Starting up")
# Index Input
option_list = list(
  optparse::make_option(c("--vcf_list"), type="character", default="No Default",
                        help="The R object RData file name into which the vcf data was stored with readVcf.",metavar = "character"),
  optparse::make_option(c("--Output_Base_Name"), type="character", default="Filtered_Vcf",
                        help="The base output name for the filtered vcf R object. Saving to vcf is not yet supported",metavar = "character"),
  optparse::make_option(c("--General_filter"), type="character", default="No Default",
                        help="The filter criteria for general parameters, defined as a character string per criterium. View examples for details",metavar = "character"),
  optparse::make_option(c("--Sample_Specific_Filter"), type="character", default=0,
                        help="The filter criteria for sample parameters, defined as a character string per criterium. View examples for details",metavar = "character"),
  optparse::make_option(c("--min_dist"), type="integer", default=FALSE,
                        help="The minimum variant free space around any given variant in basepairs. Use this when the presence of nearby variants is undesirable, for example during variant genotyping primer construction. This filtering criterium takes all variants of the vcf into account. If only high quality variants should be taken into account, perform the filtering step twice.",metavar="integer"),
  optparse::make_option(c("--min_AD"), type="integer", default=FALSE,
                        help="The minimum summed allelic depth allowed. This is a measure of coverage used to discard variant with little evidende.",metavar="integer"),
  optparse::make_option(c("--max_AD"), type="integer", default=FALSE,
                        help="The maximum summed allelic depth allowed. This is a measure of oversequencing often observed for repeats and overabundant transposons.",metavar="integer")
)
opt_parser = optparse::OptionParser(option_list=option_list)
opt= optparse::parse_args(opt_parser)

# Load RData
  message("Loading data...")
load(opt$vcf_list)

# Build Function
VariantFilter<-function(vcf_list=x, General_filter=opt$General_filter,Sample_Specific_Filter=opt$Sample_Specific_Filter,min_dist=opt$min_dist,min_AD=opt$min_AD,max_AD=opt$max_AD){

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
  }
  names(pass)<-names(vcf_list)
  # End of function
  return(pass)
}

# Execute Filter
x<-VariantFilter()

# Save Filtered object
opt$Output_Base_Name<-paste(opt$Output_Base_Name,".RData",sep="")
save(x,file=opt$Output_Base_Name)
