#
message("This is the ScanVcf() R Server Script")
message("Starting up")
# Index Input
option_list = list(
  optparse::make_option(c("--vcf_file"), type="character", default="No Default",
                        help="The Variant Call Format (VCF) file name to scan (including .vcf extension). This function is incompatible with .gz compressed file.",metavar = "character"),
  optparse::make_option(c("--Output_Base_Name"), type="character", default="ScanVcf_Statistics",
                        help="Base name of the output file containing the statistics collected with Scan_Vcf",metavar="character")
)
opt_parser = optparse::OptionParser(option_list=option_list)
opt= optparse::parse_args(opt_parser)

# Build Function
Scan_Vcf<-function(vcf_file=opt$vcf_file){
# Check input
if(file.exists(vcf_file)==F){return(message("No input vcf file has been submitted or the submitted file does not exist in the working directory."))
}else if(grepl(".gz",vcf_file)==1){return(message("Please decompress the vcf file before submission."))}

# Find chromosome names
hdr<-suppressWarnings(VariantAnnotation::scanVcfHeader(vcf_file))

# Store as dataframe
Chr_list<-data.frame(hdr@header@listData$contig@rownames)
names(Chr_list)<-"Contig Name"

# Find sample names and numbers
Samples<-as.data.frame(hdr@samples);colnames(Samples)<-("Sample_Name")
Sample_Nr<-as.numeric(rownames(Samples))+9
Samples<-cbind(Samples,Sample_Nr)
# Output the dataframe
output<-list(Chr_list,Samples)
names(output)<-c("Chromosomes","Samples")
return(output)
}

# Execute
x<-Scan_Vcf()

# Save
opt$Output_Base_Name<-paste(opt$Output_Base_Name,".txt",sep="")
cat(capture.output(print(x), file=opt$Output_Base_Name))
