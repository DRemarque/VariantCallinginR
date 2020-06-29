#
message("This is the SaveVcfObject() R Server Script")
message("Starting up")
# Index Input
option_list = list(
  optparse::make_option(c("--vcf_list"), type="character", default="No Default",
                        help="The R object RData file name into which the vcf data was stored with readVcf.",metavar = "character"),
  optparse::make_option(c("--Output_Base_Name"), type="character", default="VcfObject",
                        help="The desired output base name for the created excel file. [default='VcfObject']",metavar = "character")
)
opt_parser = optparse::OptionParser(option_list=option_list)
opt= optparse::parse_args(opt_parser)

# Load RData
load(opt$vcf_list)

# Build Function
SaveVcfObject<-function(vcf_list=x,Output_Base_Name=opt$Output_Base_Name){
  # Check input
  if(file.exists(vcf_list)=="No Default"){return(message("No vcf R object has been submitted"))
  }
  Output_Base_Name<-paste(Output_Base_Name,".xlsx",sep="")
  # General parameters
  message("Saving General Parameters")
  xlsx::write.xlsx2(vcf_list[[1]], Output_Base_Name, sheetName = "General Parameters",
                    col.names = TRUE, row.names = FALSE, append = FALSE)
  # Sample Specific parameters
  message("Saving Sample Specific Parameters")
  i<-2
  while(i<=length(vcf_list)){
    xlsx::write.xlsx2(vcf_list[[i]], Output_Base_Name, sheetName = names(vcf_list)[i],
                      col.names = TRUE, row.names = FALSE, append = TRUE)
    i<-i+1
  }

  # End of function
  message("The vcf R object has been saved as ",Output_Base_Name," in ",getwd())
}

# Execute saving
SaveVcfObject()
