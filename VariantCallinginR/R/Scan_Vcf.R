#' Scan_Vcf()
#'
#' Based on VariantAnnotation, this function allows the user to find the starting row and number of variants of all chromosomes in a vcf file. This function is meant to provide an overview of the number of variants called, as well as the nrows and skip arguments required when reading in a vcf per chromosome.
#' @param vcf_file The Variant Call Format (VCF) file name to scan (including .vcf extension). This function is incompatible with .gz compressed file.
#' @keywords exploration scan vcf
#' @export
#' @examples
#' # Scan a vcf file:
#' Scan_Vcf("file.vcf")

Scan_Vcf<-function(vcf_file="No Default"){
# Check input
if(file.exists(vcf_file)==F){return(message("No input vcf file has been submitted or the submitted file does not exist in the working directory."))
}else if(grepl(".gz",vcf_file)==1){return(message("Please decompress the vcf file before submission."))}

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
