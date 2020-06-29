#' VariantToolsCall()
#'
#' This function allows the user to send commands to the GATK, running either natively on linux and Mac or within a docker for Windows, from within R.
#' @param Reference_Genome The reference genome unto which the sequence data has been mapped (including .fasta or .fna file extension).
#' @param Alignment_File_Name The file name of the alignment data from which variants are to be called (including .bam )
#' @param Output_Base_Name The desired output base name for the variant call format (VCF) file. [default="VariantTools_Call"]
#' @param minBaseQuality Minimum nucleotide quality below which a variant will be masked. [default=10]
#' @param minMapQuality Minimum mapping quality below which a variant will be masked. [default=10]
#' @param minDepth Minimum read depth below which a variant will be masked [default=2]
#' @param maxDepth Maximum read depth above which a variant will be masked [default=0.2]
#' @param p_lower From VariantTools: 'The lower bound on the binomial probability for a true variant.'
#' @param p_error From VariantTools: 'The binomial probability for a sequencing error (default is reasonable for Illumina data with the default quality cutoff).' [default=0.001]
#' @param read_count From VariantTools: 'Require at least this many high quality reads with the alternate base. The default value is designed to catch sequencing errors where coverage is too low to rely on the LRT. Increasing this value has a significant negative impact on power.' [default=2]
#' @keywords variant call VariantTools
#' @export
#' @examples
#' # Perform a simple call
#' VariantToolsCall(Reference_Genome="Ref.fasta",Alignment_File_Name="my_alignment.bam",Output_Base_Name="Variants")


VariantToolsCall<-function(Reference_Genome="No Default",Alignment_File_Name="No Default",Output_Base_Name="VariantToolsCall.vcf",minBaseQuality=10,minMapQuality=10,minDepth=2,maxDepth=Inf,p_lower=0.2,p_error=0.001,read_count=2){
# Check input
if(file.exists(Reference_Genome)==FALSE){return(message("No existing Reference genome has been submitted"))
}else if(file.exists(Alignment_File_Name)==FALSE){return(message("No existing alignment file has been submitted"))
}
  message("Locating Reference Bounderies")
# Fetch Reference
Ref<-Biostrings::readDNAStringSet(Reference_Genome)
Ref<-Biostrings::getSeq(Ref,names=Ref[1]@ranges@NAMES)

# Set Genomic Range to the entire genome
Genomic_Range<-GenomicRanges::GRanges(seqnames=Ref[1]@ranges@NAMES,
                                      ranges=IRanges(Ref[1]@ranges@start,Ref[1]@ranges@width))

# Register filter criteria
  message("Establishing Filtering Criteria")
param<-Rsamtools::ApplyPileupsParam(minBaseQuality=minBaseQuality,
                                    minMapQuality=minMapQuality,
                                    minDepth=minDepth,
                                    maxDepth=maxDepth,
                                    which=Genomic_Range,
                                    yieldAll=FALSE)

  message("Locating relevant information")
# Fetch BAM File
BAM<-Rsamtools::BamFile(Alignment_File_Name)
BAM<-Rsamtools::PileupFiles(Alignment_File_Name,param=param)

# Create a pileup file of those positions that passed the filters
  message("Filtering BAM File. This may take a while")
Pileup<-VariantTools::pileupVariants(BAM,Ref, param=param, baseOnly=TRUE)

# Filter
Calling_Filters<-VariantTools::VariantCallingFilters(p.lower=p_lower,
                                                     p.error=p_error,
                                                     read.count=read_count)
# Call variants
  message("Calling variants")
Called_Variants<-VariantTools::callVariants(Pileup,
                                            calling.filters=Calling_Filters)

# Save to a vcf file
  message("Saving to VCF File format")
vcfR::writeVcf(Called_Variants,Output_Base_Name)
  message("The variants have been saved as ",Output_Base_Name ," in ",getwd())
}
