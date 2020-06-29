#
message("This is the VariantToolsCall() R Server Script")
message("Starting up")
# Index Input
option_list = list(
  optparse::make_option(c("--Reference_Genome"), type="character", default="No Default",
                        help="The reference genome unto which the sequence data has been mapped (including .fasta or .fna file extension).",metavar = "character"),
  optparse::make_option(c("--Alignment_File_Name"), type="character", default="No Default",
                        help="The file name of the alignment data from which variants are to be called (including .bam )",metavar = "character"),
  optparse::make_option(c("--Output_Base_Name"), type="character", default="VariantTools_Call",
                        help="The desired output base name for the variant call format (VCF) file. [default='VariantTools_Call']'",metavar = "character"),
  optparse::make_option(c("--minBaseQuality"), type="integer", default=10,
                        help="Minimum nucleotide quality below which a variant will be masked. [default=10]",metavar="integer"),
  optparse::make_option(c("--minMapQuality"), type="integer", default=10,
                        help="Minimum mapping quality below which a variant will be masked. [default=10]",metavar = "integer"),
  optparse::make_option(c("--minDepth"), type="integer", default=2,
                        help="Minimum read depth below which a variant will be masked [default=2]",metavar = "integer"),
  optparse::make_option(c("--maxDepth"), type="double", default=Inf,
                        help="Maximum read depth above which a variant will be masked [default=Inf]",metavar = "double"),
  optparse::make_option(c("--p_lower"), type="double", default=0.2,
                        help="From VariantTools: 'The lower bound on the binomial probability for a true variant.'",metavar = "double"),
  optparse::make_option(c("--p_error"), type="double", default=0.001,
                        help="From VariantTools: 'The binomial probability for a sequencing error (default is reasonable for Illumina data with the default quality cutoff).' [default=0.001]",metavar = "double"),
  optparse::make_option(c("--read_count"), type="double", default=2,
                        help="From VariantTools: 'Require at least this many high quality reads with the alternate base. The default value is designed to catch sequencing errors where coverage is too low to rely on the LRT. Increasing this value has a significant negative impact on power.' [default=2]",metavar = "double")
)
opt_parser = optparse::OptionParser(option_list=option_list)
opt= optparse::parse_args(opt_parser)

# Build function
VariantToolsCall<-function(Reference_Genome=opt$Reference_Genome,Alignment_File_Name=opt$Alignment_File_Name,Output_Base_Name=opt$Output_Base_Name,minBaseQuality=opt$minBaseQuality,minMapQuality=opt$minMapQuality,minDepth=2,maxDepth=opt$maxDepth,p_lower=opt$p_lower,p_error=opt$p_error,read_count=opt$read_count){
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
