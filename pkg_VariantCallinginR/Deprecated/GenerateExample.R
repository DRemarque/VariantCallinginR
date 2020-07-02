#' GenerateExample()
#'
#' This function allows you to generate a small example dataset of either single or paired-end sequenced DNA
#' #' @param Number_of_reads The amount of reads sampled from the entire example genome, stored in a fastq file. [default=1000]
#' #' @param Number_of_variants The amount of variants hidden in the example data. This will be rounded to a number dividable by 4 in order to equally distribute the variants per nucleotide type. Half of the variants will be heterozygous. [default=40]
#' #' @param Length_of_genome The basepair length of the example genome. Choosing a large number may cause the function to become extremely slow. [default=1000]
#' #' @param Read_length The basepair length of a single read. The default corresponds newer with Illumina read lengths [default=150]
#' #' @param Paired_End_Data Whether to generate single or paired end data, where FALSE corresponds with single-ends [default=TRUE]
#' @keywords example data
#' @export
#' @examples
#' To create old paired end Illumina reads of a single 1500 bp gene with 8 hidden variants (10/4 rounded gives 8):
#' GenerateExample(Number_of_reads=500,Number_of_variants=10,Length_of_genome=1500,Read_length=100,Paired_End_Data=T)

GenerateExample<-function(Number_of_reads=1000,Number_of_variants=40,Length_of_genome=1000,Read_length=150,Paired_End_Data=TRUE){
  # Construct S4 Class object to save output to
  setClass(Class="Example_Data",
           representation(reference="character",
                          alternative="character"))
  # Make read length even
  if(grepl("\\.",Read_length/2)==T){Read_length<-Read_length+1}
  # Construct a Reference
  Reference_Genome<-paste(sample(Biostrings::DNA_ALPHABET[1:4], size=Length_of_genome, replace=TRUE),collapse="")
  # Save Reference
  seqinr::write.fasta(Reference_Genome,"example_reference",file.out = "example.fasta")
    message("The example reference genome has been saved as example.fasta in ",getwd())
    message("")
  # Construct and Alternative Genome
  Alt_Genome<-unlist(strsplit(Reference_Genome,""))
  # Select random postions for the variants and convert half into heterozygous variants
  Number_of_variants<-round(Number_of_variants/4)*4
  Variant_position<-c(sample(which(Alt_Genome!="A")[which(which(Alt_Genome!="A")>=Read_length+1 &                                       which(Alt_Genome!="A")<=Length_of_genome-Read_length+1)],
                             Number_of_variants/4),
                      sample(which(Alt_Genome!="T")[which(which(Alt_Genome!="T")>=Read_length+1 &                                       which(Alt_Genome!="T")<=Length_of_genome-Read_length+1)],
                             Number_of_variants/4),
                      sample(which(Alt_Genome!="G")[which(which(Alt_Genome!="G")>=Read_length+1 &                                       which(Alt_Genome!="G")<=Length_of_genome-Read_length+1)],
                             Number_of_variants/4),
                      sample(which(Alt_Genome!="C")[which(which(Alt_Genome!="C")>=Read_length+1 &                                       which(Alt_Genome!="C")<=Length_of_genome-Read_length+1)],
                             Number_of_variants/4))
  Heterozygous_variants<-sample(Variant_position,Number_of_variants/2)
  # Insert the variants into the chosen basepair positions
    message("Inserting variants...")
    message("")
  Alt_Genome[Variant_position[1:(Number_of_variants/4)]]<-"A"
  Alt_Genome[Variant_position[((Number_of_variants/4)+1):((Number_of_variants/4)*2)]]<-"T"
  Alt_Genome[Variant_position[((Number_of_variants/4)*2+1):((Number_of_variants/4)*3)]]<-"G"
  Alt_Genome[Variant_position[((Number_of_variants/4)*3+1):Number_of_variants]]<-"C"
  # Store the variant positions in Example_Variants.txt
  ALT<-c(replicate(Number_of_variants/4, "A"),replicate(Number_of_variants/4, "T"),
         replicate(Number_of_variants/4, "G"),replicate(Number_of_variants/4, "C"))
  Type<-replicate(Number_of_variants, "homozygous")

  Example_Variants<-cbind(Variant_position, ALT, Type)
  Example_Variants[which(is.element(Example_Variants[,1],Heterozygous_variants)),3]<-"heterozygous"
  colnames(Example_Variants)<-c("Pos","ALT","Type")

  write.table(Example_Variants,"Example_Variants.txt",sep="\t",row.names=FALSE)
    message("The variant postions and types have been saved in Example_Variants.txt")
    message("")
  # Collapse the alternative genome
  Alt_Genome<-paste(Alt_Genome,collapse="")
  # Generate read qualities between 10 and 40
  fw_quality<-replicate(Number_of_reads,{paste(sample(unlist(strsplit(rawToChar(as.raw(64:73)),"")),Read_length,replace=TRUE), collapse="")})
  rv_quality<-replicate(Number_of_reads,{paste(sample(unlist(strsplit(rawToChar(as.raw(64:73)),"")),Read_length,replace=TRUE), collapse="")})
  # Sample forward sequence reads
    message("Creating forward reads...")
    message("")
  start<-sample(1:(Length_of_genome-Read_length),(length(fw_quality)-Number_of_variants),replace=TRUE)
  end<-start+Read_length-1 ; fw_reads<-list() ; i<-1 ; z<-1

  repeat{
    if(i<=(length(fw_quality)-Number_of_variants)){
      fw_reads[i]<-substr(Alt_Genome,start=start[i],stop=end[i])
      i<-i+1
    } else if(i>=(length(fw_quality)-Number_of_variants) & i<=Number_of_reads){
      fw_reads[i]<-substr(Reference_Genome,start=Heterozygous_variants[z]-Read_length/2,
                          stop=Heterozygous_variants[z]+(Read_length/2)-1)
      fw_reads[i+1]<-substr(Reference_Genome,start=Heterozygous_variants[z]-Read_length/2,
                            stop=Heterozygous_variants[z]+(Read_length/2)-1)
      i<-i+2;z<-z+1
    } else{rm(i,z,start,end)
      fw_reads<-unlist(fw_reads)
      break}
    }
  # write a fastq record
    message("Saving Reads...")
  names<-paste("example", 1:Number_of_reads)
  suppressWarnings(ChIPsim::writeFASTQ(fw_reads, fw_quality,name=names, file="example_fw.fastq"))
    message("Forward reads have been saved in example_fw.fastq")

  if(Paired_End_Data==TRUE){
      message("")
      message("Generating reverse reads...")
      message("")
    i<-1
    rv_reads<-list()
    # Reverse Reads
    while (i<=Number_of_reads) {
      rv_reads[[i]]<-paste(unlist(rev(strsplit(fw_reads[i],"")[[1]])),collapse="")
      i<-i+1}
    # write a fastq record
      message("Saving Reads...")
    suppressWarnings(ChIPsim::writeFASTQ(rv_reads, rv_quality,"example", file="example_rv.fastq"))
      message("Reverse reads have been saved in example_rv.fastq")

    # Load example files
    Alternative_Genome<-c("example_fw.fastq","example_rv.fastq")
    Reference_Genome<-"example.fasta"
  }else{
    # Load example files
    Alternative_Genome<-c("example_fw.fastq")
    Reference_Genome<-"example.fasta"
  }
  # Save output in S4 Class object
  return(new("Example_Data",
             reference=Reference_Genome,
             alternative=Alternative_Genome))
}
