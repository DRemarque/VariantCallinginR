#
message("This is the FastqQual() R Server Script")
message("Starting up")
# Index Input
option_list = list(
  optparse::make_option(c("--Forward_Fastq"), type="character", default="No Default",
                        help="The single-read or forward paired-end fastq file name (including .fastq file extension)",metavar = "character"),
  optparse::make_option(c("--Reverse_Fastq"), type="character", default="No Default",
                        help="When Paired_End_Data=TRUE, the reverse paired-end fastq file name (including .fastq file extension). This argument is ignored when Paired_End_Data=FALSE.",metavar = "character"),
  optparse::make_option(c("--Paired_End_Data"), type="logical", default=FALSE,
                        help="The input data type, where FALSE corresponds with single-end reads. [default=FALSE]",metavar = "logical"),
  optparse::make_option(c("--Adapter_Check"), type="logical", default=TRUE,
                        help="Whether to check for the presence of adapters when Paired-End_Data=TRUE. [default=TRUE]",metavar="logical"),
  optparse::make_option(c("--Create_Graph"), type="logical", default=TRUE,
                        help="Whether to create a quality graph when Paired-End_Data=TRUE. [default=TRUE]",metavar="logical"),
  optparse::make_option(c("--Output_Base_Name"), type="character", default="FastqQual",
                        help="The base name of the generated quality plots. No file extension is required. [default='FastqQual']",metavar="logical"),
  optparse::make_option(c("--server"), type="logical", default=TRUE,
                        help="To save the quality plots on a server without gpu access, set this variable to TRUE to use an alternative save method.",metavar="logical")
)
opt_parser = optparse::OptionParser(option_list=option_list)
opt= optparse::parse_args(opt_parser)

# Build function
FastqQual<-function(Forward_Fastq=opt$Forward_Fastq, Reverse_Fastq=opt$Reverse_Fastq, Paired_End_Data=opt$Paired_End_Data, Adapter_Check=opt$Adapter_Check, Create_Graph=opt$Create_Graph, Output_Base_Name=opt$Output_Base_Name,server=opt$server){
  # Check Input
  if(Forward_Fastq=="No Default"){return(message("No Fastq file has been submitted"))
  }else if(Reverse_Fastq=="No Default" & Paired_End_Data==TRUE){return(message("No reverse fastq file has been submitted, but paired_end_data=TRUE"))
  }else if(Paired_End_Data==TRUE & Adapter_Check==FALSE & Create_Graph==FALSE){return(message("Both adapter search and graph creation have been turned off. No analysis will be performed."))}

  # Create Single Read Quality Plots
  if(Paired_End_Data==FALSE){
    message("Generating Single Read Quality Plots. This may take a while")
    message("")
  QA_Graph<-ShortRead::report(ShortRead::qa(Forward_Fastq,type="fastq"))
  QA_Graph<-gsub("/index.html","",QA_Graph)
  file.rename(from=QA_Graph, to=file.path(getwd(),Output_Base_Name))
    message("Quality plots have been saved in ", getwd()," in the ", Output_Base_Name, " directory")

    }else{

  # Create Paired End Read Quality Plots
  if(Create_Graph==TRUE){
  # Make a named character vector and fetch relevant information
    message("Collecting data. This may take a while")
    message("")
  faq_list<-c(file.path(getwd(),Forward_Fastq),file.path(getwd(),Reverse_Fastq))
  names(faq_list)<-c("M1A","M1B")
  fqlist<-systemPipeR::seeFastq(fastq=faq_list, batchsize = 100)

  if(server==F){
  # Plot the relevant information
    message("Generating Quality Plots. This may take a while")
    message("")
  Output_Base_NameX<-paste(Output_Base_Name,".pdf",sep="")
  pdf(file=file.path(getwd(),Output_Base_NameX), width=11.7, height=11.7)
  systemPipeR::seeFastqPlot(fqlist)
  dev.off()
  rm(faq_list,fqlist)
    message("Quality plots have been saved in ", getwd(), " as ", Output_Base_NameX)
    message("")
  }else if(server==T){
    message("Generating Quality Plots. This may take a while")
    message("")
  Output_Base_NameX<-paste(Output_Base_Name,".pdf",sep="")
  bitmap(file=file.path(getwd(),Output_Base_NameX), width=11.7, height=11.7, type="pdfwrite", res=100)
  systemPipeR::seeFastqPlot(fqlist)
  dev.off()
  rm(faq_list,fqlist)
    message("Quality plots have been saved in ", getwd(), " as ", Output_Base_NameX)
    message("")

    }
  }

  # Find Adapters in Unix and Mac OS
  if(Adapter_Check==TRUE & .Platform$OS.type!="windows"){
    message("Searching for potential adapter sequences")
  Output_Base_NameX<-paste(Output_Base_Name,".adapter1",sep="")
  file.create(Output_Base_NameX)
  Output_Base_NameX<-paste(Output_Base_Name,".adapter2",sep="")
  file.create(Output_Base_NameX)
  Adapters<-suppressWarnings(Rbowtie2::identify_adapters(file1=Forward_Fastq,
                                                         file2=Reverse_Fastq,
                                                         basename=file.path(getwd(),Output_Base_Name),
                                                         overwrite=TRUE))

  # Save Info in a .txt file (optional)
  Adapter<-as.data.frame(Adapters)
  Output_Base_NameX<-paste(Output_Base_Name,"_Adapters.txt",sep="")
  row.names(Adapter)<-c("Forward Adapter", "Reverse Adapter")
  write.table(Adapter,Output_Base_NameX,sep="\t")
  message("Adapter sequences have been saved in ",Output_Base_NameX)
  }else if(Adapter_Check==TRUE & .Platform$OS.type=="windows"){

  # Find Adapters in Windows OS
    message("Searching for potential adapter sequences")
  Output_Base_NameX<-paste(Output_Base_Name,".adapter1",sep="")
  file.create(Output_Base_NameX)
  Output_Base_NameX<-paste(Output_Base_Name,".adapter2",sep="")
  file.create(Output_Base_NameX)

  command<-paste(sep="","'",.libPaths()[1],"/Rbowtie2/AdapterRemoval' --identify-adapters --file1 ",Forward_Fastq," --file2 ",Reverse_Fastq," --basename ",Output_Base_Name)
  Adapters<-system(command, intern=TRUE)

  # Save Info in a .txt file (optional)
  Adapter<-as.data.frame(Adapters)
  Output_Base_NameX<-paste(Output_Base_Name,"_Adapters.txt",sep="")
  row.names(Adapter)<-c("Forward Adapter", "Reverse Adapter")
  write.table(Adapter,Output_Base_NameX,sep="\t")
    message("Adapter sequences have been saved in ",Output_Base_NameX)

  }

  # Clean up temporary files
  Output_Base_NameX<-paste(Output_Base_Name,".adapter1",sep="")
  file.remove(Output_Base_NameX)
  Output_Base_NameX<-paste(Output_Base_Name,".adapter2",sep="")
  file.remove(Output_Base_NameX)
  rm(Output_Base_NameX)
}
}

#Execute
FastqQual()
