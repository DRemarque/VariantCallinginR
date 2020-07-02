#' FastqQual()
#'
#' Based on ShortRead, SystempipeR and Rbowtie2, this function allows you to create quality plots and potential adapter sequences of both single and paired-end sequence data.
#' @param Forward_Fastq The single-read or forward paired-end fastq file name (including .fastq file extension)
#' @param Reverse_Fastq When Paired_End_Data=TRUE, the reverse paired-end fastq file name (including .fastq file extension). This argument is ignored when Paired_End_Data=FALSE.
#' @param Paired_End_Data The input data type, where FALSE corresponds with single-end reads. [default=TRUE]
#' @param Create_Graph Whether to create a quality graph when Paired-End_Data=TRUE. [default=TRUE]
#' @param Adapter_Check Whether to check for the presence of adapters when Paired-End_Data=TRUE. [default=TRUE]
#' @param Output_Base_Name The base name of the generated quality plots. No file extension is required. [default="FastqQual"]
#' @param server To save the quality plots on a server without gpu access, set this variable to TRUE to use an alternative save method.
#' @keywords fastq quality
#' @export
#' @examples
#' To do a general check of paired-end sequence data:
#' FastqQual(Forward_Fastq="example_fw.fastq", Reverse_Fastq="example_rv.fastq")
#' To do a general check of single-end sequence data:
#' FastqQual(Forward_Fastq="example.fastq", Paired_End_Data=FALSE)
#' To search only for paired-end sequence adapters:
#' FastqQual(Forward_Fastq="example_fw.fastq", Reverse_Fastq="example_rv.fastq",Create_Graph=FALSE)
#' To create only a quality graph of paired-end sequence data without searching for adapters:
#' FastqQual(Forward_Fastq="example_fw.fastq", Reverse_Fastq="example_rv.fastq",Adapter_Check=FALSE)

FastqQual<-function(Forward_Fastq="No Default", Reverse_Fastq="No Default", Paired_End_Data=TRUE, Adapter_Check=TRUE, Create_Graph=TRUE, Output_Base_Name="FastqQual",server=F){
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
