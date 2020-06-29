#' VariantPlot()
#'
#' Based on ggplot2, this function is used to plot relevant quality parameters of variants stored in an R object with ReadVcf()
#' @param vcf_list R object into which vcf data was stored with ReadVcf()
#' @param plot If only a specific parameter has to be plotted, insert the name with plot='parameter name'. This only works for the parameters plotted when plot='No Default'. Other parameters give an empty plot.
#' @param save_name The generated quality plot can be be saved as a pdf file by setting save to the desired base save name (without .pdf extension)
#' @param server To save the quality plots on a server without gpu access, set this variable to TRUE to use an alternative save method.
#' @keywords vcf vcfR filter
#' @export
#' @examples
#' # Plot only the quality to depth ratio
#' VariantPlot(vcf_list=vcf_object,plot="QD")
#' # Plot all graphs and save the plot on a server
#' VariantPlot(vcf_list=vcf_object,save_name="Graph",server=T)


VariantPlot<-function(vcf_list="No Default",plot="No Default",save_name="No Default",server=F){

# Check input
if(exists("vcf_list")==F){return(message("No vcf data has been submitted"))
  }else if(class(vcf_list)!="list"){return(message("The given object does not exist or is not of type list constructed with readVcf()"))}


# Plot available data and store into objects to later arrange
i<-1
y<-list()
# General information on position
if(plot=="No Default"){
  y[[i]]<-(ggplot2::ggplot(data=vcf_list[[1]],ggplot2::aes(x=POS))
           +ggplot2::geom_density(stat="bin",ggplot2::aes(),fill="dodgerblue3")
           +ggplot2::labs(title="Variant Distribution", y="Number of Variants",x="Position (bp)")
           +ggplot2::coord_cartesian(expand = FALSE))
  i<-i+1
  y[[i]]<-(ggplot2::ggplot()
           +ggplot2::geom_density(alpha=0.5,stat="bin",ggplot2::aes(x=vcf_list[[1]]$POS[which(vcf_list[[1]]$Type=="SNP")],fill="SNP"))
           +ggplot2::geom_density(alpha=0.5,stat="bin",ggplot2::aes(x=vcf_list[[1]]$POS[which(vcf_list[[1]]$Type=="MNP")],fill="MNP"))
           +ggplot2::scale_fill_manual(labels=c("INDELs","SNPs"),values=c("SNP"="dodgerblue2","MNP"="blue3"))
           +ggplot2::labs(title="Variant Distribution Per Type", y="Number of Variants",x="Position (bp)",fill="Type")
           +ggplot2::theme(legend.title = ggplot2::element_blank(),legend.position=c(0.9,0.86),legend.background=ggplot2::element_rect(fill=NA))
           +ggplot2::coord_cartesian(expand = FALSE))
  i<-i+1
  y[[i]]<-(ggplot2::ggplot()
           +ggplot2::geom_bar(stat="count",ggplot2::aes(x=as.factor(vcf_list[[1]]$Type)),fill=c("dodgerblue2","blue3"))
           +ggplot2::geom_text(stat="count",data=vcf_list[[1]], ggplot2::aes(x=as.factor(Type),colour="text",label=..count..), vjust=1.5)
           +ggplot2::scale_colour_manual(values = "white")
           +ggplot2::labs(title="Variant Type Numbers", y="Number of Variants",x="Type")
           +ggplot2::theme(legend.position = "none")
           )
}

# Non standard information
if((length(grep("MQ",colnames(vcf_list[[1]])))!=0&plot=="No Default") | plot=="MQ"){
  # Mapping Quality of the variant
  i<-i+1
  y[[i]]<-(ggplot2::ggplot(data=vcf_list[[1]],ggplot2::aes(x=MQ))
           +ggplot2::geom_density(stat="bin",ggplot2::aes(),fill="darkorchid4")
           +ggplot2::labs(title="Mapping Quality", y="Number of Variants",x="MQ score")
           +ggplot2::coord_cartesian(expand = FALSE))}
if((grep("MQRankSum",colnames(vcf_list[[1]]))!=0&plot=="No Default") | plot=="MQRankSum"){
  # Mapping Quality Ranksum
  i<-i+1
  y[[i]]<-(ggplot2::ggplot(data=vcf_list[[1]],ggplot2::aes(x=MQRankSum))
           +ggplot2::geom_density(stat="bin",ggplot2::aes(),fill="darkorchid3")
           +ggplot2::labs(title="Mapping Quality of REF vs ALT Allele Reads", y="Number of Variants",x="MQRanksum")
           +ggplot2::coord_cartesian(expand = FALSE))}
if((grep("QD",colnames(vcf_list[[1]]))!=0&plot=="No Default") | plot=="QD"){
  # Quality normalised to depth
  i<-i+1
  y[[i]]<-(ggplot2::ggplot(data=vcf_list[[1]],ggplot2::aes(x=QD))
           +ggplot2::geom_density(stat="bin",ggplot2::aes(),fill="darkorchid2")
           +ggplot2::labs(title="QD distribution", y="Number of Variants",x="QD score")
           +ggplot2::coord_cartesian(expand = FALSE))}
if((grep("FS",colnames(vcf_list[[1]]))!=0&plot=="No Default") | plot=="FS"){
  # Fisher Strand Bias
  i<-i+1
  y[[i]]<-(ggplot2::ggplot(data=vcf_list[[1]],ggplot2::aes(x=FS))
           +ggplot2::geom_density(stat="bin",ggplot2::aes(),fill="deeppink3")
           +ggplot2::labs(title="Fisher Strand Bias", y="Number of Variants",x="FS score")
           +ggplot2::coord_cartesian(expand = FALSE))}
if((grep("ReadPosRankSum",colnames(vcf_list[[1]]))!=0&plot=="No Default") | plot=="ReadPosRankSum"){
  # Read Position Ranksum
  i<-i+1
  y[[i]]<-(ggplot2::ggplot(data=vcf_list[[1]],ggplot2::aes(x=ReadPosRankSum))
           +ggplot2::geom_density(stat="bin",ggplot2::aes(),fill="deeppink2")
           +ggplot2::labs(title="Allele Position Within Reads", y="Number of Variants",x="ReadPosRankSum")
           +ggplot2::coord_cartesian(expand = FALSE))}
if((grep("SOR",colnames(vcf_list[[1]]))!=0&plot=="No Default") | plot=="SOR"){
  # Strand Odds Ratio
  i<-i+1
  y[[i]]<-(ggplot2::ggplot(data=vcf_list[[1]],ggplot2::aes(x=ReadPosRankSum))
           +ggplot2::geom_density(stat="bin",ggplot2::aes(),fill="deeppink1")
           +ggplot2::labs(title="Strand Odds Ratio", y="Number of Variants",x="Strand Odds Ratio score")
           +ggplot2::coord_cartesian(expand = FALSE))}

# Arrange plots:
z<-suppressMessages(ggpubr::ggarrange(plotlist = y))
# Save plots if desired
if(save_name!="No Default"&server==F){
  save_name<-paste(save_name,".pdf",sep="")
  pdf(file=file.path(getwd(),save_name), width=11.7, height=11.7)
  z
  dev.off()
}else if(save_name!="No Default"&server==T){
 save_name<-paste(save_name,".pdf",sep="")
 bitmap(file=file.path(getwd(),save_name), width=11.7, height=11.7, type="pdfwrite", res=100)
 z
 dev.off()
}
if(server==F){
# Print plots
suppressMessages(z)}

# End of function
}
