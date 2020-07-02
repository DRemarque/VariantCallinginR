#
message("This is the AlignQual() R Server Script")
message("Starting up")
# Index Input
option_list = list(
  optparse::make_option(c("--vcf_list"), type="character", default="No Default",
                        help="The R object RData file name into which the vcf data was stored with readVcf.",metavar = "character"),
  optparse::make_option(c("--plot"), type="character", default="No Default",
                        help="If only a specific parameter has to be plotted, insert the name with plot='parameter name'. This only works for the parameters plotted when plot='No Default'. Other parameters give an empty plot.",metavar = "character"),
  optparse::make_option(c("--save_name"), type="character", default="VariantPlot",
                        help="The generated quality plot can be be saved as a pdf file by setting save to the desired base save name (without .pdf extension)",metavar = "character"),
  optparse::make_option(c("--server"), type="logical", default=FALSE,
                        help="To save the quality plots on a server without gpu access, set this variable to TRUE to use an alternative save method.",metavar="logical")
)
opt_parser = optparse::OptionParser(option_list=option_list)
opt= optparse::parse_args(opt_parser)

# Load RData
load(opt$vcf_list)

# Build function
VariantPlot<-function(vcf_list=x,plot=opt$plot,save_name=opt$save_name,server=opt$server){

# Check input
if(exists("vcf_list")==F){return(message("No vcf data has been submitted"))
  }else if(class(vcf_list)!="list"){return(message("The given object does not exist or is not of type list constructed with readVcf()"))}
# Load ggplot2
library(ggplot2)

# Plot available data and store into objects to later arrange
i<-1
x<-list()
# General information on position
if(plot=="No Default"){
  x[[i]]<-(ggplot(data=vcf_list[[1]],aes(x=POS))
           +geom_density(stat="bin",aes(),fill="dodgerblue3")
           +labs(title="Variant Distribution", y="Number of Variants",x="Position (bp)")
           +coord_cartesian(expand = FALSE))
  i<-i+1
  x[[i]]<-(ggplot()
           +geom_density(alpha=0.5,stat="bin",aes(x=vcf_list[[1]]$POS[which(vcf_list[[1]]$Type=="SNP")],fill="SNP"))
           +geom_density(alpha=0.5,stat="bin",aes(x=vcf_list[[1]]$POS[which(vcf_list[[1]]$Type=="MNP")],fill="MNP"))
           +scale_fill_manual(labels=c("INDELs","SNPs"),values=c("SNP"="dodgerblue2","MNP"="blue3"))
           +labs(title="Variant Distribution Per Type", y="Number of Variants",x="Position (bp)",fill="Type")
           +theme(legend.title = element_blank(),legend.position=c(0.9,0.86),legend.background=element_rect(fill=NA))
           +coord_cartesian(expand = FALSE))
  i<-i+1
  x[[i]]<-(ggplot()
           +geom_bar(stat="count",aes(x=as.factor(vcf_list[[1]]$Type)),fill=c("dodgerblue2","blue3"))
           +geom_text(stat="count",data=vcf_list[[1]], aes(x=as.factor(Type),colour="text",label=..count..), vjust=1.5)
           +scale_colour_manual(values = "white")
           +labs(title="Variant Type Numbers", y="Number of Variants",x="Type")
           +theme(legend.position = "none")
           )
}

# Non standard information
if((grep("QD",colnames(vcf_list[[1]]))!=0&plot=="No Default") | plot=="QD"){
  # Quality normalised to depth
  i<-i+1
  x[[i]]<-(ggplot(data=vcf_list[[1]],aes(x=QD))
           +geom_density(stat="bin",aes(),fill="darkorchid4")
           +labs(title="QD distribution", y="Number of Variants",x="QD score")
           +coord_cartesian(expand = FALSE))}
if((grep("FS",colnames(vcf_list[[1]]))!=0&plot=="No Default") | plot=="FS"){
  # Fisher Strand Bias
  i<-i+1
  x[[i]]<-(ggplot(data=vcf_list[[1]],aes(x=FS))
           +geom_density(stat="bin",aes(),fill="darkorchid3")
           +labs(title="Fisher Strand Bias", y="Number of Variants",x="FS score")
           +coord_cartesian(expand = FALSE))}
if((length(grep("MQ",colnames(vcf_list[[1]])))!=0&plot=="No Default") | plot=="MQ"){
  # Mapping Quality of the variant
  i<-i+1
  x[[i]]<-(ggplot(data=vcf_list[[1]],aes(x=MQ))
           +geom_density(stat="bin",aes(),fill="darkorchid2")
           +labs(title="Mapping Quality", y="Number of Variants",x="MQ score")
           +coord_cartesian(expand = FALSE))}
if((grep("MQRankSum",colnames(vcf_list[[1]]))!=0&plot=="No Default") | plot=="MQRankSum"){
  # Mapping Quality Ranksum
  i<-i+1
  x[[i]]<-(ggplot(data=vcf_list[[1]],aes(x=MQRankSum))
           +geom_density(stat="bin",aes(),fill="deeppink3")
           +labs(title="Mapping Quality of REF vs ALT Allele Reads", y="Number of Variants",x="MQRanksum")
           +coord_cartesian(expand = FALSE))}
if((grep("ReadPosRankSum",colnames(vcf_list[[1]]))!=0&plot=="No Default") | plot=="ReadPosRankSum"){
  # Read Position Ranksum
  i<-i+1
  x[[i]]<-(ggplot(data=vcf_list[[1]],aes(x=ReadPosRankSum))
           +geom_density(stat="bin",aes(),fill="deeppink2")
           +labs(title="Allele Position Within Reads", y="Number of Variants",x="ReadPosRankSum")
           +coord_cartesian(expand = FALSE))}
if((grep("SOR",colnames(vcf_list[[1]]))!=0&plot=="No Default") | plot=="SOR"){
  # Strand Odds Ratio
  i<-i+1
  x[[i]]<-(ggplot(data=vcf_list[[1]],aes(x=ReadPosRankSum))
           +geom_density(stat="bin",aes(),fill="deeppink1")
           +labs(title="Strand Odds Ratio", y="Number of Variants",x="Strand Odds Ratio score")
           +coord_cartesian(expand = FALSE))}

# Arrange plots:
y<-suppressMessages(ggpubr::ggarrange(plotlist = x))
# Save plots if desired
if(save_name!="No Default"&server==F){
  save_name<-paste(save_name,".pdf",sep="")
  pdf(file=file.path(getwd(),save_name), width=11.7, height=11.7)
  y
  dev.off()
}else if(save_name!="No Default"&server==T){
 save_name<-paste(save_name,".pdf",sep="")
 bitmap(file=file.path(getwd(),save_name), width=11.7, height=11.7, type="pdfwrite", res=100)
 y
 dev.off()
}
# Print plots
suppressMessages(y)

# End of function
}

# Execute plotting
VariantPlot()
