#
message("This is the ReadVcf() R Server Script")
message("Starting up")
# Index Input
option_list = list(
  optparse::make_option(c("--vcf_file"), type="character", default="No Default",
                        help="The file to be read in (including .vcf file extension)",metavar = "character"),
  optparse::make_option(c("--contig"), type="character", default="No Default",
                        help="The name of the contig that has to be read in. For example 'C1' or 'NC_044370.1'. Note that nothing is read in if the contig name is misspelled and only one contig name can be submitted",metavar = "character"),
  optparse::make_option(c("--nrows"), type="double", default=-1,
                        help="The number of rows. The Scan_Vcf function returns the nrows per chromosome.",metavar = "double"),
  optparse::make_option(c("--skip"), type="integer", default=0,
                        help="The number of rows to skip before starting to read. The Scan_Vcf function returns the skip number corresponding with the chromosome of interest.",metavar = "integer"),
  optparse::make_option(c("--sample_nr"), type="integer", default=0,
                        help="The number(s) of the sample(s) to be included. The first sample always has sample_nr 10 and the default sample_nr 0 reads in all samples present. Specific sample numbers are returned with the Scan_Vcf function",metavar = "integer"),
  optparse::make_option(c("--verbose"), type="logical", default=TRUE,
                        help="Enable verbose",metavar="logical"),
  optparse::make_option(c("--vcf_object_name"), type="character", default="vcf_list",
                        help="Save file name for the read in data",metavar="logical")
)
opt_parser = optparse::OptionParser(option_list=option_list)
opt= optparse::parse_args(opt_parser)


ReadVcf<-function(vcf_file=opt$vcf_file, nrows=opt$nrows,skip=opt$skip,sample_nr=opt$sample_nr,verbose=opt$verbose){
if(file.exists(vcf_file)==F){return(message("No input vcf file has been submitted or the submitted file does not exist in the working directory."))
}

  ## Read in vcfR object
  # Collect Stats
  stats<-vcfR::.vcf_stats_gz(vcf_file,nrows=nrows,skip=skip,verbose=verbose)
  # Inventorise how many columns have to be read in
  if(sample_nr==0){
    cols<-sort(unique((1:stats[4])))
  }else{
    cols<-sort(unique(c(1:9,sample_nr)))}
  # Read the data
  message("Reading in data...")
  body<-vcfR::.read_body_gz(vcf_file,stats=stats,nrows=nrows,skip=skip,cols=cols,convertNA=as.numeric(1),verbose=as.numeric(verbose))
  # Select relevant contigs
  if(contig!="No Default"){
    c<-grep(contig,body[,1])
    c<-c[length(c)]
    body<-body[c[1]:c[length(c)],]
  }

  ## Order the vcf data into a matrix
  # Labelling variant type
  a<-cbind(which(nchar(body[,4])==1),"SNP")
  b<-cbind(which(nchar(body[,4])>1),"MNP")
  if(ncol(a)==1){
    # No SNP
    c<-as.data.frame(b)
  }else if(ncol(b)==1){
    c<-as.data.frame(a)
    # No MNP
  } else{
    c<-as.data.frame(rbind(a,b));c[,1]<-as.numeric(as.character(c[,1]));c[,2]<-as.character(c[,2])}

  Type<-c[order(c[,1]),2]
  # Info column sorting
  message("Sorting Info columns")
  info_sort<-list()
  if(length(which(grepl("END=",body[,8])!=T))>0){
    # Info column is not empty
    tX<-as.data.frame(body[,8]);tX[,1]<-as.character(tX[,1])
    ref<-body[which((stringr::str_count(body[,8], ";") + 1) == max(stringr::str_count(body[,8], ";") + 1))[1],8]
    ref<-unlist(strsplit(ref,";"))
    ref<-gsub("=.*","",ref)

    if(length(body[which((stringr::str_count(body[,8], ";") + 1) != max(stringr::str_count(body[,8], ";") + 1)),8])==0){
      # No missing parameters
      tX[,1]<-gsub(";",":",tX[,1])
    }else{
      # Some parameters are missing
      if(length(tX[which(stringr::str_count(body[,8], ref[1])<1),1])>0){
        # Replace first missing parameter
        tX[which(stringr::str_count(body[,8], ref[1])<1),1]<-sub("","NA;",tX[which(stringr::str_count(body[,8], ref[1])<1),1])}
      z<-2
      while(z<=length(ref)){
        # Replace non-first missing parameters
        tX[which(stringr::str_count(body[,8], ref[z])<1),1]<-sub(";",":NA;",tX[which(stringr::str_count(body[,8], ref[z])<1),1])
        # Log what parameter has been checked
        tX[which(stringr::str_count(body[,8], ref[z])>=1),1]<-sub(";",":",tX[which(stringr::str_count(body[,8], ref[z])>=1),1])
        z<-z+1
      }
    }
    # Separate
    tX<-tidyr::separate(tX, col=1, into=ref,sep=":")
    # General removal of parameter names
    z<-1
    while(z<=length(ref)){
      tX[,z]<-gsub(paste(ref[z],"=",sep=""),"",tX[,z])
      z<-z+1
    }
    tX[,length(ref)]<-gsub(";","",tX[,length(ref)])
    # Addition of new info column to the body
    info_sort[[1]]<-cbind(body[,1:3],Type,body[,4:7],tX)
    names(info_sort)[[1]]<-"General_Parameters"

  }else{info_sort<-cbind(body[,1:3],Type,body[,4:8])}

  # Sample Specific ordering:
  # Prepare an output
  message("Sorting Sample Specific columns")
  Sample_Specific_Parameters<-list()
  if(ncol(body)>9){
    # Get full parameter list
    ref<-body[which((stringr::str_count(body[,9], ":") + 1) == max(stringr::str_count(body[,9], ":") + 1))[1],9]
    ref<-unlist(strsplit(ref,":"))
    # Process each sample
    i<-1
    col<-10
    while(col<=ncol(body)){
      # Process sample specific data
      tX<-as.data.frame(body[,col]);tX$`body[, col]`<-as.character(tX$`body[, col]`)
      if(length(body[which((stringr::str_count(tX[,1], ":") + 1) != max(stringr::str_count(body[,9], ":") + 1))[1],9])==0){
        # No parameters are missing
        tX<-tidyr::separate(tX, col=1, into=ref,sep=":")
      }else{
        # Some parameters are missing
        if(length(tX[which(stringr::str_count(body[,9], ref[1])<1),1])>0){
          # Replace first missing parameter
          tX[which(stringr::str_count(body[,9], ref[1])<1),1]<-sub("","NA:",tX[which(stringr::str_count(body[,9], ref[1])<1),1])}
        z<-2
        while(z<=length(ref)){
          # Replace non-first missing parameters
          tX[which(stringr::str_count(body[,9], ref[z])<1),1]<-sub(":",";NA:",tX[which(stringr::str_count(body[,9], ref[z])<1),1])
          # Log what parameter has been checked
          tX[which(stringr::str_count(body[,9], ref[z])>=1),1]<-sub(":",";",tX[which(stringr::str_count(body[,9], ref[z])>=1),1])
          z<-z+1
        }
        # Replace the remainder of the parameter separators and split into columns
        tX[,1]<-gsub("(;\\.)",";NA",tX[,1])
        suppressWarnings(tX<-tidyr::separate(tX, col=1, into=ref,sep=";"))
      }

      # Store work and prepare for the next sample
      Sample_Specific_Parameters[[i]]<-tX
      names(Sample_Specific_Parameters)[i]<-colnames(body)[col]
      message("Processed sample ",i," out of ",ncol(body)[1]-9)
      i<-i+1
      col<-col+1}

    # End of samples loop
  }

  # Compile
  output<-c(info_sort,Sample_Specific_Parameters)
  # Change parameters into R numerics and characters
  message("Changing into R characters and numericals")
  z<-1
  while(z<=length(output)){
    i<-1
    while(i<=ncol(output[[z]])){
      if(length(output[[z]][which(grepl("[B|C|D|F|G|H|I|J|K|L|M|O|P|Q|R|S|T|U|V|W|X|Y|Z|:|,|//|//|]",unname(output[[z]][,i]),ignore.case = T)==T),i])>0){
        # The data contains non-numerics
        output[[z]][,i]<-as.character(output[[z]][,i])
        i<-i+1
      }else{
        # The data contains only numerics
        output[[z]][,i]<-suppressWarnings(as.numeric(as.character(output[[z]][,i])))
        i<-i+1
      }}
    z<-z+1
  }

  # End of function:
  return(output)

}
# Execute
x<-ReadVcf()
# Store
opt$vcf_object_name<-paste(opt$vcf_object_name,".RData",sep="")
save(x,file=opt$vcf_object_name)
