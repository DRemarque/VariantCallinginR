#' SaveVcfObject()
#'
#' This function allows the user to save a vcf R object created by readVcf() into an excel file where sample specific parameters are stored in separate sheets. Note that this saving method is only functional for small datasets and not whole genome data.
#' @param vcf_list The R object into which the vcf data was stored by readVcf.
#' @param Output_Base_Name The desired output base name for the created excel file. [default="VcfObject"]
#' @keywords save vcf xlsx excel
#' @export
#' @examples
#' # Simple Save
#' SaveVcfObject(vcf_list=vcf_object,Output_Base_Name="My_Variants")

SaveVcfObject<-function(vcf_list="No Default",Output_Base_Name="VcfObject"){
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
