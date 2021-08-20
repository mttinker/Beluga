rm(list = ls())
require(svDialogs)
require(rJava)
require(rChoiceDialogs)
require(rmarkdown)

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}
rspnse = dlg_message(c("This script is used to generate an html report that summarizes results of  ",
                       "an age-structured population model fit to St. Lawrence Beluga survey data. ",
                       "You can begin by selecting a Results file. ",
                       "Continue?"), "okcancel")$res
if(rspnse == "cancel"){
  stop_quietly()
}
Sys.sleep(.1)
file_list = list.files(path = "./results",pattern="beluga_base_rslt.*\\.rdata$",full.names = F)
rslt_list = grep("Report", file_list, value = TRUE, invert = TRUE)
rslt_list = rslt_list[order(rslt_list,decreasing = T)]
rdata_file = rselect.list(rslt_list, preselect = NULL, multiple = FALSE,
                          title = "Select results file" ,
                          graphics = getOption("menu.graphics")) 
if(length(rdata_file)==0){
  dlg_message(c("No data file selected"), "ok")
  stop_quietly()
}else{
  file.copy(paste0("./results/",rdata_file),"Results.rdata",overwrite = TRUE)
  resultsfilename = substr(rdata_file,1,nchar(rdata_file)-6)
  fitobj_file = paste0(resultsfilename,"_fit.RDS")
  file.copy(paste0("./results/",fitobj_file),"Results_fit.RDS",overwrite = TRUE)
}
YearT = substr(resultsfilename,18,21)
vers = read.csv("./data/Version.csv"); vers = vers$Version
title = paste0("St. Lawrence Beluga Population Model, ", YearT)
subtitle = paste0("Model ", vers,", Results file: ", resultsfilename)
Daterun = Sys.Date()
render("./Beluga_model_summary.Rmd",
       output_dir = "./results",
       output_file = paste0("Report_",resultsfilename,".html"),
       params = list(rep_title = title, rep_subtitle = subtitle, rep_date = Daterun)) # 
dlg_message(c("The results can be viewed by opening the approproiate 'Report' html file in the results folder"),
            "ok")
