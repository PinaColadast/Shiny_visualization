if(Sys.info()[[1]] == "Windows"){
  #===================================================================================================
  # THis is running on my local windows
  #===================================================================================================
  #identify where is pandoc in Rstrudio
  # Sys.getenv("RSTUDIO_PANDOC")
  #then set the position of R pandoc -> Sys.setenv(RSTUDIO_PANDOC="--Location--")
  Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/bin/pandoc")
  
  #Define the App directory
  APP_DIRECTORY <- getwd()
  
  #set the working folder
  setwd(APP_DIRECTORY)
  
  #define a temp directory
  TEMP_DIRECTORY <- "H:/Temp/"
  
  #Calculate the number of core of the Machine
  NB_CORES       <- 1   # detectCores()
  
  #Reticulate for windows 
  #==========================================================================
  # Sys.setenv(PATH= paste("C:/Users/sfontaine/AppData/Local/anaconda/Library/bin",
  #                        Sys.getenv()["PATH"], 
  #                        sep = ";"))
  # Sys.setenv(RETICULATE_PYTHON = "C:/Users/sfontaine/AppData/Local/anaconda/python.exe")
  # 
  # 
  # library(reticulate)



}else if(Sys.info()[[1]] == "Linux"){
  #===================================================================================================
  # THis is running on the shiny app server
  #===================================================================================================
  
  #Calculate the number of core of the Machine
  NB_CORES       <- detectCores()/3
  
  Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio/bin/pandoc")
  
  #Define the App directory
  APP_DIRECTORY <- getwd()
  
  #set the working folder
  setwd(APP_DIRECTORY)
  
  #define a temp directory
  TEMP_DIRECTORY <- "/var/tmp/"
  
  #Calculate the number of core of the Machine
  NB_CORES <- detectCores()/3
  
  library(reticulate)
  use_python("/home/administrator/anaconda3/bin/python", required = T)
  
}