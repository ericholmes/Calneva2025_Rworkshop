## Zooper module 1 of 2: data retrieval tutorial ----------------------------
## Written for the Cal-Neva AFS workshop 2025-05-30 in Lodi
## Author: Eric Holmes, contact: eric.holmes@water.ca.gov

## This is a data exploration script working with data from the Zooper R package created by folks at the IEP.
## More details on the Zooper package can be found here: https://github.com/InteragencyEcologicalProgram/zooper

# Load libraries ----------------------------------------------------------

##Run these lines (without the '#' sign) to install the zooper package. Only needs to be done once.
# install.packages("devtools")
# devtools::install_github("InteragencyEcologicalProgram/zooper")
# devtools::install_github("InteragencyEcologicalProgram/deltamapr")
library(zooper) # this package is not available from CRAN, if installing for the first time use the code above
library(deltamapr) # this package contains spatial data for the Delta. Not available from CRAN

# Define functions ----------------------------------------------------------

## Function to download delta outflow data from cdec.water.ca.gov
downloadCDEC <- function(site_no, parameter_code, interval, start_date, end_date){
  temp <- read.table(paste("http://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", site_no, "&SensorNums=",
                           parameter_code, "&dur_code=", interval, "&Start=", start_date, "&End=", end_date, sep=""),
                     header=FALSE, sep=",", skip=1, stringsAsFactors = F)
  
  temp <- temp[,c(5,7)]
  colnames(temp) <- c("datetime", "param_val")
  temp$site_no <- site_no
  temp$parameter_code <- parameter_code
  temp$datetime <- as.POSIXct(temp$datetime, format = "%Y%m%d %H%M")
  
  return(temp[,c("site_no", "datetime", "parameter_code", "param_val")])
}

# Download and prepare data --------------------------------------------

### Download Delta outflow data ----------------------------------------
dto <- downloadCDEC(site_no = "DTO", parameter_code = 23, interval = "D", 
                    start_date = "1980-10-1", end_date = "2020-9-30")

##Convert flow to numeric
dto$param_val <- as.numeric(dto$param_val)

##Format date for future merge
dto$Date <- as.Date(dto$datetime)

##add month, year, and wy fields to dto
dto$month <- as.integer(format(dto$datetime, format = "%m"))
dto$year <- as.integer(format(dto$datetime, format = "%Y"))
dto$wy <- ifelse(dto$month %in% c(10:12), dto$year + 1, dto$year)

### Download CA water year type data from CDEC ------------------------------
## get vector of column names
wycolnames <- trimws(read.fwf("https://cdec.water.ca.gov/reportapp/javareports?name=WSIHIST", 
                              skip = 1192, nrow = 1, header = F,
                              widths = c(4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8)))
##adjust column names
wycolnames[2:6] <- paste0("sac_", wycolnames[2:6])
wycolnames[7:11] <- paste0("sj_", wycolnames[7:11])

wytype <- read.fwf("https://cdec.water.ca.gov/reportapp/javareports?name=WSIHIST", 
                   skip = 1194, nrow = 1318-1194, col.names = wycolnames,
                   widths = c(4, 8, 8, 8, 8, 8, 8, 8, 8, 8,6))

wytype <- janitor::clean_names(wytype)

wytype$sac_yr_typefac <- factor(trimws(wytype$sac_yr_type), levels = c("W", "AN", "BN", "D", "C"))

### Download zooplankton data ----

##Function from Zooper package to download and synthesize zoop data from multiple surveys
##Sources:  "EMP" = Environmental Monitoring Program, "20mm" = 20 millimeter survey,
##          "FMWT" = Fall midwater trawl, "FRP" = Fish restoration program,
##          "DOP" = USBR Directed Outflow Project, "YBFMP" = Yolo Bypass Fish Monitoring Program

MyZoops <- Zoopsynther(Data_type = "Taxa", 
                       Sources = c("EMP", "20mm"), 
                       Size_class = "Meso", 
                       Date_range = c("2000-10-01", "2020-09-30"))

##filter out total counts
zoop <- MyZoops[!(grepl(MyZoops$Taxname, pattern = "all_Meso")),]

##Add julien day, year, month and water year
zoop$jday <- as.numeric(format(as.Date(zoop$Date), format = "%j"))
zoop$year <- as.integer(format(zoop$Date, format = "%Y"))
zoop$month <- as.integer(format(zoop$Date, format = "%m"))
zoop$wy <- ifelse(zoop$month %in% c(10:12), zoop$year + 1, zoop$year)

## Find lowest common taxonomic grouping
## Issue: UnID taxa limits grouping to the Class level
zoop$group <- ifelse(is.na(zoop$Class) == T, zoop$Phylum, zoop$Class)
zoop$group <- ifelse(zoop$group == "Eurotatoria", "Rotifera", 
                     ifelse(zoop$group %in% c("Insecta", "Annelida", "Ostracoda", "Malacostraca"), 
                            "other", zoop$group))

# Save data -----------------------------------------------------------------
## Save downloaded and prepared dataframes in an Rdata file
## Rdata files retain the structure of the data and are compressed to take up less disk space
## Added a conditional awaiting user input deciding whether to save data

if(readline(prompt = "Save data? (y/n):") == "y"){save(zoop, wytype, dto, 
     file = "data/Calneva2025_workshop_data.Rdata")}
