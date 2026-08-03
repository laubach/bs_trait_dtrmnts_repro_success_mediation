#################################################
####     Database Variables to Lowercase     ####
####                                         ####
####            By: Zach Laubach             ####
####        last updated: 18 Jan 2018        ####
#################################################

### PURPOSE: This function formats all variables to match R style guide
### recommendations. 

### ERRORS / BUG REPORTS: Please send to zchlaubach@gmail.com

  ## FormatVarNames()---- Function to convert all camelcase variable names
  ## in a datafram to the recommended R style guide format for variables
  ## for example 'NameVarible' to name.variable
    FormatVarNames <- function (df) {
      #add a "." between lower case letters and uppercase letters
      colnames(df) <- gsub("([a-z])([A-Z])", "\\1\\.\\2", colnames(df)) 
      #add a "." between where there is a "_"
      colnames(df) <- gsub("\\_", "\\1\\.\\2", colnames(df))
      #add a "." between where there is a " "
      colnames(df) <- gsub("\\ ", "\\1\\.\\2", colnames(df))
      # convert all letters to lowercase
      colnames(df) <- tolower(colnames(df))
      # return the data frame with updated varibable names
      df = as_tibble(df)
      
    }

