#################################################
####     Database Variables to Lowercase     ####
####                                         ####
####            By: Zach Laubach             ####
####        last updated: 20 Dec 2017        ####
#################################################

### PURPOSE: This function identifies variables that are class 'character'
### and converts the text data to lower case. 

### ERRORS / BUG REPORTS: Please send to zchlaubach@gmail.com

  ## AllCharactersToLower()---- Function to check all variables in a dataframe 
  ## that are class 'character'and change the text data in these variables 
  ## to lowercase 
    AllCharactersToLower <- function (df) {
          # apply function to each column of dataframe
      data.frame (lapply (df, function (variable) {
          # Set the local system encoding for mac as inherited from rApache.
          # This handles multibyte error.
        Sys.setlocale('LC_ALL', locale ='C') 
          # if variable is class character, change text to lowercase
        if (is.character (variable)) return (tolower (variable)) 
        else return (variable)
      }))
    }