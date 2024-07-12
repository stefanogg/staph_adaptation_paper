# Column names of a gff file
library(tidyverse)
gff_colnames <- c("CHROM",
                  "SOURCE",
                  "TYPE",
                  "START",
                  "END",
                  "SCORE",
                  "STRAND",
                  "PHASE",
                  "ATTRIBUTES")

# Global function to read a gff file
read_gff <- function(gff_file, # character with path to gff file
                     n_skip = NULL, # integer of number of initial lines to skip (default: 0)
                     supp_colnames = NULL, # character vector with names of columns at positions 10+ (e.g. positions of interval that was intersected with bed file, sequence id, and so on)
                     parse_gff = TRUE){ # should read_gff parse the attributes columns of the raw gff?
  
  my_colnames <- c(gff_colnames,
                   supp_colnames)
  
  raw_gff <- read_tsv(file = gff_file, 
                      skip = n_skip,
                      col_names = my_colnames)
  
  if (parse_gff) {
    out <- gff_to_tab(raw_gff)
  } else {
    out <- raw_gff
  }
  
  return(out)
  
}


# Function to parse the attributes column of a gff file 
# It generates a dataframe with the attributes

parse_gff_attributes <- function(gff_attributes) { # a character with the attributes column
    
    features <- gff_attributes
    
    list_features <- list()
    
    for (i in 1:length(features)) {
      
      data <- str_split(features[i], pattern = ";")[[1]]
      
      message(glue::glue("Extracting features for {data[1]}"))
      
      name <- c()
      value <- c()
      
      for (j in 1:length(data)) {
        name[j] <- str_split_fixed(data[j], pattern = "=", n = 2)[, 1]
        value[j] <- str_split_fixed(data[j], pattern = "=", n = 2)[, 2]
      }
      
      m <- matrix(value, nrow = 1)
      colnames(m) <- name
      df <- as_tibble(m)
      
      list_features[[i]] <- df
      
    }
    
    df_features <- bind_rows(list_features)
    
    return(df_features)
}

# Function to reformat a gff file in tabular format (this is a wrapper round parse_gff_attributes)

gff_to_tab <- function(gff_dataframe) { # dataframe with unformatted gff columns. The attributes column must be named "ATTRIBUTES"
  
  df <- gff_dataframe
  
  gff_attributes <- df$ATTRIBUTES
  
  df_features <- parse_gff_attributes(gff_attributes)
  
  df_gff <- bind_cols(df, df_features) %>%
    select(-ATTRIBUTES)
  
  return(df_gff)
}