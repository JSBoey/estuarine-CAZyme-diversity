# Environment ----
library(XML)
library(tidyverse)

# Functions ----
# Create HTML link based on CAZyme family input
makeAddress <- function(cazyme_family) {
  paste0("http://www.cazy.org/", cazyme_family, "_characterized.html")
}

# Are tables paginated?
# If yes, function will return vector of page links
# Otherwise, return an empty vector

makePageLinks <- function(html_link) {
  # Read HTML as text lines
  html_vector <- readLines(html_link)
  
  # Find pagination lines
  match_vector <- str_subset(html_vector, '\\?debut_FUNC')
  
  # If there are HTML lines indicating pagination, create list of links, otherwise return
  # HTML link
  if (is_empty(match_vector)) {
    result <- html_link
  } else {
    page_link <- unique(
      unlist(
        str_extract_all(match_vector, "\\?debut_FUNC=\\d+#pagination_FUNC")
      )
    )
    page_link_vector <- paste0(html_link, page_link)
    result <- append(html_link, page_link_vector)
  }
  
  return(result)
}

# Download and concatenate tables
downloadTables <- function(page_links) {
  # Download HTML and parse into table
  tbs <- map(page_links, ~ {
    html_lines <- readLines(.x)
    
    # Replace breaks '<br>' with spaces
    parsed_breaks <- str_replace_all(html_lines, "<br>", " ")
    
    # Parse downloaded HTML lines
    parsed_html <- readHTMLTable(parsed_breaks)[[2]]
  })
  
  # Bind rows
  if (length(tbs) > 1) {
    bind_rows(tbs)
  } else {
    tbs[[1]]
  }
}

# Re-format tables
formatTable <- function(tb) {
  headers <- c("protein_name", "ec_number", "reference", "organism", "genbank", "uniprot", "pdb3d")
  
  tb0 <- tb %>% 
    # Remove unnecessary rows
    filter(!str_detect(V1, "suivante|Protein Name|Top")) %>% 
    # Clean up cells
    mutate(
      across(everything(), ~ str_trim(.x))
    )
  
  # Append headers
  if (ncol(tb0) > 7) {
    colnames(tb0) <- append(headers, "subfamily")
  } else {
    colnames(tb0) <- headers
  }
  
  return(tb0)
}

# Append taxonomic domain as column
appendTaxaDomain <- function(formatted_table) {
  # Create vector of where to split
  split_indicator <- is.na(formatted_table[["organism"]])
  
  # Split table by taxonomic domain
  split_table <- split(formatted_table, cumsum(split_indicator))
  
  # Get taxonomic domains
  taxa_domain <- map_chr(split_table, ~ .x[1, 1])
  split_table <- set_names(split_table, taxa_domain)
  taxa <- unique(taxa_domain) %>% 
    set_names(., .)
  
  # Append taxonomic domain to each row
  domain_tables <- map(taxa, ~ {
    keep(split_table, names(split_table) == .x) %>% 
      bind_rows() %>% 
      filter(!is.na(organism))
  }) %>% 
    bind_rows(.id = "taxa_domain")
}

# Wrapper function ----
getCAZyTables <- function(cazyme_family) {
  
  if (!is.character(cazyme_family) || length(cazyme_family) > 1) {
    stop("Input must be a character vector of length 1.")
  }
  
  output <- makeAddress(cazyme_family) %>%
    makePageLinks() %>%
    downloadTables() %>%
    formatTable() %>%
    appendTaxaDomain()
  
  gc()
  
  output
}

# Example usage ----
# GH5 <- getCAZyTables("GH5")
