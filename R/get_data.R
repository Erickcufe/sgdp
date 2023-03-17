download_sgd_snps <- function(output_file) {
  # URL for SNP data
  url <- "https://storage.googleapis.com/simons-dogs-project-public/release/2022-02-11/sgdp_snps_filt.bed"

  # Check for HTTP errors
  resp <- httr::GET(url)
  if (resp$status_code != 200) {
    stop(paste0("Failed to download SNP data. HTTP status code: ", resp$status_code))
  }

  # Download SNP data and save to file
  download.file(url, destfile = output_file)

  # Convert binary BED file to tab-separated text file
  system(paste0("bedtools", " ", "convert2 -i", " ", output_file, " ", "-o", " ", "tsv"), intern = TRUE)
}




download_sgd_snps("Simons_Data.tsv")
