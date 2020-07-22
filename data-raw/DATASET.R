library(rmongodb)
library(mongolite)
library(future)
library(dplyr)

# Cargamos la Coleccion HGDP de MongoDB
HGDP_collection <- mongo(collection = "human", db = "HGDP")


# CARGANDO FUNCIONES
search_mongo <- function(IDs){
  a <- c('"V1"', IDs)
  a <- paste0(a, " : 1", collapse = ",")
  search <- paste0('{', a, '}')
  df <- HGDP_collection$find('{}', fields = search)
  return(df)

}

primer_limpia <- function(x){
  c <- list()
  for(i in 1:ncol(x)){
    vect <- paste(x[,i], collapse = "")
    vect <- stringr::str_remove_all(vect, "-")
    names(vect) <- colnames(x)[i]
    c[[i]] <- vect

  }

  return(c)
}

segunda_limpia <- function(letters_n, total_pop){

  while(letters_n != ""){
    t <- stringr::str_count(letters_n, "T")
    c <- stringr::str_count(letters_n, "C")
    g <- stringr::str_count(letters_n, "G")
    a <- stringr::str_count(letters_n, "A")

    all_nucl <- c(t,c,g,a)
    names(all_nucl) <- c("T", "C", "G", "A")
    all_nucl <- all_nucl[all_nucl > 0]

    if(length(all_nucl) == 1 ){
      freq_ancestral <- (all_nucl[which.max(all_nucl)] / 2) / total_pop
      freq_minor <- 0
      ancestral <- names(freq_ancestral)
      minor <- names(freq_ancestral)
      freq_ancestral <- unname(freq_ancestral)
      # freq_minor <- unname(freq_minor)
      data_freqs <- data.frame(SNP = names(letters_n),
                               Ancestral_freq = freq_ancestral,
                               MAF = 0,
                               Ancestral_allel = ancestral,
                               Minor_allel = minor)

    } else {

      freq_ancestral <- (all_nucl[which.max(all_nucl)] / 2) / total_pop
      freq_minor <- (all_nucl[which.min(all_nucl)] / 2) / total_pop
      ancestral <- names(freq_ancestral)
      minor <- names(freq_minor)
      freq_ancestral <- unname(freq_ancestral)
      freq_minor <- unname(freq_minor)
      data_freqs <- data.frame(SNP = names(letters_n),
                               Ancestral_freq = freq_ancestral,
                               MAF = freq_minor,
                               Ancestral_allel = ancestral,
                               Minor_allel = minor)
    }
    return(data_freqs)

  }
}

tercer_limpia <- function(contents_request){

  contents_request[sapply(contents_request, is.null)] <- NULL
  prueba4 <- rlist::list.stack(contents_request)
  return(prueba4)

}

limpieza_total <- function(df){
  total_pop <- ncol(df) -2
  prueba1 <- t(df[,-1])
  colnames(prueba1) <- (prueba1[1,])
  prueba1 <- prueba1[-1,]
  prueba2 <- primer_limpia(prueba1)
  names(prueba2) <- colnames(prueba1)

  future::plan(multiprocess)
  content <- furrr::future_map(prueba2, purrr::safely(segunda_limpia),
                               total_pop = total_pop,
                               .progress = TRUE)
  contents_1 <- purrr::transpose(content)
  contents_request <- contents_1[["result"]]
  prueba3 <- tercer_limpia(contents_request)

  return(prueba3)
}

# Metadata
ver <- read.delim("metadata_project.txt")




# Adygei in Caucasus, Russia SGDP

IDs_adygei <- c('"HGDP01401"', '"HGDP01402"')


Adygei_russia <- search_mongo(IDs_adygei) %>%
  limpieza_total()

usethis::use_data(Adygei_russia, overwrite = TRUE)


# Albanian in Albania SGDP

zapotec <- c('"zapo0099"', '"zapo0098"')
Zapotec_Mexico <- search_mongo(zapotec) %>%
  limpieza_total()

usethis::use_data(Albanian_Albania, overwrite = TRUE)

