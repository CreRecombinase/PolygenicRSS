library(ldmap)
library(dplyr)
library(purrr)
library(readr)
library(vroom)
input_f <- snakemake@input[["bimf"]]
stopifnot(file.exists(input_f),!is.null(input_f))
output_f <- snakemake@output[["snplistf"]]
host <- snakemake@params[["host"]] %||% "gardner"
has_chr <- as.logical(snakemake@params[["has_chr"]] %||% FALSE)

col_chromosome <- function (prefix_chr = TRUE, ...) 
{
  if (prefix_chr) {
    return(vroom::col_factor(levels = ldmap:::chromosome_levels(TRUE)))
  }
  return(vroom::col_factor(levels = ldmap:::chromosome_levels(FALSE)))
}


bim_cols <- function (chrom = col_chromosome(prefix_chr = TRUE), rsid = vroom::col_character(), 
                      map = vroom::col_double(), pos = vroom::col_integer(), alt = vroom::col_character(), 
                      ref = vroom::col_character(), ...) 
{
  return(vroom::cols(chrom = chrom, rsid = rsid, map = map, 
                     pos = pos, alt = alt, ref = ref))
}


if(host!="gardner"){

rpb <- function(file, compact = TRUE,has_chr=FALSE, cols = bim_cols(chrom=col_chromosome(has_chr)),read_fun=vroom::vroom) {
  ret_df <- dplyr::filter(vroom(file, col_names = names(cols$cols),col_types=cols,delim="\t"),nchar(alt)==1,nchar(ref)==1)
  if (compact)
    return(compact_snp_struct(ret_df))
  return(ret_df)
}



  walk2(input_f, output_f,
        function(input,output){
          rpb(input,has_chr=has_chr) %>%
            count(rsid) %>%
            filter(n==1) %>%
            pull(rsid) %>% 
            write_lines(output)
        })
}else{
  input_df <- readr::read_tsv(input_f,col_names=FALSE)
  ctp <- count(input_df,X3)
  sct <- semi_join(input_df,
                   filter(ctp,n==1))
  sct_snp <- sct %>% filter(nchar(X4)==1,
                            nchar(X5)==1)
  filter(sct_snp,
         X6>0.01) %>% pull(X2) %>% write_lines(output_f)
}
