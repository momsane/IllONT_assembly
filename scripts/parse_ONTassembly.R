library(tidyverse)

bases_tabs<-unlist(snakemake@input[["map_bases"]])

bases_concat<- function(bases_file){
  sam_name<- unlist(strsplit(bases_file, split = "/"))
  sam_name<-unlist(strsplit(sam_name[length(sam_name)], split = "_"))[1]

  df<- read.csv(bases_file, header=F, sep="\t") %>% mutate ("sample"=sam_name)
  colnames(df)<- c("contig", "base", "cov", "sample")
  head(df)
  return(df)
}

a<- do.call("rbind", lapply(bases_tabs, bases_concat))

write.table(a, file=snakemake@output[["all_bases"]], sep = "\t", quote=F, , row.names = F)

#

flagstats<-unlist(snakemake@input[["map_reads"]])

format_flagstat<- function(file){
  sam_name<- unlist(strsplit(file, split = "/"))
  sam_name<-unlist(strsplit(sam_name[length(sam_name)], split = "_"))[1]

  df<-read.csv(file, header=F, sep=" ") %>% drop_na() %>%
    dplyr::summarize("data"=as.numeric(V1),
                     "header"=c("Total_reads","primary", "secondary", "supplementary", "duplicates", "primary_duplicates",
                                                   "mapped", "primary_mapped", "paired_in_seq", "read1", "read2", "properly_paired", "with_itself_and_mate_mapped",
                                                   "singletons", "with_mate_mapped_to_a_different_chr", "mapQ" )) %>%
    mutate("sample"=sam_name) %>%
    pivot_wider(names_from = header, values_from = data)

  print(df)

  df<- df %>% mutate(perc_mapped=as.numeric(mapped)/as.numeric(Total_reads),
                   perc_prop=as.numeric(properly_paired)/as.numeric(Total_reads),
                   perc_singleton=as.numeric(singletons)/as.numeric(Total_reads))



  return(df)
}

b<- do.call("rbind", lapply(flagstats, format_flagstat))
write.table(b, file=snakemake@output[["all_reads"]], sep = "\t", quote=F, row.names = F)
