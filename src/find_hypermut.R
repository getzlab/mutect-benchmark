cancer_dirs<- list.dirs("/xchip/cga_home/jhess/for_Qing/results2/", recursive = F, full.names = T)
df <- list()
for (i in seq_along(cancer_dirs)){
    cancer <- cancer_dirs[i]
    ppath <- paste(cancer, "patient_counts_and_rates.txt", sep = "/")
    if(file.exists(ppath)){
        
        patients.counts_and_rates <- read.delim(ppath)
        patients.counts_and_rates$ttype <- basename(cancer)
        nmuts <- patients.counts_and_rates$nmut
        thres <- median(nmuts)*100
        df[[i]] <- patients.counts_and_rates[patients.counts_and_rates$nmut > thres,]
        print(paste0(nrow(df[[i]]), "hypermutation patients identified"))
    }
    else (print(paste("cancer", ": patient count file not found")))
}
df_combined <- do.call(rbind, df)
write.table(df_combined,file = "hypermut_patients.txt", sep = "\t")
print("finished!")