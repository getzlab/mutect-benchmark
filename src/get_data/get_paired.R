# with dalmatian generated paths, get the uuid from GDC api with the help of
# TCGAbiolinks r package.


library(dplyr)
library(tidyr)
library(TCGAbiolinks)

maxchar <- 50

mpairs = read_csv("mpairs_update.csv")

get_uuid <- function(cohort_value){
    
    # project id for GDC
    project_id <- paste("TCGA", cohort_value, sep = "-")
    
    # subset mpairs to one cohort
    dalmatian <- mpairs %>% filter(cohort == cohort_value) %>% 
        select(patient, file_name, sample_type) %>% 
        mutate(file_name = gsub("TCGA_MC3.", "", file_name))
    
    
    # query GDC
    suppressMessages({
        q <- GDCquery(project = project_id,
                      data.category = "Raw sequencing data",
                      legacy = TRUE,
                      data.format = "BAM",
                      #platform = c("Illumina GA", "Illumina HiSeq"),
                      experimental.strategy = "WXS",
                      access = "controlled", 
                      data.type = "Aligned reads")
        
        manifest <- getResults(q,  cols = c("platform", 
                                           "sample_type", 
                                           "project", 
                                           "file_name",
                                           "sample.submitter_id",
                                           "id")) 
        
        manifest_processed <- manifest %>% rename(uuid = id) %>% 
            mutate( patient = substring(sample.submitter_id, 6, 12)) %>% 
            mutate( cohort = cohort_value) %>% 
            mutate( sample_type = ifelse( grepl("normal", tolower(sample_type)), "normal", "tumor")) %>% 
            mutate( sample_name_manifest = substr(file_name, nchar(file_name) -maxchar, nchar(file_name)))
        tt <- left_join(dalmatian, manifest_processed, by = c("patient", "sample_type")) %>% 
            
            select(cohort, patient, uuid, sample_type, sample_name_manifest) 
    })
    
    uuid_completion_rate <- mean(!is.na(tt$uuid))
    print(sprintf("the non-missing rate for UUID is %.03f", uuid_completion_rate))
    print(table(is.na(tt$uuid),tt$sample_type))
    return(tt)
}

if (FALSE){
    res <- get_uuid("PAAD")
}


##### batch request uuid ######
i <- 0
res <- data_frame()
for (cohort in unique(mpairs$cohort)){
    i <- i+1
    print(sprintf("current cohort requested: %s", cohort))
    temp <- get_uuid(cohort)
    res <- rbind(res, temp)
    if (i %% 10 == 0) gc()
}


#### see if the mpairs is also duplicated ####
mtab <- mpairs %>% 
    select(cohort, patient, sample_type, file_name) %>% 
    pivot_wider( values_from = file_name, names_from =  sample_type, values_fn = list(uuid = length)) %>% 
    filter(!is.na(tumor) & !is.na(normal)) %>% 
    mutate(tumor = substr(tumor, nchar(tumor)-maxchar, nchar(tumor)),
           normal = substr(normal, nchar(normal)-maxchar, nchar(normal))) 

#### correct the duplicated entries by firecloud bam #####
unique_df <- res %>% select(-sample_name_manifest) %>% 
    pivot_wider(values_from = uuid, names_from =  sample_type, values_fn = list(uuid = length)) %>% 
    filter(tumor == 1 & normal == 1) %>% select(cohort, patient)  %>% 
    left_join(res) %>% select(cohort, patient, sample_type, uuid) %>% 
    pivot_wider(values_from = uuid, names_from =  sample_type)

# identify the duplicated ones
res_dups <- res %>% select(-sample_name_manifest) %>% 
    pivot_wider(values_from = uuid, names_from =  sample_type, values_fn = list(uuid = length)) %>% 
    filter(tumor != 1 | normal != 1) %>% select(cohort, patient) 

correct_df <- left_join(res_dups, mtab) %>% 
    pivot_longer(cols = c(tumor, normal), names_to = "sample_type", values_to = "sample_name_manifest") %>% 
    inner_join(res, by = c("cohort", "patient", "sample_type", "sample_name_manifest")) %>%
    select(cohort, patient, sample_type, uuid) %>% 
    pivot_wider(values_from = uuid, names_from =  sample_type) %>% 
    rbind(unique_df)
    


write_tsv(correct_df, "paired_uuid.tsv")

