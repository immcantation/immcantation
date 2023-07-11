library(airr)
library(dplyr)

# Script to extract fasta file and the necessary filtered_contig_annotations.csv file from the BCR.data.rds file for the tutorials
tab <- readRDS("BCR.data.rds") %>% ungroup()

tab_lightchain <- tab %>% filter(locus %in% c("IGL", "IGK")) 
tab_heavychain <- tab %>% filter(locus %in% c("IGH")) %>%
                          mutate( c_gene = case_when( isotype == "IgM" ~ "IGHM",
                                                      isotype == "IgG" ~ "IGHG",
                                                      isotype == "IgA" ~ "IGHA",
                                                      isotype == "IgD" ~ "IGHD"))

tab_lightchain_fixed <- tab_lightchain %>%
                          mutate( sequence_id = paste(sequence_id, sample, sep = "_") ) 

tab_merged <- bind_rows(tab_heavychain, tab_lightchain_fixed)
write_rearrangement(tab_merged %>% select(-c(c_gene, clone_id, clone_count, mu_freq)), "BCR_data.tsv")

colnames(tab)

# Convert AIRR rearrangment format to filtered_contig_annotations.csv format for example data
tab_10x <- tab_merged %>% ungroup() %>%
            mutate(
            barcode = sequence_id,
            is_cell = "true",
            contig_id = sequence_id,
            sample_id = sample,
            high_confidence = "true",
            length = length(sequence),
            chain = locus,
            v_gene = v_call,
            d_gene = d_call,
            j_gene = j_call,
            full_length = "true",
            productive = "true",
            fwr1_nt = fwr1,
            cdr1_nt = cdr1,
            fwr2_nt = fwr2,
            cdr2_nt = cdr2,
            fwr3_nt = fwr3,
            cdr3_nt = cdr3,
            cdr3 = cdr3_aa,
            fwr4_nt = fwr4,
            reads = consensus_count,
            umis = duplicate_count
        ) %>%
        select( barcode, is_cell, contig_id, high_confidence, length, chain, 
                v_gene, d_gene, j_gene, c_gene, full_length, productive, fwr1_nt, cdr1_nt,
                fwr2_nt, cdr2_nt, fwr3_nt, cdr3, cdr3_nt, fwr4_nt, reads, umis)

write.csv(tab_10x, "filtered_contig_annotations.csv", row.names = FALSE, quote = FALSE)

# Then run ChangeO ConvertDB to get the fasta file
# ConvertDb.py fasta -d /data/assets/BCR_data.tsv -o /data/assets/BCR_data_sequences.fasta --if sequence_id --sf sequence