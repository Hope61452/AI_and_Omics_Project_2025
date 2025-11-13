
deg_results <- read.csv("C:/Users/USER/Desktop/AI_Omics_Internship_2025/preprocessing/Results/DEGs_Results.csv")

gc()
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("clusterProfiler", type = "binary")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("GO.db")
install.packages("msigdbr")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

# Convert gene symbols to Entrez IDS
symbols <- as.character(rownames(deg_results))
symbols <- deg_results$X
entrezid <- mapIds(x = org.Hs.eg.db,
                   keys = symbols,
                   column = "ENTREZID",
                   keytype = "SYMBOL",
                   multiVals = "first") %>%
stack() %>%
  dplyr::rename(Entrezid = values, Symbol = ind)

#Add Entrez IDS back to DEG reslts
deg_results <- deg_results %>%
  +   tibble::rownames_to_column("Symbol") %>%
  +   left_join(entrezid, by = "Symbol") %>% 
  +   relocate(Entrezid, .before = logFC )

deg_results <- deg_results %>%
  drop_na(Entrezid, adj.P.Val) %>%
  arrange(adj.P.Val) %>%
  distinct(Entrezid, .keep_all = TRUE)

#Subset Gene Groups
upregulated <- deg_results %>% filter(threshold == "Upregulated")
downregulated <- deg_results %>% filter(threshold == "Downregulated")
deg_updown <- bind_rows(upregulated, downregulated)

# Background and foreground sets
background_genes <- deg_results$Entrezid
sig_genes <- deg_updown$Entrezid


# GO Overrepresentation Analysis (ORA)#
## Biological Process ont = "BP"
## Molecular Functions = "MF"
## Cellular Components = "  CC"
## ont = "ALL"

go_ora <-enrichGO(
  gene = as.character(sig_genes),
  universe = as.character(background_genes),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  readable = TRUE
)

go_ora
GO <- as.data.frame(go_ora)

# visualization
dotplot(go_ora, showCategory = 20, title = "GO Overrepresentation (BP, MF, CC)")
barplot(go_ora, showCategory = 20, title = "GO Overrepresentation")

# Split Ontology

barplot(go_ora,
        split = "ONTOLOGY",
        showCategory = 5,
        font.size = 10,
        color = "qvalue",
        label_format = 60) +
  facet_grid(ONTOLOGY ~ ., scales = "free")

#KEGG Overrepresentation Analysis (ORA)
kegg_ora <- enrichKEGG(
  gene = as.character(sig_genes),
  organism = "hsa",
  keyType = "kegg",
  universe = as.character(background_genes)
)

kegg <- as.data.frame(kegg_ora)
dotplot(kegg_ora, title = "KEGG Pathway Enrichment")
upsetplot(kegg_ora) 

#Gene Set Enrichment Analysis (GSEA)
gene_ranks <- deg_results$logFC
names(gene_ranks) <- deg_results$Entrezid
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

#Retriee Hallmark gene sets from MSigDB
hallmark <- msigdbr(species = "Homo sapiens", category = "H")

hallmark <- hallmark %>%
  dplyr::select(gs_name, entrez_gene)


#Run GSEA
gsea_res <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = hallmark,
  eps =0
)
gsea <-as.data.frame(gsea_res)

# Visualization of enriched pathways

gseaplot(gsea_res, geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

p1 <- gseaplot(gsea_res,
               geneSetID = "HALLMARK_GLYCOLYSIS",
               by = "runningScore",
               title = "GLYCOLYSIS")

p2 <- gseaplot(gsea_res,
               geneSetID = "HALLMARK_INFLAMMATORY_RESPONSE",
               by = "runningScore",
               title = "INFLAMMATORY RESPONSE")

cowplot::plot_grid(p1, p2, ncol = 1)

#Display multiple pathways
gseaplot2(gsea_res, geneSetID = 1:3)

# display p value table
gseaplot2(gsea_res, geneSetID = 1:3,
          pvalue_table = TRUE)

# Summary and Reporting#

# Convert enrichment results to data frames
GO_df <- as.data.frame(go_ora)
KEGG_df <- as.data.frame(kegg_ora)
GSEA_df <- as.data.frame(gsea_res@result)

# Filter sigificant results
# Gene Ontology (GO)
GO_sig <- GO_df %>% filter(p.adjust < 0.05)

# separate significant GO results into BP, MF, and CC categories
top_GO_BP <- GO_sig %>% 
  filter(ONTOLOGY == "BP") %>% 
  arrange(p.adjust) %>% 
  select(ID, Description, p.adjust) %>% 
  head(10)


top_GO_MF <- GO_sig %>% 
  filter(ONTOLOGY == "MF") %>% 
  arrange(p.adjust) %>% 
  select(ID, Description, p.adjust) %>% 
  head(10)

top_GO_CC <- GO_sig %>% 
  filter(ONTOLOGY == "CC") %>% 
  arrange(p.adjust) %>% 
  select(ID, Description, p.adjust) %>% 
  head(10)

# Combine for summary Table
top_GO <- bind_rows(
top_GO_BP %>% mutate(Category = "Biological Process"),
top_GO_MF %>% mutate(Category = "Molecular Function"),
top_GO_CC %>% mutate(Category = "Celluar Component")
)

# KEGG
KEGG_sig <- KEGG_df %>% filter(p.adjust <0.05)

# GSEA
GSEA_sig <- GSEA_df %>% filter(p.adjust <0.05)

# Separate unregulated and downregulated GSEA sets

GSEA_up <- GSEA_sig %>% filter(NES > 0) %>% arrange(desc(NES))
GSEA_down <- GSEA_sig %>% filter(NES < 0) %>% arrange(NES)

# Extract gene symbols

GSEA_sig$core_genes <- sapply(GSEA_sig$core_enrichment, function(x){
  ids <- unlist(strsplit(x, "/"))
  symbols <- mapIds(org.Hs.eg.db,
                    keys = ids,
                    column = "SYMBOL",
                    keytype = "ENTREZID",
                    multiVals = "first")
  paste(unique(symbols[!is.na(symbols)]), collapse = ",")
})


#Extract top pathways
top_GO <- GO_sig %>% arrange(p.adjust) %>% 
  head(10) %>% 
  select(ID, Description, p.adjust)

top_KEGG <- KEGG_sig %>% 
  arrange(p.adjust) %>% 
  head(10) %>% 
  select(ID, Description, p.adjust)

top_GSEA_up <- GSEA_up %>%
  head(10) %>% 
  select(ID, NES, p.adjust)
  
top_GSEA_down <- GSEA_down %>%
  head(10) %>% 
  select(ID, NES, p.adjust)

# SUMMARY TABLE

summary_overall <- data.frame(
  Analysis = c("GO(ORA)", "KEGG(ORA)", "GSEA"),
  Total_Significant = c(nrow(GO_sig), nrow(KEGG_sig), nrow(GSEA_sig)),
  Upregulated = c(NA, NA, nrow(GSEA_up)),
  Downregulated = c(NA, NA, nrow(GSEA_down))
)

#View overview summary
summary_overall

#Export Results
