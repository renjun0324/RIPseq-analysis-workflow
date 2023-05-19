
#-----------------------------------------------------------------
#                                                                
#                          enrichment                     
#                                                                 
#-----------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)

load("chipseeker_annotatedPeak.rda")
# entrezid <- unique(df$ENTREZID[which(df$pvalAdj < 0.1 | df$eFDR < 0.1)])
entrezid <- unique(df$ENTREZID[which(df$peak_qvalue < 0.1)])
# entrezid <- unique(df$ENTREZID)
entrezid <- entrezid[!is.na(entrezid)]

### go enrich
y <- enrichGO(gene = entrezid,
              OrgDb = org.Hs.eg.db, 
              ont = "ALL",
              pAdjustMethod = "BH",
              minGSSize = 10,
              maxGSSize = 500,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.1,
              readable = TRUE) 
y <- y %>% data.frame
gr <- stringr::str_split_fixed(y$GeneRatio,"/",2)
br <- stringr::str_split_fixed(y$BgRatio,"/",2)
fe <- (as.numeric(gr[,1])/as.numeric(gr[,2])) / (as.numeric(br[,1])/as.numeric(br[,2]))
y$fold_enrichment <- fe
y$generatio <- as.numeric(gr[,1])/as.numeric(gr[,2])
y$bgratio <- as.numeric(br[,1])/as.numeric(br[,2])
# go <- y[which(y$Count>=10),]
go <- y

### kegg enrich
y <- enrichKEGG(entrezid,
                organism = 'human',
                keyType = 'kegg',
                pAdjustMethod = "BH",
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.1,
                use_internal_data = FALSE)
y <- as.data.frame(y)
# y = setReadable(y, OrgDb = org.Hs.eg.db, keyType="ENTREZID") %>% data.frame
gr <- stringr::str_split_fixed(y$GeneRatio,"/",2)
br <- stringr::str_split_fixed(y$BgRatio,"/",2)
fe <- (as.numeric(gr[,1])/as.numeric(gr[,2])) / (as.numeric(br[,1])/as.numeric(br[,2]))
y$fold_enrichment <- fe
y$generatio <- as.numeric(gr[,1])/as.numeric(gr[,2])
y$bgratio <- as.numeric(br[,1])/as.numeric(br[,2])
# kegg <- y[which(y$Count>=10),]
kegg <- y

### save
enrichment <- list(go = go, kegg = kegg)
save(enrichment, file = "chipseeker_enrichment.rda")
system("rm chipseeker_enrichment.xlsx")
for(i in c("go","kegg")){
  tryCatch(
    {
      xlsx::write.xlsx(enrichment[[i]],
                       file = paste0("chipseeker_enrichment.xlsx"),
                       sheetName = i,
                       col.names = T,
                       row.names = T,
                       append = T)
    },
    error = function(e){
      message("nothing")
    }
  )
}

