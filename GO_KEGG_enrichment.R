
# KEGG enrichment analysis
suppressPackageStartupMessages({
options(stringsAsFactors = F) 
Sys.setenv(R_MAX_NUM_DLLS=999) 
library(Seurat)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(stringr) 
library(ggplot2)
})

setwd("pathway/to/analysis")
outdir="./"
pro="silkworm_deg_kegg"
rds="sample10_ASG_MPSG_umap_CellType.rds"
myColor = c("#e11159", "#b43e94", "#16b384", "#e8be6a", "#f3eae5")


## load seurat object and set idents
sc=readRDS(rds)
Idents(sc) = "sample"
sc$sample = factor(sc$sample, levels = c("D_L5D1","D_L5D3","D_L4M", "D_L5D5", "W_L5D5"))
names(myColr) = c("D_L5D1","D_L5D3","D_L4M", "D_L5D5", "W_L5D5")

## load silkbase information
gene_keggID_all = read.csv("./Silkbase_KEGGID_description_GO_KEGG_final.csv",header = T) 
head(gene_keggID_all)

## KEGG analysis
ident1="ident1"
ident2="ident2"

deg = FindMarkers(object = sc, 
                  ident.1 = ident1,
                  ident.2 = ident2,
                  logfc.threshold = 0.15,
                  est.use='wilcox' )
deg = deg[order(deg$avg_log2FC, decreasing = T),]
deg = deg %>% mutate(Silkbase_id = row.names(deg))

deg = merge(deg, gene_keggID_all, by.x = 'Silkbase_id', by.y = 'Silkbase_id')  
row.names(deg) = deg$Silkbase_id
head(deg)

## set avg_log2FC greater than 0.5 as significantly up and down gene sets
deg_up = deg[deg$avg_log2FC>0.5,]
gene_up= rownames(deg_up)
length(gene_up)
head(gene_up)

deg_down = deg[deg$avg_log2FC < -0.5,]
gene_down= rownames(deg_down)
head(gene_down)
length(gene_down)


write.csv(deg, file = paste0(outdir, pro, ident1, "_VS_", ident2, "_DEG_description.csv"))
write.csv(deg_up, file = paste0(outdir, pro, ident1, "_VS_", ident2,"_DEG_description_up.csv"))
write.csv(deg_down, file = paste0(outdir, pro, ident1, "_VS_", ident2,"_DEG_description_down.csv"))

write.csv(gene_up, file =   paste0(outdir, pro, ident1, "_VS_", ident2, "_DEG_up_geneid.csv"))
write.csv(gene_down, file = paste0(outdir, pro, ident1, "_VS_", ident2, "_DEG_down_geneid.csv"))

#####################################################################
# KEGG analysis
## Up genes KEGG analysis

# deg_up = read.csv(paste0(ident1, "_VS_", ident2,  "_DEG_up_geneid.csv"), header = T, row.names = 1) 
colnames(deg_up) = "GENE"
head(deg_up)

deg_up_keggID = gene_keggID_all$KEGG_id[match(deg_up$GENE, gene_keggID_all$Silkbase_id)] 
deg_up_keggID = na.omit(deg_up_keggID) # delete NA
length(deg_up_keggID)

kk.up <-enrichKEGG(deg_up_keggID, 
           organism = "bmor",  ## http://www.genome.jp/kegg/catalog/org_list.html. bmor	KGB	Bombyx mori (domestic silkworm)
           keyType = "kegg", #  “kegg”, ‘ncbi-geneid’, ‘ncib-proteinid’ and ‘uniprot’
           pAdjustMethod = "BH", 
           pvalueCutoff = 0.9,
           qvalueCutoff =0.9
           )
p_kegg_up = dotplot(kk.up)
p_kegg_up

write.csv(kk.up, file = paste0(ident1, "_VS_", ident2, "_kk.up.csv"))
ggsave(p_kegg_up, file= paste0(ident1, "_VS_", ident2, "_kk.up.pdf")) 

## Down genes KEGG analysis
# deg_down = read.csv(paste0(outdir, pro, ident1, "_VS_", ident2, "_DEG_down_geneid.csv"), header = T, row.names = 1) 
colnames(deg_down) = "GENE"
deg_down_keggID = gene_keggID_all$KEGG_id[match(deg_down$GENE, gene_keggID_all$Silkbase_id)] 
deg_down_keggID = na.omit(deg_down_keggID) # delete NA
length(deg_down_keggID)

kk.down <-enrichKEGG(deg_down_keggID, 
           organism = "bmor",  #  http://www.genome.jp/kegg/catalog/org_list.html. bmor	KGB	Bombyx mori (domestic silkworm)
           keyType = "kegg", #  “kegg”, ‘ncbi-geneid’, ‘ncib-proteinid’ and ‘uniprot’
           pAdjustMethod = "BH", 
           pvalueCutoff = 0.9,
           qvalueCutoff =0.9
                  )
p_kegg_down = dotplot(kk.down)
p_kegg_down

write.csv(kk.down, file = paste0(outdir, pro, ident1, "_VS_", ident2, "_kk.down.csv"))
ggsave(p_kegg_down, file= paste0(outdir, pro, ident1, "_VS_", ident2, "_kk.down.pdf")) 

## Up and down KEGG analysis
gene_diff = unique(c(deg_up_keggID,deg_down_keggID))
kk.diff <-enrichKEGG(gene_diff, 
           organism = "bmor",  # http://www.genome.jp/kegg/catalog/org_list.html. bmor	KGB	Bombyx mori (domestic silkworm)
           keyType = "kegg", #  “kegg”, ‘ncbi-geneid’, ‘ncib-proteinid’ and ‘uniprot’
           pAdjustMethod = "BH", 
           pvalueCutoff = 0.9,
           qvalueCutoff =0.9
                  )
 
head(kk.diff)[,1:6]
p_kegg_diff = dotplot(kk.diff)
p_kegg_diff

kegg_diff_dt = as.data.frame(kk.diff)
kegg_down_dt = as.data.frame(kk.down)
kegg_up_dt = as.data.frame(kk.up)
    
down_kegg = kegg_down_dt[kegg_down_dt$pvalue<0.01,];down_kegg$group=-1
up_kegg = kegg_up_dt[kegg_up_dt$pvalue<0.01,];up_kegg$group=1
    
kegg_plot <- function(up_kegg,down_kegg){
      dat=rbind(up_kegg,down_kegg)
      colnames(dat)
      dat$pvalue = -log10(dat$pvalue)
      dat$pvalue=dat$pvalue*dat$group 
      dat=dat[order(dat$pvalue,decreasing = F),]
      g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
        geom_bar(stat="identity") + 
        scale_fill_gradient(low="blue",high="red",guide = "none") + 
        scale_x_discrete(name ="Pathway names") +
        scale_y_continuous(name ="log10P-value") +
        coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
        ggtitle("Pathway Enrichment") 
    }
    
g_kegg=kegg_plot(up_kegg,down_kegg)
g_kegg
ggsave(g_kegg, file = paste0(outdir, pro, ident1, "_VS_", ident2, "_UP_down_bar.png"))
ggsave(g_kegg, file = paste0(outdir, pro, ident1, "_VS_", ident2,"_UP_down_bar.pdf")) 
      

#######################################################################
# GO analysis
## Up gene GO analysis
   go.up <-enrichGO(deg_up_keggID, 
           OrgDb = bmor,  
           pvalueCutoff = 0.05,
           qvalueCutoff =0.05,
           ont = "all"
    head(go.up)[,1:6]

    p_go_up_bar = barplot(go.up, showCategory=20)+    
                  scale_y_discrete(labels=function(x) str_wrap(x, width = 100))
    p_go_up_bar
    ggsave(p_go_up_bar, file= paste0(outdir, pro, ident1, "_VS_", ident2, "_GO.up_bar_top.png")) 
    ggsave(p_go_up_bar, file= paste0(outdir, pro, ident1, "_VS_", ident2, "_GO.up_bar_top.pdf"))                
                                   

p_go_up_bar_facet = barplot(go.up,
        split="ONTOLOGY")+
  facet_grid(ONTOLOGY~., scale="free")+
  scale_y_discrete(labels=function(x) str_wrap(x, width = 100))+
  theme(
    axis.text = element_text(size = 1, angle = 0, hjust = 1, vjust = 0.5),
  ) 
p_go_up_bar_facet
ggsave(p_go_up_bar_facet, file= paste0(outdir, pro, ident1, "_VS_", ident2, "_GO.up_bar_facet.png")) 
ggsave(p_go_up_bar_facet, file= paste0(outdir, pro, ident1, "_VS_", ident2, "_GO.up_bar_facet.pdf")) 

## Down gene GO analysis
go.down <-enrichGO(deg_down_keggID, 
           OrgDb = bmor,  
           pvalueCutoff = 0.05,
           qvalueCutoff =0.05,
           ont = "all"
                  )
head(go.down)[,1:6]

p_go_down_bar = barplot(go.down, showCategory=20)+    
                  scale_y_discrete(labels=function(x) str_wrap(x, width = 100))
p_go_down_bar
ggsave(p_go_down_bar, file= paste0(outdir, pro, ident1, "_VS_", ident2, "_GO.down_bar_top.png")) 
ggsave(p_go_down_bar, file= paste0(outdir, pro, ident1, "_VS_", ident2, "_GO.down_bar_top.pdf"))                
                                   
p_go_down_bar_facet = barplot(go.down,
        split="ONTOLOGY")+
  facet_grid(ONTOLOGY~., scale="free")+
  scale_y_discrete(labels=function(x) str_wrap(x, width = 100))+
  theme(
  axis.text = element_text(size = 1, angle = 0, hjust = 1, vjust = 0.5),
  ) 
p_go_down_bar_facet
ggsave(p_go_down_bar_facet, file= paste0(outdir, pro, ident1, "_VS_", ident2, "_GO.down_bar_facet.png")) 
ggsave(p_go_down_bar_facet, file= paste0(outdir, pro, ident1, "_VS_", ident2, "_GO.down_bar_facet.pdf")) 
   
