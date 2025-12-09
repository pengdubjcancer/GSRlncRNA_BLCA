##01热图
getwd() 
setwd("/home/pengdu/gem/01pRRophetic")
clin <- read.table("clinical.txt",sep = "\t",header =T,row.names = 1)
exp <- read.table("gem.txt",header = T,row.names = 1)
##样本排序 
identical(rownames(clin),rownames(exp))
exp <- as.matrix(exp)
exp <- t(exp)
annColors <- list() 
annColors[['Sex']] <- c('Male'='lightseagreen','Female'='plum')
annColors[['Subtype']] <- c('Non-Papillary'='paleturquoise','Papillary'='red','Not Available'='white')
annColors[['T']] <- c('T1'='cornsilk','T2'='paleturquoise','T3'='goldenrod','T4'='firebrick','TX'='white')
annColors[['N']] <- c('N0'='cornsilk','N1'='skyblue','N2'='steelblue','N3'='firebrick','NX'='white')
annColors[['M']] <- c('M0'='cornsilk','M1'='red','MX'='white')
annColors[['Stage']] <- c('I'='cornsilk','II'='lightseagreen',
                          'III'='blueviolet','IV'='red',
                          'Unknown'='white')
annColors[['Grade']] <- c('High'='red','Low'='cornsilk','Unknown'='white')
annColors[['Smoke']] <- c('No'='lightseagreen','Yes'='goldenrod')
annColors[['Lymphovascular_invasion']] <- c('No'='lightseagreen','Yes'='goldenrod','Unknown'='white')
annColors[['Status']] <- c('Alive'='steelblue','Dead'='goldenrod')
annColors[['Recurrence']] <- c('No'='steelblue','Yes'='goldenrod')
annColors[['Radiotherapy']] <- c('Yes'='yellowgreen','No'='pink')
annColors[['Chemotherapy']] <- c('Yes'='yellowgreen','No'='pink')

library(pheatmap)
#pdf(file="pheatmap.pdf",width = 10,height = 12)
p <- pheatmap(exp,
              color = colorRampPalette(c("#5bc0eb",'white',"#ECE700"))(100),
              annotation_col = clin,
              annotation_colors = annColors,
              show_rownames = T,
              show_colnames = F,
              cluster_rows = F,
              cluster_cols = F,
              cellheight = 15)
p
dev.off()
#结束出图 

##02degs
setwd("/home/pengdu/gem/02degs")
sample <- read.table("TcgaTargetGTEX_phenotype.txt",header = T,sep = "\t")
normal <- subset(sample,X_sample_type == "Solid Tissue Normal"|X_sample_type == "Normal Tissue")
tumor <- subset(sample,X_sample_type == "Primary Tumor")

rt = data.table::fread("TcgaTargetGtex_expected_count.gz", header=T, sep="\t", check.names=F)
rt <- as.data.frame.matrix(rt)
geneid <- rt[,1]
colnames(rt)
rt <- rt[,sample$sample]
rt <- cbind(geneid,rt)
# ID转换
library(dplyr)
probeMap <- read.delim("probeMap_gencode.v23.annotation.transcript.probemap")
rt1 <- rt %>%
  inner_join(probeMap, by = c("geneid" = "id")) %>%
  select(gene, starts_with("TCGA"))
rt2 <- rt %>%
  inner_join(probeMap, by = c("geneid" = "id")) %>%
  select(gene, starts_with("GTEX"))
rt <- cbind(rt1,rt2[,2:10])

library(limma)
rt <- as.data.frame(avereps(rt[, -1], ID = rt$gene)) # 取平均
dim(rt)
# 对于重复基因的处理还可以使用其他处理方式，如取表达值最大的一个
# 可以看到去除重复之后还剩58387个基因
# 读取并探索gtf文件
gtf <- rtracklayer::import("gencode.v23.annotation.gtf.gz")
gtf <- as.data.frame(gtf)
colnames(gtf)
gtf <- dplyr::select(gtf,c("gene_id","gene_name","gene_type"))
gtf <- unique(gtf)
library(tidyverse)
table(gtf$gene_type)
gtf_mrna <- gtf[gtf$gene_type=="protein_coding",] #提取mrna
lnc = c("3prime_overlapping_ncRNA", "antisense", "bidirectional_promoter_lncRNA", "lincRNA", "macro_lncRNA", "non_coding", "processed_transcript", "sense_intronic" , "sense_overlapping")
gtf_lncRNA <- gtf[gtf$gene_type %in% lnc,]
gtf_miRNA <- gtf[gtf$gene_type=="miRNA",] #提取mrna
#排序
sample <- rbind(tumor,normal)
rt <- rt[,sample$sample]
colnames(rt)
table(sample$X_sample_type)
rt <- as.data.frame(rt)
#lncrna矩阵
lncrna <- rt[gtf_lncRNA$gene_name,]
lncrna <- na.omit(lncrna)
write.table(lncrna,"lncrna.txt")
#mirna矩阵
mirna <- rt[gtf_miRNA$gene_name,]
mirna <- na.omit(mirna)
write.table(mirna,"mirna.txt")
#mrna矩阵
mrna <- rt[gtf_mrna$gene_name,]
mrna <- na.omit(mrna)
write.table(mrna,"mrna.txt")

group <- c(rep('Tumor',407),rep('Normal',28))
group <- factor(group)
table(group)

# lncrna构建DESeq2中的对象
collncrna <- data.frame(row.names = colnames(lncrna),
                          group = group)
collncrna$group <- factor(collncrna$group, levels = c("Tumor", "Normal"))
head(collncrna)
lncrna_int <- 2^(lncrna) - 1
lncrna_int <- apply(lncrna_int, 2, as.integer)
rownames(lncrna_int) <- rownames(lncrna)
library(DESeq2)
dds <-  DESeqDataSetFromMatrix(
  countData = lncrna_int,
  colData = collncrna,
  design = ~ group)
# 进行差异表达分析
dds <- DESeq(dds)
# 查看结果的名称
resultsNames(dds)
# contrast参数必须写成下面三个元素的向量格式，且顺序不能反
res <- results(dds, contrast = c("group", rev(levels(group))))
resOrdered <- res[order(res$padj), ]
DEG_DESeq2 <- as.data.frame(resOrdered)
write.table(DEG_DESeq2,file = "lncrna_DEG_DESeq2.txt",sep = "\t")
#lncrna火山图，设置筛选条件，自定义cutoff
input <- read.table("lncrna_DEG_DESeq2.txt",header = T,sep = "\t")
input$gene <- rownames(input)
input$Threshold=as.factor(ifelse(input$padj<0.05 & abs(input$log2FoldChange)>=log2(1.2),
                                 ifelse(input$log2FoldChange>log2(1.2),'Up','Down'),'No')) 
#查看log2FoldChange的最大值和最小值来确定下面xlim的参数
min(input$log2FoldChange)
max(input$log2FoldChange)

#参数说明
#geom_point(alpha=0.3, size=3) ：alpha调节点的透明度，size调节点的大小
#xlim(-10,10)：调节横坐标边界
#保存为pdf格式
#pdf(file = 'volcano.pdf',height = 5,width = 5)
##down top 10
down <- filter(input, Threshold == "Down") %>% 
  distinct(gene, .keep_all = T) %>%
  top_n(10, -log2FoldChange)
##up top 10
up <- filter(input, Threshold == "Up") %>% 
  distinct(gene, .keep_all = T) %>%
  top_n(10, log2FoldChange)

nudge_x_up = 7 - up$log2FoldChange
nudge_x_down = -6 - down$log2FoldChange
library(ggrepel)
ggplot(data = input,aes(x = log2FoldChange,y=-log10(padj),
                        colour=Threshold)) +geom_point(alpha=0.3, size=3) +
  scale_color_manual(values=c("blue", "grey","red"))+
  theme_bw()+xlim(-6,8)+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = NULL, colour = "black"))+
  geom_vline(xintercept=c(-log2(1.2),log2(1.2)),lty=2,col="black",lwd=0.4) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.4) +
  #'@添加关注的点的基因名
  #'@添加down top gene
  geom_text_repel(
    data = up,aes(x = log2FoldChange, y = -log10(padj), label = gene),
    seed = 123,color = 'black',show.legend = FALSE, 
    min.segment.length = 0,#始终为标签添加指引线段；若不想添加线段，则改为Inf
    segment.linetype = 1, #线段类型,1为实线,2-6为不同类型虚线
    segment.color = 'black', #线段颜色
    segment.alpha = 0.5, #线段不透明度
    nudge_x = nudge_x_up, #标签x轴起始位置调整
    direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
    hjust = 0, #对齐标签：0右对齐，1左对齐，0.5居中
    force = 2,#重叠标签间的排斥力
    force_pull = 2,#标签和数据点间的吸引力
    size = 4,
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    max.overlaps = Inf)+
  ##'@添加up top gene
  geom_text_repel(
    data = down,aes(x = log2FoldChange, y = -log10(padj), label = gene),
    seed = 123,color = 'black',show.legend = FALSE, 
    min.segment.length = 0,#始终为标签添加指引线段；若不想添加线段，则改为Inf
    segment.linetype = 1, #线段类型,1为实线,2-6为不同类型虚线
    segment.color = 'black', #线段颜色
    segment.alpha = 0.5, #线段不透明度
    nudge_x = nudge_x_down, #标签x轴起始位置调整
    direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
    hjust = 1, #对齐标签：0右对齐，1左对齐，0.5居中
    force = 2,#重叠标签间的排斥力
    force_pull = 2,#标签和数据点间的吸引力
    size = 4,
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    max.overlaps = Inf)

# mrna构建DESeq2中的对象
colmrna <- data.frame(row.names = colnames(mrna),
                      group = group)
colmrna$group <- factor(colmrna$group, levels = c("Tumor", "Normal"))
head(colmrna)
mrna_int <- 2^(mrna) - 1
mrna_int <- apply(mrna_int, 2, as.integer)
rownames(mrna_int) <- rownames(mrna)
library(DESeq2)
dds <-  DESeqDataSetFromMatrix(
  countData = mrna_int,
  colData = colmrna,
  design = ~ group)
# 进行差异表达分析
dds <- DESeq(dds)
# 查看结果的名称
resultsNames(dds)
# contrast参数必须写成下面三个元素的向量格式，且顺序不能反
res <- results(dds, contrast = c("group", rev(levels(group))))
resOrdered <- res[order(res$padj), ]
DEG_DESeq2 <- as.data.frame(resOrdered)
write.table(DEG_DESeq2,file = "mrna_DEG_DESeq2.txt",sep = "\t")

#mrna火山图，设置筛选条件，自定义cutoff
input <- read.table("mrna_DEG_DESeq2.txt",header = T,sep = "\t")
input$gene <- rownames(input)
input$Threshold=as.factor(ifelse(input$padj<0.05 & abs(input$log2FoldChange)>=1,
                                 ifelse(input$log2FoldChange>1,'Up','Down'),'No')) 
#查看log2FoldChange的最大值和最小值来确定下面xlim的参数
min(input$log2FoldChange)
max(input$log2FoldChange)

#参数说明
#geom_point(alpha=0.3, size=3) ：alpha调节点的透明度，size调节点的大小
#xlim(-10,10)：调节横坐标边界
#保存为pdf格式
#pdf(file = 'volcano.pdf',height = 5,width = 5)
##down top 10
down <- filter(input, Threshold == "Down") %>% 
  distinct(gene, .keep_all = T) %>%
  top_n(10, -log2FoldChange)
##up top 10
up <- filter(input, Threshold == "Up") %>% 
  distinct(gene, .keep_all = T) %>%
  top_n(10, log2FoldChange)

nudge_x_up = 7 - up$log2FoldChange
nudge_x_down = -5 - down$log2FoldChange
library(ggrepel)
ggplot(data = input,aes(x = log2FoldChange,y=-log10(padj),
                        colour=Threshold)) +geom_point(alpha=0.3, size=3) +
  scale_color_manual(values=c("blue", "grey","red"))+
  theme_bw()+xlim(-5.5,7.5)+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = NULL, colour = "black"))+
  geom_vline(xintercept=c(-1,1),lty=2,col="black",lwd=0.4) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.4) +
  #'@添加关注的点的基因名
  #'@添加down top gene
  geom_text_repel(
    data = up,aes(x = log2FoldChange, y = -log10(padj), label = gene),
    seed = 123,color = 'black',show.legend = FALSE, 
    min.segment.length = 0,#始终为标签添加指引线段；若不想添加线段，则改为Inf
    segment.linetype = 1, #线段类型,1为实线,2-6为不同类型虚线
    segment.color = 'black', #线段颜色
    segment.alpha = 0.5, #线段不透明度
    nudge_x = nudge_x_up, #标签x轴起始位置调整
    direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
    hjust = 0, #对齐标签：0右对齐，1左对齐，0.5居中
    force = 2,#重叠标签间的排斥力
    force_pull = 2,#标签和数据点间的吸引力
    size = 4,
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    max.overlaps = Inf)+
  ##'@添加up top gene
  geom_text_repel(
    data = down,aes(x = log2FoldChange, y = -log10(padj), label = gene),
    seed = 123,color = 'black',show.legend = FALSE, 
    min.segment.length = 0,#始终为标签添加指引线段；若不想添加线段，则改为Inf
    segment.linetype = 1, #线段类型,1为实线,2-6为不同类型虚线
    segment.color = 'black', #线段颜色
    segment.alpha = 0.5, #线段不透明度
    nudge_x = nudge_x_down, #标签x轴起始位置调整
    direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
    hjust = 1, #对齐标签：0右对齐，1左对齐，0.5居中
    force = 2,#重叠标签间的排斥力
    force_pull = 2,#标签和数据点间的吸引力
    size = 4,
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    max.overlaps = Inf) 

##03wgcna判断gem敏感性相关mrna和lncrna
setwd("/home/pengdu/gem/03wgcna")
library(WGCNA)
#mrna的wgcna
Gemcitabine <- read.table("gem.txt",header = T)
sample <- Gemcitabine$sample 
sample <- sample[which(sample != "TCGA.ZF.A9RD.01")] 
mrna <- read.table("mrna.txt",header = T)
mrna <- as.data.frame.matrix(mrna)
write.table(colnames(mrna),"a.txt",sep = "\t")
mrna <- mrna[,sample]
write.table(mrna,"mrna.txt",sep = "\t")
# 提取表达矩阵并行列转置，将行名变为样本名
datExpr0 <- as.data.frame(t(mrna))
rownames(datExpr0) <- names(mrna)
# 判断矩阵的样本是否都合格？
gsg <- goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,9) #视图
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
traitData = read.table("gem.txt",header = T)
rownames(traitData) <- traitData[,1]
traitData = traitData[sample,]
traitData = as.data.frame(traitData)
dim(traitData)  #每行是一个样本，每列是一种信息
sensitivity = rownames(datExpr0)
traitRows = match(sensitivity, traitData$sample);
datTraits = traitData[traitRows, -1]
datTraits <- as.numeric(datTraits)
datTraits <- as.matrix(datTraits)
rownames(datTraits) = traitData[traitRows, 1];
colnames(datTraits) = "gemcitabine sensitivity"
collectGarbage()
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) #用颜色代表关联度
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(datTraits),
                    main = "Sample dendrogram and trait heatmap")
powers = c(1:20)
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.9,col="red")  #查看位于0.9以上的点，可以改变高度值
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
sft$powerEstimate
#结果为11
net = blockwiseModules(datExpr0, power = 11,
                       TOMType = "unsigned", minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.1,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "lncrnagemTOM",
                       verbose = 3)
table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[2]], mergedColors[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
# 重新计算带有颜色标签的模块
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 通过相关值对每个关联进行颜色编码
sizeGrWindow(10,6)
# 展示模块与表型数据的相关系数和 P值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
colnames(datTraits)
# 用热图的形式展示相关系数
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#colors = greenWhiteRed(50)不适用于红绿色盲患者，建议用 blueWhiteRed代替.
#基因与表型数据的关系、重要模块：基因显著性和模块成员
sensitivity = as.data.frame(datTraits);
names(sensitivity) = "gemcitabine sensitivity";
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, sensitivity, use = "p"));#和体重性状的关联
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(sensitivity), sep="");
names(GSPvalue) = paste("p.GS.", names(sensitivity), sep="");
#运行以下代码可视化GS和MM
#module = c("pink","red","blue","turquoise","brown","green","black","yellow","grey") 
module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
table(moduleGenes)
brown_module <-as.data.frame(dimnames(data.frame(datExpr0))[[2]][moduleGenes]) 
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
MM <- abs(geneModuleMembership[moduleGenes, column])
GS <- abs(geneTraitSignificance[moduleGenes, 1])
brown <-as.data.frame(cbind(MM,GS)) #包含了MM和GS的数据，可以保留一下
rownames(brown)=brown_module[,1]
#green_hub <-abs(c$MM)>0.8&abs(c$GS)>0.2 #筛选hub基因 
write.csv(brown, "hubgene_MMGS_mrna_brown.csv") 
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for sensitivity",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "brown")
abline(h=0.3,v=0.5,col="black",lwd=1.5) 
names(datExpr0)#会返回所有在分析中的基因ID
brown <- names(datExpr0)[moduleColors=="brown"]#返回属于棕色模块的基因ID

module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;
table(moduleGenes)
green_module <-as.data.frame(dimnames(data.frame(datExpr0))[[2]][moduleGenes]) 
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
MM <- abs(geneModuleMembership[moduleGenes, column])
GS <- abs(geneTraitSignificance[moduleGenes, 1])
green <-as.data.frame(cbind(MM,GS)) #包含了MM和GS的数据，可以保留一下
rownames(green)=green_module[,1]
#green_hub <-abs(c$MM)>0.8&abs(c$GS)>0.2 #筛选hub基因 
write.csv(green, "hubgene_MMGS_mrna_green.csv") 
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for sensitivity",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "green")
abline(h=0.3,v=0.5,col="black",lwd=1.5) 
names(datExpr0)#会返回所有在分析中的基因ID
green <- names(datExpr0)[moduleColors=="green"]#返回属于棕色模块的基因ID

#lncrna的wgcna
Gemcitabine <- read.table("gem.txt",header = T)
sample <- Gemcitabine$sample 
sample <- sample[which(sample != "TCGA.ZF.A9RD.01")] 
lncrna <- read.table("lncrna.txt",header = T)
lncrna <- as.data.frame.matrix(lncrna)
write.table(colnames(lncrna),"a.txt",sep = "\t")
lncrna <- lncrna[,sample]
write.table(lncrna,"lncrna.txt",sep = "\t")
# 提取表达矩阵并行列转置，将行名变为样本名
datExpr0 <- as.data.frame(t(lncrna))
rownames(datExpr0) <- names(lncrna)
# 判断矩阵的样本是否都合格？
gsg <- goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,9) #视图
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
traitData = read.table("gem.txt",header = T)
rownames(traitData) <- traitData[,1]
traitData = traitData[sample,]
traitData = as.data.frame(traitData)
dim(traitData)  #每行是一个样本，每列是一种信息
sensitivity = rownames(datExpr0)
traitRows = match(sensitivity, traitData$sample);
datTraits = traitData[traitRows, -1]
datTraits <- as.numeric(datTraits)
datTraits <- as.matrix(datTraits)
rownames(datTraits) = traitData[traitRows, 1];
colnames(datTraits) = "gemcitabine sensitivity"
collectGarbage()
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) #用颜色代表关联度
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(datTraits),
                    main = "Sample dendrogram and trait heatmap")
powers = c(1:20)
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.8;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.8,col="red")  #查看位于0.9以上的点，可以改变高度值
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
sft$powerEstimate
#结果为11
net = blockwiseModules(datExpr0, power = 3,
                       TOMType = "unsigned", minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.1,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "lncrnagemTOM",
                       verbose = 3)
table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[2]], mergedColors[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
# 重新计算带有颜色标签的模块
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 通过相关值对每个关联进行颜色编码
sizeGrWindow(10,6)
# 展示模块与表型数据的相关系数和 P值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
colnames(datTraits)
# 用热图的形式展示相关系数
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#colors = greenWhiteRed(50)不适用于红绿色盲患者，建议用 blueWhiteRed代替.
#基因与表型数据的关系、重要模块：基因显著性和模块成员
sensitivity = as.data.frame(datTraits);
names(sensitivity) = "gemcitabine sensitivity";
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, sensitivity, use = "p"));#和体重性状的关联
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(sensitivity), sep="");
names(GSPvalue) = paste("p.GS.", names(sensitivity), sep="");
module = "magenta"
column = match(module, modNames);
moduleGenes = moduleColors==module;
table(moduleGenes)
magenta_module <-as.data.frame(dimnames(data.frame(datExpr0))[[2]][moduleGenes]) 
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
MM <- abs(geneModuleMembership[moduleGenes, column])
GS <- abs(geneTraitSignificance[moduleGenes, 1])
magenta <-as.data.frame(cbind(MM,GS)) #包含了MM和GS的数据，可以保留一下
rownames(magenta)=magenta_module[,1]
#magenta_hub <-abs(c$MM)>0.8&abs(c$GS)>0.2 #筛选hub基因 
write.csv(magenta, "hubgene_MMGS_lncrna_magenta.csv") 
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for sensitivity",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "magenta")
abline(h=0.3,v=0.5,col="black",lwd=1.5) 
names(datExpr0)#会返回所有在分析中的基因ID
magenta <- names(datExpr0)[moduleColors=="magenta"]#返回属于棕色模块的基因ID

##04mrna的go/kegg富集分析
setwd("/home/pengdu/gem/04gokegg")
#取交集lncrna
options(stringsAsFactors = F)
dat <- read.csv("hubgene_MMGS_lncrna_magenta.csv",header=T)
lncRNA_magenta <- dat$lncRNA_magenta[1:188]
DEGs <- dat$DEGs[1:11100]
library(VennDiagram)
venn_list<-list(lncRNA_magenta = lncRNA_magenta,DEGs = DEGs)
venn.diagram(venn_list,filename = 'venn.png',imagetype = 'png',scaled = FALSE,
             fill = c('orange','blue'),alpha =0.50,
             cat.col = c('orange','blue'),cat.cex = 1.5,cat.pos = 1,cat.fontfamily = 'serif',
             col = "black",cex = 1.5, fontfamily = 'serif',
             margin = 0.1)
inter<- get.venn.partitions(venn_list)
for(i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']],collapse = '\t')
write.table(inter[,c(5,6)],'venn_inter_lncrna.txt',row.names = FALSE , sep = '\t',quote = FALSE)
#取交集mrna
dat <- read.csv("mrna_DEG_DESeq2.csv",header=T)
mRNA_brown_green <- dat$mRNA_brown_green[1:620]
DEGs <- dat$DEGs[1:4603]
library(VennDiagram)
venn_list<-list(mRNA_brown_green = mRNA_brown_green,DEGs = DEGs)
venn.diagram(venn_list,filename = 'venn_mrna.png',imagetype = 'png',scaled = FALSE,
             fill = c('orange','blue'),alpha =0.50,
             cat.col = c('orange','blue'),cat.cex = 1.5,cat.pos = 1,cat.fontfamily = 'serif',
             col = "black",cex = 1.5, fontfamily = 'serif',
             margin = 0.1)
inter<- get.venn.partitions(venn_list)
for(i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']],collapse = '\t')
write.table(inter[,c(5,6)],'venn_inter_mrna.txt',row.names = FALSE , sep = '\t',quote = FALSE)

library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
info <- read.table("mrna_DEG_DESeq2.txt",header = T,sep = "\t")
#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html
#gene ID转换
gene <- bitr(info$gene_symbol,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
GO<-enrichGO(gene$ENTREZID,#GO富集分析
             OrgDb = GO_database,
             keyType = "ENTREZID",#设定读取的gene ID类型
             ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
             pvalueCutoff = 0.05,#设定p值阈值
             qvalueCutoff = 0.05,#设定q值阈值
             readable = T)
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
dotplot(GO, x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
        showCategory =5, #只显示前5
        split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") #点状图
dotplot(KEGG, x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
        showCategory =5)

##05 101种机器学习方法 
getwd()
setwd("/home/pengdu/gem/05ml/Figures")
bmt <-read.csv('survivalplot.csv',header=T,sep = ",")
#TCGA_BLCA
blca <- subset(bmt,Cohort == "TCGA_BLCA")
library(survminer)
#OS
# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(blca, #数据集
                         time = "OS.time", #生存状态
                         event = "OS", #生存时间
                         variables = c("Riskscore") #需要计算的数据列名
)
summary(res.cut) #查看数据最佳截断点及统计量
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS == 1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="TCGA_BLCA_cohort",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)
#RFS
# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(blca, #数据集
                         time = "RFS.time", #生存状态
                         event = "Recurrence", #生存时间
                         variables = c("Riskscore") #需要计算的数据列名
)
summary(res.cut) #查看数据最佳截断点及统计量
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(RFS.time,Recurrence == 1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Recurrence free survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="TCGA_BLCA_cohort",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)
#GSE31684 
GSE31684 <- subset(bmt,Cohort == "GSE31684")
library(survminer)
#OS
# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(GSE31684, #数据集
                         time = "OS.time", #生存状态
                         event = "OS", #生存时间
                         variables = c("Riskscore") #需要计算的数据列名
)
summary(res.cut) #查看数据最佳截断点及统计量
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(OS.time,OS == 1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="GSE31684_cohort",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)
#RFS
# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(GSE31684, #数据集
                         time = "RFS.time", #生存状态
                         event = "Recurrence", #生存时间
                         variables = c("Riskscore") #需要计算的数据列名
)
summary(res.cut) #查看数据最佳截断点及统计量
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(RFS.time,Recurrence == 1)~Riskscore,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Recurrence free survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="GSE31684_cohort",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

#roc曲线绘制
bmt <-read.csv('survivalplot.csv',header=T,sep = ",") 
library(survival)
library(survivalROC)
library(timeROC)
blca <- subset(bmt,Cohort == "TCGA_BLCA")
#tcga os
tROC <-timeROC(T = blca$OS.time,delta = blca$OS,marker = blca$Riskscore,
               cause = 1,times = c(365,730,1095),ROC=T)
#开始画图
plot(tROC,time=365,col="red",title=F,lwd=2) #1年ROC
plot(tROC,time=730,col="green",add=T,title=F,lwd=2) #3年ROC
plot(tROC,time=1095,col="darkblue",add=T,title=F,lwd=2) #5年ROC
legend(0.5,0.5, # 这里设置legend的位置坐标：横坐标0,纵坐标1。也可以像前面一样设置为"bottomright"等
       c(paste0("OS at 365 days  ", round(tROC$AUC[1], 2)),
         paste0("OS at 730 days  ", round(tROC$AUC[2], 2)),
         paste0("OS at 1095 days  ", round(tROC$AUC[3], 2))),
       col=c("red","green","darkblue"),lwd=2,cex=1.2,bty="n",title = "TCGA-BLCA")
#tcga rfs
tROC <-timeROC(T = blca$RFS.time,delta = blca$Recurrence,marker = blca$Riskscore,
               cause = 1,times = c(365,730,1095),ROC=T)
#开始画图
plot(tROC,time=365,col="red",title=F,lwd=2) #1年ROC
plot(tROC,time=730,col="green",add=T,title=F,lwd=2) #3年ROC
plot(tROC,time=1095,col="darkblue",add=T,title=F,lwd=2) #5年ROC
legend(0.5,0.5, # 这里设置legend的位置坐标：横坐标0,纵坐标1。也可以像前面一样设置为"bottomright"等
       c(paste0("RFS at 365 days  ", round(tROC$AUC[1], 2)),
         paste0("RFS at 730 days  ", round(tROC$AUC[2], 2)),
         paste0("RFS at 1095 days  ", round(tROC$AUC[3], 2))),
       col=c("red","green","darkblue"),lwd=2,cex=1.2,bty="n",title = "TCGA-BLCA")
GSE31684 <- subset(bmt,Cohort == "GSE31684")
#GSE31684 os
tROC <-timeROC(T = GSE31684$OS.time,delta = GSE31684$OS,marker = GSE31684$Riskscore,
               cause = 1,times = c(365,730,1095),ROC=T)
#开始画图
plot(tROC,time=365,col="red",title=F,lwd=2) #1年ROC
plot(tROC,time=730,col="green",add=T,title=F,lwd=2) #3年ROC
plot(tROC,time=1095,col="darkblue",add=T,title=F,lwd=2) #5年ROC
legend(0.5,0.5, # 这里设置legend的位置坐标：横坐标0,纵坐标1。也可以像前面一样设置为"bottomright"等
       c(paste0("OS at 365 days  ", round(tROC$AUC[1], 2)),
         paste0("OS at 730 days  ", round(tROC$AUC[2], 2)),
         paste0("OS at 1095 days  ", round(tROC$AUC[3], 2))),
       col=c("red","green","darkblue"),lwd=2,cex=1.2,bty="n",title = "GSE31684")
#GSE31684 rfs
tROC <-timeROC(T = GSE31684$RFS.time,delta = GSE31684$Recurrence,marker = GSE31684$Riskscore,
               cause = 1,times = c(365,730,1095),ROC=T)
#开始画图
plot(tROC,time=365,col="red",title=F,lwd=2) #1年ROC
plot(tROC,time=730,col="green",add=T,title=F,lwd=2) #3年ROC
plot(tROC,time=1095,col="darkblue",add=T,title=F,lwd=2) #5年ROC
legend(0.5,0.5, # 这里设置legend的位置坐标：横坐标0,纵坐标1。也可以像前面一样设置为"bottomright"等
       c(paste0("RFS at 365 days  ", round(tROC$AUC[1], 2)),
         paste0("RFS at 730 days  ", round(tROC$AUC[2], 2)),
         paste0("RFS at 1095 days  ", round(tROC$AUC[3], 2))),
       col=c("red","green","darkblue"),lwd=2,cex=1.2,bty="n",title = "GSE31684")

##06 突变组对比
getwd()
setwd("/home/pengdu/gem/06mutation")
library(maftools)
blca.maf <- TCGAmutations::tcga_available()
blca.maf <- TCGAmutations::tcga_load(study = "blca")
blca <- blca.maf@data
write.table(blca,"blca_mutation.txt",sep = "\t")
risk <- read.table("clinical.txt",sep = "\t",header = T)
#high组
high <- risk[risk$group == "high",]
high <- high$sample
high <- subset(blca, Tumor_Sample_Barcode_min %in% high)
maf.high <- read.maf(maf = high,vc_nonSyn=names(tail(sort(table(high$Variant_Classification)))))
##read.maf读取以制表符作为分隔的maf文件
vc_cols <- RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) <- c(
  'In_Frame_Ins',
  'Translation_Start_Site',
  'Nonstop_Mutation',
  'Frame_Shift_Ins',
  'Frame_Shift_Del',
  'Splice_Site',
  'In_Frame_Del',
  'Missense_Mutation'
)
oncoplot(maf=maf.high,colors = vc_cols,top=10,fontSize=1,showTumorSampleBarcodes=F) 
plotmafSummary(maf = maf.high,
               rmOutlier = T,
               addStat = "median",
               dashboard = T,
               titvRaw = F)
#low组
low <- risk[risk$group == "low",]
low <- low$sample
low <- subset(blca, Tumor_Sample_Barcode_min %in% low)
maf.low <- read.maf(maf = low,vc_nonSyn=names(tail(sort(table(low$Variant_Classification)))))
##read.maf读取以制表符作为分隔的maf文件
vc_cols <- RColorBrewer::brewer.pal(n = 8, name = 'Accent')
names(vc_cols) <- c(
  'In_Frame_Ins',
  'Translation_Start_Site',
  'Nonstop_Mutation',
  'Frame_Shift_Ins',
  'Frame_Shift_Del',
  'Splice_Site',
  'In_Frame_Del',
  'Missense_Mutation'
)
oncoplot(maf=maf.low,colors = vc_cols,top=10,fontSize=1,showTumorSampleBarcodes=F) 
plotmafSummary(maf = maf.low,
               rmOutlier = T,
               addStat = "median",
               dashboard = T,
               titvRaw = F)
#total组
high$group <- "high"
low$group <- "low"
total <- rbind(high,low)
maf.total <- read.maf(maf = total,clinicalData = total,
                      vc_nonSyn=names(tail(sort(table(total$Variant_Classification)))))
##read.maf读取以制表符作为分隔的maf文件
col_group=c("#E7B800","#2E9FDF") 
assign_group=setNames(col_group,unique(total$group))
oncoplot(maf=maf.total,top=10,fontSize=1,showTumorSampleBarcodes=F,
         clinicalFeatures = 'group',
         annotationColor = list(group=assign_group),
         sortByAnnotation = TRUE) 
plotmafSummary(maf = maf.total,
               rmOutlier = T,
               addStat = "median",
               dashboard = T,
               titvRaw = F)
#使用mafCompare比较差异突变基因
fvsm <- mafCompare(m1=maf.high, m2=maf.low, m1Name="high", m2Name="low", minMut=10)
result <- fvsm$results
write.table(result,"highvslow.txt",sep = "\t",na = "",row.names = F)
forestPlot(mafCompareRes=fvsm, pVal=0.01, color=c("#2E9FDF","#E7B800"), geneFontSize=0.8)
#共突变和互斥突变
Interact <- somaticInteractions(maf = maf.total, top = 20, pvalue = c(0.05, 0.1))
write.csv(Interact,"Interact.csv")
#tmb
tmb_table = tmb(maf = blca.maf)   #默认以log10转化的TMB绘图
write.table(tmb_table,"tmb.txt",sep = "\t")
tmb <- read.table("tmb.txt",sep = "\t",header = T,row.names = 1)
colnames(tmb)
library(ggpubr)
my_comparisons <- list(c("low","high"))
ggviolin(tmb, x="group", y="total_perMB_log", color = "group", 
         ylab="total_perMB_log",
         xlab="group",
         # legend.title=x,
         add.params = list(fill="white"),
         palette = c("#ECE700","#5bc0eb"),    
         width=1, add = "boxplot") +
  stat_compare_means(aes(group=group),comparisons = my_comparisons,
                     method="wilcox.test", # 可换其他统计方法
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif") + 
  ylim(-2,2.5) 

##07 high vs low
#gsea
setwd("/home/pengdu/gem/07highlow/gsea")
rt = data.table::fread("TCGA-BLCA.star_counts.tsv.gz", header=T, sep="\t", check.names=F)
rt <- as.data.frame.matrix(rt)
geneid <- rt[,1]
# ID转换
library(dplyr)
probeMap <- read.delim("gencode.v36.annotation.gtf.gene.probemap")
rt <- rt %>%
  inner_join(probeMap, by = c("Ensembl_ID" = "id")) %>%
  select(gene, starts_with("TCGA"))
library(limma)
rt <- as.data.frame(avereps(rt[, -1], ID = rt$gene)) # 取平均
dim(rt)
write.table(rt,"counts.txt",sep = "\t")
#high vs low
group <- read.table("group.txt",header = T, sep="\t", check.names=F)
sample <- intersect(group$sample,colnames(rt))
rt <- rt[,sample]
rownames(group) <- group$sample
group <- group[sample,]
table(group$group)
group <- c(rep('high',181), rep('low',219))
group <- factor(group)
table(group)
# 构建DESeq2中的对象
library(DESeq2)
colData <- data.frame(row.names = colnames(rt),group = group)
colData$group <- factor(colData$group, levels = c("high", "low"))
head(colData)
data_int <- 2^(rt) - 1
data_int <- apply(data_int, 2, as.integer)
rownames(data_int) <- rownames(rt)
dds <- DESeqDataSetFromMatrix(countData = data_int,colData = colData,design = ~ group)
dds <- DESeq(dds) # 进行差异表达分析
resultsNames(dds) # 查看结果的名称
# contrast参数必须写成下面三个元素的向量格式，且顺序不能反
res <- results(dds, contrast = c("group", rev(levels(group))))
resOrdered <- res[order(res$padj), ]
DEG_DESeq2 <- as.data.frame(resOrdered)
write.table(DEG_DESeq2,file = "low_vs_high_DEG_DESeq2.txt",sep = "\t")
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
info <- read.table("low_vs_high_DEG_DESeq2.txt",header = T,sep = "\t",row.names = 1)
info$SYMBOL <- row.names(info)
gene <- bitr(info$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
info_merge <- merge(info,gene,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange
GSEA_input <- info_merge$log2FoldChange
names(GSEA_input) = info_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- gseKEGG(GSEA_input, organism = "hsa", pvalueCutoff = 0.05)#GSEA富集分析
result <- GSEA_KEGG@result
write.csv(result,"GSEA_result.csv")
GSEA_result <- read.csv("GSEA_result.csv",header = T,sep = ",",row.names = 1)
GSEA_KEGG@result <- GSEA_result 
gseaplot2(GSEA_KEGG,geneSetID = 1:5)
gseaplot2(GSEA_KEGG,geneSetID = 6:10)

#infiltration 
setwd("/home/pengdu/gem/07highlow/infiltration")
immune <- read.csv("infiltration_estimation_for_tcga.csv",header = T)
clin <- immune[,1:2]
rownames(clin) <- clin[,1]
clin <- clin[,-1]
clin <- as.data.frame(clin)
rownames(clin) <- immune[,1]
colnames(clin) <- c("group")
exp <- immune[,c(1,3:121)]
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp <- t(exp)
identical(rownames(clin),colnames(exp)) #计算显著性
table(clin$group)
#把之前的结果填到这里
library(limma)
group_list <- factor(rep(c('low','high'),c(222,183)))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exp)
contrast.matrix <- makeContrasts('low-high',levels = design)
fit <- lmFit(exp,design = design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
#输出全部差异结果
alldiff=topTable(fit2,coef = 1,n = Inf) 
#判断星号
alldiff$lab = as.factor(ifelse(alldiff$P.Value>=0.05,"",
                               
                               ifelse(alldiff$P.Value>=0.01&alldiff$P.Value<0.05,"*",
                                      
                                      ifelse(alldiff$P.Value>=0.001&alldiff$P.Value<0.01,"**", "***"
                                             
                                      ))))
#变量名与显著性合在一起
alldiff$new <- paste(rownames(alldiff),alldiff$lab)
#热图数据替换
alldiff <- alldiff[rownames(exp),]
rownames(exp) <- alldiff$new 
#画图
pre_heatdata <- t(scale(t(exp)))
pre_heatdata[pre_heatdata> 1] <- 1###这里的数值自己调，哪个好看就怎么来
pre_heatdata[pre_heatdata< -1] <- -1###这里的数值自己调，哪个好看就怎么来
annColors <- list()
annColors[['group']] <- c('low'="#5bc0eb",'high'="#ECE700")
library(pheatmap)
#pdf(file="pheatmap.pdf",width = 10,height = 12)
annotation_row = data.frame(
  infiltration = c(rep("TIMER",6), rep("CIBERSORT", 22), rep("CIBERSORT-ABS", 22),
                   rep("QUANTISEQ", 11), rep("MCPCOUNTER", 11), rep("XCELL", 39),
                   rep("EPIC", 8))
)
row.names(annotation_row) <- rownames(pre_heatdata)
pheatmap(pre_heatdata,
         color = colorRampPalette(c("#5bc0eb",'white',"#ECE700"))(1000),
         annotation_col = clin,
         annotation_colors = annColors,
         treeheight_row = 50,
         gaps_col = 222,
         gaps_row = c(6,28,50,61,72,111),
         show_rownames = T,
         show_colnames = F,
         annotation_row = annotation_row,
         cluster_rows = F,
         cluster_cols = F)
dev.off() 
#immunecheckpoint
checkpoint <- rt[c("PDCD1","CD274","CTLA4"),]
group <- read.table("group.txt",header = T, sep="\t", check.names=F)
sample <- intersect(group$sample,colnames(checkpoint))
checkpoint <- checkpoint[,sample]
rownames(group) <- group$sample
group <- group[sample,]
table(group$group)
checkpoint <- t(checkpoint)
checkpoint <- as.data.frame(checkpoint)
checkpoint <- cbind(group,checkpoint)
checkpoint <- checkpoint[,-1]
checkpoint <- reshape2::melt(checkpoint,id.vars=c("group"))
colnames(checkpoint)=c("group","checkpoint","Relative expression")
library(ggpubr)
my_comparisons <- list(c("low","high"))
ggviolin(checkpoint, x="checkpoint", y="Relative expression", color = "group", 
         ylab="Relative expression",
         xlab="checkpoint",
         # legend.title=x,
         add.params = list(fill="white"),
         palette = c("#ECE700","#5bc0eb"),    
         width=1, add = "boxplot") +
  stat_compare_means(aes(group=group), #comparisons = my_comparisons,
                     method="wilcox.test", # 可换其他统计方法
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif") +rotate_x_text(60)
#临床信息
setwd("/home/pengdu/gem/07highlow/clinical")
library("ggpubr")
# Tstage
Tstage <- read.csv("Tstage.csv",header = T,sep = ",",check.names = F)
mycom <- list(c("T1","T3"),c("T2","T3"),c("T1","T4"))
ggviolin(Tstage, x = "T stage",y = "Riskscore",color = "T stage", size = 1,
         palette = c("#B0E854","#DEA7CC","#79ACD5","#D99264","#ECE700"),
         add = "boxplot") +
  stat_compare_means(comparisons = mycom, label.y = c(2.85,2.5,3.2),
                     label = "p.signif",size = 6)+
  stat_compare_means(label.y = 3.6)+
  ylim(0.5,3.9) 
# AJCCstage
AJCCstage <- read.csv("AJCCstage.csv",header = T,sep = ",",check.names = F)
mycom <- list(c("I/II","III/IV"))
ggviolin(AJCCstage, x = "AJCC stage",y = "Riskscore",color = "AJCC stage", size = 1.2,
         palette = c("#B0E854","#DEA7CC","#79ACD5","#D99264","#ECE700"),
         add = "boxplot") +
  stat_compare_means(comparisons = mycom, label.y = c(2.8),
                     label = "p.signif",size = 6)+
  stat_compare_means(label.y = 3.6)+
  ylim(0.5,3.9) 
# grade
grade <- read.csv("grade.csv",header = T,sep = ",",check.names = F)
mycom <- list(c("low","high"))
ggviolin(grade, x = "grade",y = "Riskscore",color = "grade", size = 1.2,
         palette = c("#B0E854","#DEA7CC","#79ACD5","#D99264","#ECE700"),
         add = "boxplot") +
  stat_compare_means(comparisons = mycom, label.y = c(2.8),
                     label = "p.signif",size = 6)+
  stat_compare_means(label.y = 3.6)+
  ylim(0.5,3.9) 
# gross
gross <- read.csv("gross.csv",header = T,sep = ",",check.names = F)
table(gross$`gross subtype`)
mycom <- list(c("Non-Papillary","Papillary"))
ggviolin(gross, x = "gross subtype",y = "Riskscore",color = "gross subtype", size = 1.2,
         palette = c("#B0E854","#DEA7CC","#79ACD5","#D99264","#ECE700"),
         add = "boxplot") +
  stat_compare_means(comparisons = mycom, label.y = c(2.8),
                     label = "p.signif",size = 6)+
  stat_compare_means(label.y = 3.6)+
  ylim(0.5,3.9) 
# molecular subtype
molecular <- read.csv("molecular subtype.csv",header = T,sep = ",",check.names = F)
table(molecular$`molecular subtype`)
mycom <- list(c("BLCA.1","BLCA.3"),c("BLCA.2","BLCA.3"),c("BLCA.2","BLCA.4"))
ggviolin(molecular, x = "molecular subtype",y = "Riskscore",color = "molecular subtype", size = 1.2,
         palette = c("#B0E854","#DEA7CC","#79ACD5","#D99264","#ECE700"),
         add = "boxplot") +
  stat_compare_means(comparisons = mycom, label.y = c(2.5,2.8,3.1),
                     label = "p.signif",size = 6)+
  stat_compare_means(label.y = 3.6)+
  ylim(0.5,3.9) 
# Subtype_Immune_Model_Based
subtype <- read.csv("Subtype_Immune_Model_Based.csv",header = T,sep = ",",check.names = F)
table(subtype$`immune model subtype`)
mycom <- list(c("Wound Healing (Immune C1)","IFN-gamma Dominant (Immune C2)"),
              c("IFN-gamma Dominant (Immune C2)","Inflammatory (Immune C3)"),
              c("IFN-gamma Dominant (Immune C2)","Lymphocyte Depleted (Immune C4)"))
ggviolin(subtype, x = "immune model subtype",y = "Riskscore",color = "immune model subtype", size = 1.2,
         palette = c("#B0E854","#DEA7CC","#79ACD5","#D99264","#ECE700"),
         add = "boxplot") +
  stat_compare_means(comparisons = mycom, label.y = c(2.5,2.8,3.1),
                     label = "p.signif",size = 6)+
  stat_compare_means(label.y = 3.6)+
  ylim(0.5,3.9) 
#dnass
dnass <- read.csv("dnass.csv",header = T,sep = ",",row.names = 1)
dnass <- reshape2::melt(dnass,id.vars=c("group"))
colnames(dnass)=c("group","stemness","stemness score")
library(ggpubr)
my_comparisons <- list(c("low","high"))
ggviolin(dnass, x="stemness", y="stemness score", color = "group", 
         ylab="stemness score",
         xlab="stemness",
         # legend.title=x,
         add.params = list(fill="white"),
         palette = c("#ECE700","#5bc0eb"),    
         width=1, add = "boxplot") +
  stat_compare_means(aes(group=group), #comparisons = my_comparisons,
                     method="wilcox.test", # 可换其他统计方法
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif") +rotate_x_text(60)
#ImmuneSigs 
setwd("/home/pengdu/gem/07highlow/infiltration")
immune <- read.csv("TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.csv",header = T)
clin <- immune[,1:2]
rownames(clin) <- clin[,1]
clin <- clin[,-1]
clin <- as.data.frame(clin)
rownames(clin) <- immune[,1]
colnames(clin) <- c("group")
exp <- immune[,c(1,3:70)]
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp <- t(exp)
identical(rownames(clin),colnames(exp)) #计算显著性
table(clin$group)
#把之前的结果填到这里
library(limma)
group_list <- factor(rep(c('low','high'),c(219,182)))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exp)
contrast.matrix <- makeContrasts('low-high',levels = design)
fit <- lmFit(exp,design = design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
#输出全部差异结果
alldiff=topTable(fit2,coef = 1,n = Inf) 
#判断星号
alldiff$lab = as.factor(ifelse(alldiff$P.Value>=0.05,"",
                               
                               ifelse(alldiff$P.Value>=0.01&alldiff$P.Value<0.05,"*",
                                      
                                      ifelse(alldiff$P.Value>=0.001&alldiff$P.Value<0.01,"**", "***"
                                             
                                      ))))
#变量名与显著性合在一起
alldiff$new <- paste(rownames(alldiff),alldiff$lab)
#热图数据替换
alldiff <- alldiff[rownames(exp),]
rownames(exp) <- alldiff$new 
#画图
pre_heatdata <- t(scale(t(exp)))
pre_heatdata[pre_heatdata> 1] <- 1###这里的数值自己调，哪个好看就怎么来
pre_heatdata[pre_heatdata< -1] <- -1###这里的数值自己调，哪个好看就怎么来
annColors <- list()
annColors[['group']] <- c('low'="#5bc0eb",'high'="#ECE700")
library(pheatmap)
#pdf(file="pheatmap.pdf",width = 10,height = 12)
pheatmap(pre_heatdata,
         color = colorRampPalette(c("#5bc0eb",'white',"#ECE700"))(1000),
         annotation_col = clin,
         annotation_colors = annColors,
         treeheight_row = 50,
         gaps_col = 219,
         #gaps_row = c(6,28,50,61,72,111),
         show_rownames = T,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F)
dev.off() 
#pathway 
setwd("/home/pengdu/gem/07highlow/infiltration")
immune <- read.csv("Pancan12_GenePrograms_drugTargetCanon_in_Pancan33.csv",header = T)
clin <- immune[,1:2]
rownames(clin) <- clin[,1]
clin <- clin[,-1]
clin <- as.data.frame(clin)
rownames(clin) <- immune[,1]
colnames(clin) <- c("group")
exp <- immune[,c(1,3:41)]
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp <- t(exp)
identical(rownames(clin),colnames(exp)) #计算显著性
table(clin$group)
#把之前的结果填到这里
library(limma)
group_list <- factor(rep(c('low','high'),c(219,182)))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exp)
contrast.matrix <- makeContrasts('low-high',levels = design)
fit <- lmFit(exp,design = design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
#输出全部差异结果
alldiff=topTable(fit2,coef = 1,n = Inf) 
#判断星号
alldiff$lab = as.factor(ifelse(alldiff$P.Value>=0.05,"",
                               
                               ifelse(alldiff$P.Value>=0.01&alldiff$P.Value<0.05,"*",
                                      
                                      ifelse(alldiff$P.Value>=0.001&alldiff$P.Value<0.01,"**", "***"
                                             
                                      ))))
#变量名与显著性合在一起
alldiff$new <- paste(rownames(alldiff),alldiff$lab)
#热图数据替换
alldiff <- alldiff[rownames(exp),]
rownames(exp) <- alldiff$new 
#画图
pre_heatdata <- t(scale(t(exp)))
pre_heatdata[pre_heatdata> 1] <- 1###这里的数值自己调，哪个好看就怎么来
pre_heatdata[pre_heatdata< -1] <- -1###这里的数值自己调，哪个好看就怎么来
annColors <- list()
annColors[['group']] <- c('low'="#5bc0eb",'high'="#ECE700")
library(pheatmap)
#pdf(file="pheatmap.pdf",width = 10,height = 12)
pheatmap(pre_heatdata,
         color = colorRampPalette(c("#5bc0eb",'white',"#ECE700"))(1000),
         annotation_col = clin,
         annotation_colors = annColors,
         treeheight_row = 50,
         gaps_col = 219,
         #gaps_row = c(6,28,50,61,72,111),
         show_rownames = T,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F)
dev.off() 

##08 nomogram and web calculator
library(autoReg)
library(survival)
library(tidyverse)
setwd("/home/pengdu/gem/08nomogram")
blca <- read.csv("nomogram.csv",header = T,row.names = 1)
#cox回归模型构建
colnames(blca)
coxmod<-coxph(Surv(blca$OS,blca$Status == 1)~Age + Sex + Subtype + T + N + M + 
                   Grade + Smoke + Lymphovascular_invasion + Chemotherapy + Riskscore,
              data=blca)
summary(coxmod)
autoReg(coxmod)
result <- autoReg(coxmod,uni=TRUE, multi = TRUE, final=TRUE)
write.csv(result, "result.csv")
#多因素
blca <- read.csv("nomogram.csv",header = T,row.names = 1)
coxmod<-coxph(Surv(blca$OS,blca$Status == 1)~Age + T + N + M + 
                Lymphovascular_invasion + Chemotherapy + Riskscore,
              data=blca)
summary(coxmod)
#森林图单因素
#install.packages("forestplot")
library(forestplot)
rs_forest <- read.csv('senlin1.csv',header = FALSE,sep = ",")

# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。

forestplot(labeltext = as.matrix(rs_forest[,1:4]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V5, #设置均值
           
           lower = rs_forest$V6, #设置均值的lowlimits限
           
           upper = rs_forest$V7, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1,2,3),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,2.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box='red',summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 5)#设置森林图的位置，此处设置为4，则出现在第四列

rs_forest <- read.csv('senlin2.csv',header = FALSE,sep = ",")

# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。

forestplot(labeltext = as.matrix(rs_forest[,1:4]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V5, #设置均值
           
           lower = rs_forest$V6, #设置均值的lowlimits限
           
           upper = rs_forest$V7, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1,2,3),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,2.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box='green',summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 5)#设置森林图的位置，此处设置为4，则出现在第四列
#roc曲线
setwd("/home/pengdu/gem/08nomogram")
blca <- read.csv("nomogram.csv",header = T,row.names = 1)
#cox回归模型构建
coxmod <- coxph(Surv(blca$OS,blca$Status == 1)~Age + Sex + Subtype + T + N + M + 
                Grade + Smoke + Lymphovascular_invasion + Chemotherapy + Riskscore,
                data=blca)
blca$risk_score <- predict(coxmod,newdata = blca)
#校准曲线
library(rms)
rm(list = ls()) #清理环境
blca <- read.csv("nomogram.csv",header = T,row.names = 1)
nomo<-datadist(blca)
options(datadist='nomo')
nomo1 <- cph(Surv(blca$OS,blca$Status == 1)~risk_score,
             x=T,y=T,
             data=blca,
             surv=T,
             time.inc=365*5#示例数据time=月所以12*5就是评估5年的校准曲线
)#这里的time.inc一定要与下面画校准曲线的函数一致，不然图会出错！
p5<- calibrate(nomo1,#模型名称
              cmethod='KM',
              method='boot',#检测方法
              u=365*5,#评估的时间，注：一定要与模型的时间一致
              m=100, #每次抽样的样本量，
              B=1000)#抽样次数
#注，m值的确定：m=数据总数/3-4,即你想让最终的校准曲线有3个点，那就是m=数据总数/3
#B值一般1000，电脑配置不好可以选500,300,100等
plot(p5,
     add=F,#增加第二条线
     conf.int=T,#95%CI
     subtitles = F,#副标题
     cex.subtitles=0.8, #副标题大小
     lwd=2,#95%CI粗细
     lty=1,#95%CI实线，2=虚线
     errbar.col="blue",#95%CI颜色
     xlim=c(0.0,1),#x轴范围
     ylim=c(0.0,1),
     xlab="Prediced OS",
     ylab="Observed OS",
     col="darkblue")#曲线颜色
nomo2 <- cph(Surv(blca$OS,blca$Status == 1)~risk_score,
             x=T,y=T,
             data=blca,
             surv=T,
             time.inc=365*3)#示例数据time=月所以12*5就是评估5年的校准曲线
p3<- calibrate(nomo2,#模型名称
               cmethod='KM',
               method='boot',#检测方法
               u=365*3,#评估的时间，注：一定要与模型的时间一致
               m=100, #每次抽样的样本量，
               B=1000)#抽样次数
plot(p3,
     add=T,
     conf.int=T,#95%CI（蓝色线）
     subtitles = F,#关闭副标题
     cex.subtitles=0.8, 
     lwd=2,
     lty=1,
     errbar.col="blue",
     col="green")
nomo3 <- cph(Surv(blca$OS,blca$Status == 1)~risk_score,
             x=T,y=T,
             data=blca,
             surv=T,
             time.inc=365*1)#示例数据time=月所以12*5就是评估5年的校准曲线
p1<- calibrate(nomo3,#模型名称
               cmethod='KM',
               method='boot',#检测方法
               u=365*1,#评估的时间，注：一定要与模型的时间一致
               m=100, #每次抽样的样本量，
               B=1000)#抽样次数
plot(p1,
     add=T,
     conf.int=T,
     subtitles = F,
     cex.subtitles=0.8, 
     lwd=2,
     lty=1,
     errbar.col="blue",
     col="red")
#加上图例
legend("topleft", legend=c("1-year","3-year","5-year"), col=c("red", "green","darkblue"), lwd=2,
       inset= 0.02,box.lty = 0)
#调整对角线
abline(0,1,lty=3,lwd=1,col="grey")

#dca绘制
remotes::install_github('yikeshu0611/ggDCA')
library(rms)
coxmod<-cph(Surv(OS,Status == 1)~risk_score, data=blca)
dcacurve <- ggDCA::dca(coxmod,times=c(365,365*3,365*5))####多个模型5年后生存率
ggplot(dcacurve)

#tcga os
library(survival)
library(survivalROC)
library(timeROC)
tROC <-timeROC(T = blca$OS,delta = blca$Status,marker = blca$risk_score,
               cause = 1,times = c(365,1095,1825),ROC=T)
#开始画图
plot(tROC,time=365,col="red",title=F,lwd=2) #1年ROC
plot(tROC,time=1095,col="green",add=T,title=F,lwd=2) #3年ROC
plot(tROC,time=1825,col="darkblue",add=T,title=F,lwd=2) #5年ROC
legend(0.5,0.5, # 这里设置legend的位置坐标：横坐标0,纵坐标1。也可以像前面一样设置为"bottomright"等
       c(paste0("OS at 365 days  ", round(tROC$AUC[1], 2)),
         paste0("OS at 1095 days  ", round(tROC$AUC[2], 2)),
         paste0("OS at 1825 days  ", round(tROC$AUC[3], 2))),
       col=c("red","green","darkblue"),lwd=2,cex=1.2,bty="n",title = "TCGA-BLCA")
#网页计算器
library(DynNom)
library(shiny)
library(plotly)
library(compare)
library(stargazer)
library(survival)
library(rms)
setwd("/home/pengdu/gem/08nomogram")
training_dataset <- read.csv("calculator.txt",header = T,row.names = 1,sep = "\t")
colnames(training_dataset)
dd<-datadist(training_dataset)
options(datadist='dd')
f_cph <- cph(Surv(OS,Status == 1) ~ Riskscore+Chemotherapy+Lymphovascular_invasion+M+N+T+Age,                    ,
             x=T, y=T, surv=T,
             data=training_dataset)
print(f_cph)#查看多因素Cox分析结果
##2.0 尝试构建本地页面列线图
#原装代码如下
DynNom(f_cph,training_dataset,clevel = 0.95) #生成动态列线图
DNbuilder(model = f_cph,data = training_dataset,clevel = 0.95)

##09 immunotherapy
setwd("/home/pengdu/gem/09immunotherapy")
load("~/gem/09immunotherapy/IMvigor210CoreBiologies.Rdata")
#EMX2OS
EMX2OS <- subset(annoData,symbol == "EMX2OS")
EMX2OS <- EMX2OS$entrez_id
EMX2OS <- exprSet[EMX2OS,]
EMX2OS <- t(EMX2OS)
EMX2OS <- cbind(phenoData,EMX2OS)
colnames(EMX2OS)[26] <- "EMX2OS"
sum(is.na(EMX2OS))
write.csv(EMX2OS,"EMX2OS.csv")
EMX2OS <- read.csv("EMX2OS.csv",header = T,row.names = 1,sep = ",")
EMX2OS <- as.data.frame(EMX2OS)
library(survminer)
res.cut <- surv_cutpoint(EMX2OS, #数据集
                         time = "os", #生存状态
                         event = "censOS", #生存时间
                         variables = c("EMX2OS") #需要计算的数据列名
)
summary(res.cut) #查看数据最佳截断点及统计量
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(os,censOS == 1)~EMX2OS,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (months)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="BLCA samples in IMvigor210 cohort",ggtheme=theme_light(),xlim=c(0,24),break.x.by=6)
#LINC00930
LINC00930 <- subset(annoData,symbol == "LINC00930")
LINC00930 <- LINC00930$entrez_id
LINC00930 <- exprSet[LINC00930,]
LINC00930 <- t(LINC00930)
LINC00930 <- cbind(phenoData,LINC00930)
colnames(LINC00930)[26] <- "LINC00930"
sum(is.na(LINC00930))
write.csv(LINC00930,"LINC00930.csv")
LINC00930 <- as.data.frame(LINC00930)
library(survminer)
res.cut <- surv_cutpoint(LINC00930, #数据集
                         time = "os", #生存状态
                         event = "censOS", #生存时间
                         variables = c("LINC00930") #需要计算的数据列名
)
summary(res.cut) #查看数据最佳截断点及统计量
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(os,censOS == 1)~LINC00930,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="TCGA_BLCA_cohort",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)
#C17orf82
C17orf82 <- subset(annoData,symbol == "C17orf82")
C17orf82 <- C17orf82$entrez_id
C17orf82 <- exprSet[C17orf82,]
C17orf82 <- t(C17orf82)
C17orf82 <- cbind(phenoData,C17orf82)
colnames(C17orf82)[26] <- "C17orf82"
sum(is.na(C17orf82))
write.csv(C17orf82,"C17orf82.csv")
C17orf82 <- as.data.frame(C17orf82)
library(survminer)
res.cut <- surv_cutpoint(C17orf82, #数据集
                         time = "os", #生存状态
                         event = "censOS", #生存时间
                         variables = c("C17orf82") #需要计算的数据列名
)
summary(res.cut) #查看数据最佳截断点及统计量
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit<-survfit(Surv(os,censOS == 1)~C17orf82,res.cat)
print(fit)
ggsurvplot(fit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="TCGA_BLCA_cohort",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)
#binaryResponse
binaryResponse <- read.csv("binaryResponse.csv",header = T,sep = ",",row.names = 1)
data <- reshape2::melt(binaryResponse,id.vars=c("binaryResponse"))
data <- na.omit(data)
table(data$binaryResponse)
my_comparisons <- list(c("CR/PR","SD/PD"))
colnames(data) <- c("binaryResponse","Gene","Relative expression")
ggviolin(data, x="Gene", y="Relative expression", color = "binaryResponse", 
         ylab="Relative expression",
         xlab="Gene",
         # legend.title=x,
         add.params = list(fill="white"),
         palette = c("#ECE700","#5bc0eb","black"),    
         width=1, add = "boxplot"
         ) + rotate_x_text(45) + 
  stat_compare_means(aes(group=binaryResponse),
    #comparisons = my_comparisons,
    method="wilcox.test", # 可换其他统计方法
    #symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
    label = "p.signif"
  ) 
#immune_phenotype
immune_phenotype <- read.csv("immune_phenotype.csv",header = T,sep = ",",row.names = 1)
data <- reshape2::melt(immune_phenotype,id.vars=c("immune_phenotype"))
data <- na.omit(data)
table(data$immune_phenotype)
colnames(data) <- c("immune_phenotype","Gene","Relative expression")
ggviolin(data, x="Gene", y="Relative expression", color = "immune_phenotype", 
         ylab="Relative expression",
         xlab="Gene",
         # legend.title=x,
         add.params = list(fill="white"),
         palette = c("#ECE700","#5bc0eb","red"),    
         width=1, add = "boxplot"
) + rotate_x_text(45) + 
  stat_compare_means(aes(group=immune_phenotype),
                     #comparisons = my_comparisons,
                     method="kruskal.test", # 可换其他统计方法
                     #symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif"
  ) 
#neoantigen_burden_per_MB
neoantigen_burden_per_MB <- read.csv("neoantigen_burden_per_MB.csv",header = T,sep = ",",row.names = 1)
data <- reshape2::melt(neoantigen_burden_per_MB,id.vars=c("neoantigen_burden_per_MB"))
data <- na.omit(data)
table(data$variable)
colnames(data) <- c("neoantigen_burden_per_MB","Gene","Relative expression")
ggviolin(data, x="Gene", y="neoantigen_burden_per_MB", color = "Relative expression", 
         ylab="neoantigen_burden_per_MB",
         xlab="Gene",
         # legend.title=x,
         add.params = list(fill="white"),
         palette = c("#ECE700","#5bc0eb","red"),    
         width=1, add = "boxplot"
) + rotate_x_text(45) + 
  stat_compare_means(aes(group=neoantigen_burden_per_MB),
                     #comparisons = my_comparisons,
                     method="kruskal.test", # 可换其他统计方法
                     #symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif"
  ) 
#Lund
Lund <- read.csv("Lund.csv",header = T,sep = ",",row.names = 1)
data <- reshape2::melt(Lund,id.vars=c("Lund"))
data <- na.omit(data)
table(data$Lund)
colnames(data) <- c("Lund","Gene","Relative expression")
ggviolin(data, x="Gene", y="Relative expression", color = "Lund", 
         ylab="Relative expression",
         xlab="Gene",
         # legend.title=x,
         add.params = list(fill="white"),
         #palette = c("#ECE700","#5bc0eb","red"),    
         width=1, add = "boxplot"
) + rotate_x_text(45) + 
  stat_compare_means(aes(group=Lund),
                     #comparisons = my_comparisons,
                     method="kruskal.test", # 可换其他统计方法
                     #symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif"
  ) 
#Lund2
Lund2 <- read.csv("Lund2.csv",header = T,sep = ",",row.names = 1)
data <- reshape2::melt(Lund2,id.vars=c("Lund2"))
data <- na.omit(data)
table(data$Lund2)
colnames(data) <- c("Lund2","Gene","Relative expression")
ggviolin(data, x="Gene", y="Relative expression", color = "Lund2", 
         ylab="Relative expression",
         xlab="Gene",
         # legend.title=x,
         add.params = list(fill="white"),
         #palette = c("#ECE700","#5bc0eb","red"),    
         width=1, add = "boxplot"
) + rotate_x_text(45) + 
  stat_compare_means(aes(group=Lund2),
                     #comparisons = my_comparisons,
                     method="kruskal.test", # 可换其他统计方法
                     #symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif"
  ) 
#TCGAsubtype
TCGAsubtype <- read.csv("TCGAsubtype.csv",header = T,sep = ",",row.names = 1)
data <- reshape2::melt(TCGAsubtype,id.vars=c("TCGAsubtype"))
data <- na.omit(data)
table(data$TCGAsubtype)
colnames(data) <- c("TCGAsubtype","Gene","Relative expression")
ggviolin(data, x="Gene", y="Relative expression", color = "TCGAsubtype", 
         ylab="Relative expression",
         xlab="Gene",
         # legend.title=x,
         add.params = list(fill="white"),
         palette = c("#ECE700","#5bc0eb","red","brown"),    
         width=1, add = "boxplot"
) + rotate_x_text(45) + 
  stat_compare_means(aes(group=TCGAsubtype),
                     #comparisons = my_comparisons,
                     method="kruskal.test", # 可换其他统计方法
                     #symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif"
  ) 

##10cerna
setwd("/home/pengdu/gem/10cerna")
#lncrna和mirna
lncrna_mirna <- read.table("starBaseV3_hg38_CLIP-seq_lncRNA_all.txt",sep = "\t",header = T)
lncrna <- read.table("Training_expr.txt",header = T,row.names = 1)
lncrna <- row.names(lncrna)
lncrna_mirna <- subset(lncrna_mirna,lncrna_mirna$geneName %in% lncrna)
write.csv(lncrna_mirna,"lncrna_mirna.csv")
mirna <- lncrna_mirna$miRNAname
#mirna和mrna
#BiocManager::install("multiMiR")
library(multiMiR)
db.ver = multimir_dbInfoVersions()
mrna <- read.csv("hubgene_MMGS_mrna_brown.csv",sep = ",",header = T,row.names = 1)
mrna <- row.names(mrna)
mrna_mirna <- get_multimir(org     = 'hsa',
                         target  = mrna,
                         table   = 'validated',
                         summary = TRUE,
                         predicted.cutoff.type = 'n',
                         predicted.cutoff      = 500000)
table(mrna_mirna@data[["support_type"]])
table(mrna_mirna@data[["experiment"]])
table(mrna_mirna@data[["mature_mirna_id"]])
table(mrna_mirna@data$database)
example1_result <- mrna_mirna@data 
example1_sum <- mrna_mirna@summary               # 提取summary数据
head(example1_sum)
apply(example4@summary[, 6:10], 2, sum)
example1_Lucifer <- example1_result[grep("Luciferase", example1_result[, "experiment"]), ]
head(example1_Lucifer)

# 提取3个数据检索得到的target_symbol，绘制Veen图
library(VennDiagram)            
library(tidyverse)
db <- unique(example1_result$database)
row1 <- example1_result %>% filter(database == db[1])
row2 <- example1_result %>% filter(database == db[2])
row3 <- example1_result %>% filter(database == db[3])
venn.plot <-venn.diagram(
  x = list(target1=unique(row1$target_symbol),
           target2=unique(row2$target_symbol),
           target3=unique(row3$target_symbol)),
  filename ="Venn1.tif",
  lty ="dotted",
  lwd =0.5,
  col ="black",
  fill =c("dodgerblue", "goldenrod1", "darkorange1"),
  alpha =0.60,
  cat.col =c("dodgerblue", "goldenrod1", "darkorange1"),
  cat.cex =1,
  cat.fontface="bold",
  margin =0.05,
  cex =1)
# 获取交集
venn_list <- list(target1=unique(row1$target_symbol),
                  target2=unique(row2$target_symbol),
                  target3=unique(row3$target_symbol))
for (i in 1:length(venn_list)) {
  if(i == 1){interGenes <- venn_list[[1]]}
  else{interGenes <- intersect(interGenes, venn_list[[i]])}
}
interGenes
write.csv(example1_Lucifer,"mrna_mirna.csv")
#mirna的交集
mirna2 <- example1_Lucifer$mature_mirna_id
mirna <- intersect(mirna,mirna2)
lncrna_mirna <- subset(lncrna_mirna,lncrna_mirna$miRNAname %in% mirna)
write.csv(lncrna_mirna,"lncrna_mirna.csv")
example1_Lucifer <- subset(example1_Lucifer,example1_Lucifer$mature_mirna_id %in% mirna)
write.csv(example1_Lucifer,"mrna_mirna.csv")
#出图
library(igraph)
library(dplyr)
library(magrittr)
network_data <- read.csv("cerna.csv", header = TRUE)
# 创建空的网络对象
g <- graph.empty(n =length(c(unique(network_data$mirna),unique(network_data$lncrna),unique(network_data$mrna))), directed = TRUE)
# 添加节点
g <- g %>%
  set_vertex_attr("name", value = c(unique(network_data$lncrna), unique(network_data$mirna), unique(network_data$mrna))) %>%
  set_vertex_attr("type", value = c(rep("lncrna", length(unique(network_data$lncrna))), 
                                    rep("mirna", length(unique(network_data$mirna))), 
                                    rep("mrna", length(unique(network_data$mrna)))))
g <- set_vertex_attr(g,"color", value = ifelse(V(g)$type == "lncrna", "#fb8072", ifelse(V(g)$type == "mirna", "yellow3", "#80b1d3")))
# 添加边与边长
afedge <- c()
aflength <- c()
for(i in 1:nrow(network_data)) {
  lncrna_node <- which(V(g)$name == network_data[i,1])
  mirna_node <- which(V(g)$name == network_data[i,2])
  mrna_node <- which(V(g)$name == network_data[i,3])
  aflength <- c(aflength,20,10)
  afedge <- c(afedge,lncrna_node,mirna_node,mirna_node,mrna_node)
  
}
g <- g %>% add_edges(afedge) %>% set_edge_attr("edge.length", value = aflength)
# 添加节点大小
lncrna.size=as.vector(scale(as.vector(table(network_data$lncrna)),center = F))+15
mirna.size=as.vector(scale(as.vector(table(network_data$mirna)),center = F))+6
mrna.size=as.vector(scale(as.vector(table(network_data$mrna)),center = F))+10
V(g)$size=c(lncrna.size,mirna.size,mrna.size)
# 使用Graphopt算进行布局，保存为ceRNA.net.pdf文件
pdf(file="ceRNA.net.pdf",height=10,width=10)
plot(g, 
     layout=layout.graphopt(g),  
     vertex.label=V(g)$name,
     vertex.label.family="sans",
     vertex.label.cex=ifelse(V(g)$type == "lncRNA", 0.8, ifelse(V(g)$type == "mirna", 0.8, 0.8)),
     vertex.size=V(g)$size, 
     vertex.color=V(g)$color,
     vertex.label.color="black", 
     edge.arrow.size=0.5, 
     edge.width=0.5
)
dev.off()

##11单细胞测序
setwd("/home/pengdu/gem/11singlecell")
samples=list.files("GSE135337_RAW/",pattern = 'gz')
library(data.table)
library(limma)
ctList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)
  ct=fread(file.path("GSE135337_RAW/",pro),data.table = F,header = T)
  ct[1:4,1:4]
  ct <- as.data.frame(avereps(ct[,-1], ID = ct$Symbol)) # 取平均
  rownames(ct)=ct[,1]
  colnames(ct) = paste(gsub('.txt.gz','',pro),
                       colnames(ct) ,sep = '_')
  ct=ct[,-1] 
  return(ct)
})
#合并样品
lapply(ctList, dim)
tmp =table(unlist(lapply(ctList, rownames)))
cg = names(tmp)[tmp==length(samples)]
bigct = do.call(cbind,
                lapply(ctList,function(ct){ 
                  ct = ct[cg,] 
                  return(ct)
                }))
library(Seurat)
library(harmony)
data <- CreateSeuratObject(counts = bigct)
table(data@meta.data$orig.ident) 
#线粒体
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
#核糖体 
data[["percent.rp"]] <- PercentageFeatureSet(data, pattern = "^RP")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 1)
data <- NormalizeData(data,normalization.method = "LogNormalize",scale.factor = 1e4) 
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)
JackStrawPlot(data, dims = 1:20)
ElbowPlot(data, ndims = 50)
data <- RunHarmony(data, group.by.vars = "orig.ident")
data <- FindNeighbors(data, reduction = "harmony", dims = 1:20)
data <- RunUMAP(data, dims = 1:20, reduction = "harmony",reduction.name="umap")
data <- RunTSNE(data, dims = 1:20, reduction = "harmony",reduction.name="tsne")
DimPlot(data, reduction = "umap",group.by = "orig.ident", label = F,repel = TRUE)
DimPlot(data, reduction = "tsne",group.by = "orig.ident", label = F,repel = TRUE)
library(clustree)
sce_res <- data
for (i in c(0.05,seq(0.3,1.2,0.1))){
  sce_res <- FindClusters(sce_res,resolution = i)
}
clustree(sce_res,prefix = 'RNA_snn_res.')
rm(sce_res)
data <- FindClusters(data, resolution = 0.8)  #resolution可以设0.1-1之间，值越高，亚群数目越多，常规0.5
DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE, raster=FALSE)
DimPlot(data, reduction = "tsne", label = TRUE, repel = TRUE, raster=FALSE)
table(data$seurat_clusters)
library(dplyr)
diff.wilcox = FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "top10_diff_genes_wilcox.csv", row.names = F)
DotPlot(data, features = c("EPCAM","KRT8","KRT18","MUC1",  #尿路上皮0,1,2,3,4,5,6,8,9,10,13
                           "SPP1","CSF1R","AIF1","APOE","C1QB","CD14","CD68",  #巨噬细胞11
                           "CD2","CD3D","CD3E","CCL5","CD52","KLRB1",  #T细胞12
                           "COL1A1","PDPN","TAGLN"  #成纤维细胞7,14,15
)) + 
  RotatedAxis()
new.cluster.ids <- c("Urothelial cell","Urothelial cell","Urothelial cell","Urothelial cell","Urothelial cell",
                     "Urothelial cell","Urothelial cell","Fibroblast","Urothelial cell","Urothelial cell", 
                     "Urothelial cell","Macrophage","T cell","Urothelial cell","Fibroblast",
                     "Fibroblast","Urothelial cell")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
data$celltype <- Idents(data)
list_genes=list(Urothelial=c("EPCAM","KRT8","KRT18","MUC1"),
                Macrophage=c("SPP1","CSF1R","AIF1","APOE","C1QB","CD14","CD68"),
                T_cell=c("CD2","CD3D","CD3E","CCL5","CD52","KLRB1"),
                Fibroblast=c("COL1A1","PDPN","TAGLN"))
DotPlot(data,features=list_genes,cols = c("grey", "red"),cluster.idents = T,group.by = "celltype") + 
  RotatedAxis()+
  theme(
    # 面板
    panel.border = element_rect(color="black"), #面板边框
    panel.spacing = unit(1, "mm"), #面板间距
    
    # 分面标题
    #strip.background = element_rect(color="red"),
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', #
    
    # 坐标轴线
    axis.line = element_blank(),
  )+labs(x="", y="")
Idents(data) <- data$celltype
VlnPlot(data,features = c("LINC00930","EMX2-AS1"))+RotatedAxis()
DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE, raster=FALSE) 
DimPlot(data, reduction = "tsne", label = TRUE, repel = TRUE, raster=FALSE)
#美化dimplot
library(ggplot2)
library(ggsci)
library(grid)
library(ggrepel)
library(tidyverse)
library(tidydr)
library(patchwork)
DimPlot(data,reduction = 'umap',label = TRUE, repel = TRUE, raster=FALSE) +
  scale_colour_npg() +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  #NoLegend() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())
DimPlot(data,reduction = 'tsne',label = TRUE, repel = TRUE, raster=FALSE) +
  scale_colour_npg() +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  #NoLegend() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())
saveRDS(data,"data.rds")
#上皮细胞
table(data$celltype)
data <- subset(data, celltype == "Urothelial cell")
data <- NormalizeData(data,normalization.method = "LogNormalize",scale.factor = 1e4) 
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)
JackStrawPlot(data, dims = 1:20)
ElbowPlot(data, ndims = 50)
data <- FindNeighbors(data, dims = 1:20)
data <- RunUMAP(data, dims = 1:20, reduction.name="umap")
data <- RunTSNE(data, dims = 1:20, reduction.name="tsne")
saveRDS(data,"epithelial.rds")
#分群尝试 
library(clustree)
sce_res <- data 
for (i in c(seq(0.1,1,0.1))){
  sce_res <- FindClusters(sce_res,resolution = i)
}
clustree(sce_res,prefix = 'RNA_snn_res.')
rm(sce_res)
data <- FindClusters(data, resolution = 0.1) #resolution可以设0.1-1之间，值越高，亚群数目越多，常规0.5
metadata <- data@meta.data
sc <- as.matrix(data@assays[["RNA"]]@data)
metadata$LINC00930 <- sc["LINC00930",]
write.csv(metadata,"metadata.csv")
metadata <- read.csv("metadata.csv",header = T,row.names = 1)
data@meta.data$LINC00930.expression <- metadata$LINC00930.expression
#monocle2分析
library(monocle)
sample_ann <- data@meta.data 
# rownames(sample_ann) <- sample_ann$X 
# sample_ann <- sample_ann[,-1]
gene_ann <- data.frame(gene_short_name = rownames(data@assays$RNA),
                       row.names = rownames(data@assays$RNA))
pd <- new('AnnotatedDataFrame', data = sample_ann)
fd <- new('AnnotatedDataFrame', data = gene_ann)
subdata <- as.data.frame(data@assays$RNA@counts)
mycds <- newCellDataSet(as.matrix(data@assays$RNA@counts),
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size(),
                        lowerDetectionLimit = 1)
mycds <- detectGenes(mycds, min_expr = 0.1) 
mycds <- mycds[fData(mycds)$num_cells_expressed > 10,]
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds) 
disp_table <- dispersionTable(mycds)
unsup_clustering_genes <- subset(disp_table,mean_expression >= 0.1)
mycds <- setOrderingFilter(mycds,unsup_clustering_genes$gene_id)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds) 
saveRDS(mycds,"mycds.rds")
#以pseudotime值上色 (Pseudotime是monocle2基于细胞基因表达信息计算的概率，表示时间的先后。)
plot_cell_trajectory(mycds,color_by="Pseudotime",size=1,show_backbone=TRUE)
#以细胞类型上色
plot_cell_trajectory(mycds,color_by="LINC00930.expression",size=1,show_backbone=TRUE)
#以细胞状态上色
plot_cell_trajectory(mycds, color_by = "State", size=1,show_backbone=TRUE)
#按照seurat分群排序细胞
plot_cell_trajectory(mycds, color_by = "seurat_clusters", size=1,show_backbone=TRUE)
#时间分布图
library(ggpubr)
df <- pData(mycds)
ggplot(df, aes(Pseudotime, colour = LINC00930.expression, fill = LINC00930.expression)) +
       geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
#拟时序展示单个基因表达量
colnames(pData(mycds))
pData(mycds)$LINC00930 = log2(exprs(mycds)['LINC00930',]+1)
plot_cell_trajectory(mycds, color_by ="LINC00930") + scale_color_gsea()
library(patchwork)
#pseudotime对比
library(ggpubr)
my_comparisons <- list(c("low","high"))
ggviolin(df, x="LINC00930.expression", y="Pseudotime", color = "LINC00930.expression", 
         ylab="Pseudotime",
         xlab="LINC00930.expression",
         # legend.title=x,
         add.params = list(fill="white"),
         palette = c("#ECE700","#5bc0eb"),    
         width=1, add = "boxplot") +
  stat_compare_means(aes(group=group),comparisons = my_comparisons,
                     method="wilcox.test", # 可换其他统计方法
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif") + 
  ylim(0,45) 
#cytotrace
.libPaths(c("/home/pengdu/rpackage")) 
library(CytoTRACE)
data <- readRDS("~/gem/11singlecell/epithelial.rds")
expr <- as.matrix(data@assays$RNA@counts) 
# expr <- data@assays$RNA@counts
results <- CytoTRACE(mat = expr, 
                     #enableFast = F, 
                     ncores = 8) 
phenot <- data$LINC00930.expression  
phenot <- as.character(phenot)
names(phenot) <- rownames(data@meta.data)
emb <- data@reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phenot, emb = emb)
#cytotrace联合monocle
monocle <- readRDS("mycds.rds")
monocle_meta <- data.frame(t(monocle@reducedDimS), 
                             monocle$Pseudotime, 
                             monocle$State, 
                             monocle$LINC00930.expression)
colnames(monocle_meta) <- c("C1", "C2", "Pseudotime", "State", "LINC00930.expression")
phenot1 <- monocle$LINC00930.expression
phenot1 <- as.character(phenot1)
names(phenot1) <- rownames(monocle_meta)
emb_monocle <- monocle_meta[,1:2]
plotCytoTRACE(results, phenotype = phenot, emb = emb_monocle)
saveRDS(results, "results.rds")
#infercnvpy
library(AnnoProbe)
library(Seurat)
library(infercnv)
data <- readRDS("~/gem/11singlecell/data.rds")
dat <- as.data.frame(GetAssayData(data))
Idents(data) <- data$celltype 
groupinfo <- data.frame(v1=colnames(dat),v2=data@active.ident)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match(geneInfor[,1], rownames(dat) ),]
dim(dat)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
expFile='expFile_t.txt'## 注意这里的exp和R版不一样，做了转置
write.table(t(dat),file = expFile,sep = '\t',quote = F)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
#infercnvpy
setwd("/home/pengdu/gem/11singlecell")
infercnvpy <- read.table("infercnvpy.txt",header = T,sep = "\t")
my_comparisons <- list(c("LINC00930low","LINC00930high"))
ggviolin(infercnvpy, x="cell_groups", y="cnv_score", color = "cell_groups", 
         ylab="cnv_score",
         xlab="cell_groups",
         # legend.title=x,
         add.params = list(fill="white"),
         #palette = c("#ECE700","#5bc0eb"),    
         width=1, add = "boxplot") +
  stat_compare_means(#aes(group=group),
                     comparisons = my_comparisons,
                     #method="anova", # 可换其他统计方法
                     #symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     #label = "p.signif"
                     ) +
  coord_flip()

#GSE130001
library(dplyr)
library(Seurat)
library(patchwork)
setwd("/home/pengdu/gem/11singlecell/GSE130001")
sc.data<-Read10X_h5("BLCA_GSE130001_expression.h5")
data <- CreateSeuratObject(counts = sc.data, project = "GSE130001") 
###线粒体基因比例
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-") 
###质控可视化
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)  #top2000差异基因
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data, features = VariableFeatures(object = data))
data <- JackStraw(data, num.replicate = 100) #确定数据维度
data <- ScoreJackStraw(data, dims = 1:20)
JackStrawPlot(data, dims = 1:20)
ElbowPlot(data)
metadata <- read.table("BLCA_GSE130001_CellMetainfo_table.tsv", sep = "\t", header = TRUE, row.names = 1)
data@meta.data <- metadata
data <- FindNeighbors(data, dims = 1:8)
data <- RunUMAP(data, dims = 1:8)
data <- RunTSNE(data, dims = 1:8)
library(clustree)
sce_res <- data
for (i in c(0.05,seq(0.3,1.2,0.1))){
  sce_res <- FindClusters(sce_res,resolution = i)
}
clustree(sce_res,prefix = 'RNA_snn_res.')
rm(sce_res)
data <- FindClusters(data, resolution = 0.9) #resolution可以设0.1-1之间，值越高，亚群数目越多，常规0.5
DimPlot(data, reduction = "umap", label = TRUE) 
DimPlot(data, reduction = "tsne", label = TRUE) 
DotPlot(data, features = c("PTPRC",  #免疫细胞0,1,2,3,4,5,6,8,9,10,13
                           "VWF",  #内皮细胞11
                           "MYLK","COL1A2","MCAM","MYL9",  #成纤维细胞12
                           "EPCAM"  #上皮细胞7,14,15
)) + 
  RotatedAxis()
table(metadata$Celltype..minor.lineage.)
list_genes=list(Urothelial=c("EPCAM"),
                Immune=c("SPP1","CSF1R","AIF1","APOE","C1QB","CD14","CD68"),
                Endothelial=c("VWF"),
                Fibroblast=c("MYLK","COL1A2","MCAM","MYL9"))
DotPlot(data,features=list_genes,cols = c("grey", "red"),cluster.idents = T,group.by = "Celltype..minor.lineage.") + 
  RotatedAxis()+
  theme(
    # 面板
    panel.border = element_rect(color="black"), #面板边框
    panel.spacing = unit(1, "mm"), #面板间距
    
    # 分面标题
    #strip.background = element_rect(color="red"),
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', #
    
    # 坐标轴线
    axis.line = element_blank(),
  )+labs(x="", y="")
Idents(data) <- data@meta.data$Celltype..minor.lineage.
DimPlot(data, reduction = "umap", label = TRUE) 
DimPlot(data, reduction = "tsne", label = TRUE) 
VlnPlot(data,features = c("EMX2OS"))+RotatedAxis()
saveRDS(data,"gse130001.rds")
#美化dimplot
library(ggplot2)
library(ggsci)
library(grid)
library(ggrepel)
library(tidyverse)
library(tidydr)
library(patchwork)
data <- readRDS("~/gem/11singlecell/GSE130001/gse130001.rds")
Idents(data) <- data@meta.data$Celltype..minor.lineage.
DimPlot(data,reduction = 'umap',label = TRUE, repel = TRUE, raster=FALSE) +
  scale_colour_npg() +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  #NoLegend() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())
DimPlot(data,reduction = 'tsne',label = TRUE, repel = TRUE, raster=FALSE) +
  scale_colour_npg() +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  #NoLegend() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())
#上皮细胞
table(data$Celltype..minor.lineage.)
data <- subset(data, Celltype..minor.lineage. == "Epithelial")
data <- NormalizeData(data,normalization.method = "LogNormalize",scale.factor = 1e4) 
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)
JackStrawPlot(data, dims = 1:20)
ElbowPlot(data, ndims = 50)
data <- FindNeighbors(data, dims = 1:8)
data <- RunUMAP(data, dims = 1:8, reduction.name="umap")
data <- RunTSNE(data, dims = 1:8, reduction.name="tsne")
saveRDS(data,"epithelial.rds")
#分群尝试 
library(clustree)
sce_res <- data 
for (i in c(seq(0.1,1,0.1))){
  sce_res <- FindClusters(sce_res,resolution = i)
}
clustree(sce_res,prefix = 'RNA_snn_res.')
rm(sce_res)
data <- FindClusters(data, resolution = 0.1) #resolution可以设0.1-1之间，值越高，亚群数目越多，常规0.5
metadata <- data@meta.data
sc <- as.matrix(data@assays[["RNA"]]@data)
metadata$EMX2OS <- sc["EMX2OS",]
write.csv(metadata,"metadata.csv")
metadata <- read.csv("metadata.csv",header = T,row.names = 1)
data@meta.data$EMX2OS.expression <- metadata$EMX2OS.expression
#monocle2分析
library(monocle)
sample_ann <- data@meta.data 
# rownames(sample_ann) <- sample_ann$X 
# sample_ann <- sample_ann[,-1]
gene_ann <- data.frame(gene_short_name = rownames(data@assays$RNA),
                       row.names = rownames(data@assays$RNA))
pd <- new('AnnotatedDataFrame', data = sample_ann)
fd <- new('AnnotatedDataFrame', data = gene_ann)
subdata <- as.data.frame(data@assays$RNA@counts)
mycds <- newCellDataSet(as.matrix(data@assays$RNA@counts),
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size(),
                        lowerDetectionLimit = 1)
mycds <- detectGenes(mycds, min_expr = 0.1) 
mycds <- mycds[fData(mycds)$num_cells_expressed > 1,]
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds) 
disp_table <- dispersionTable(mycds)
unsup_clustering_genes <- subset(disp_table,mean_expression >= 0.1)
mycds <- setOrderingFilter(mycds,unsup_clustering_genes$gene_id)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds) 
saveRDS(mycds,"mycds.rds")
#以pseudotime值上色 (Pseudotime是monocle2基于细胞基因表达信息计算的概率，表示时间的先后。)
plot_cell_trajectory(mycds,color_by="Pseudotime",size=1,show_backbone=TRUE)
#以细胞类型上色
plot_cell_trajectory(mycds,color_by="EMX2OS.expression",size=1,show_backbone=TRUE)
#以细胞状态上色
plot_cell_trajectory(mycds, color_by = "State", size=1,show_backbone=TRUE)
#按照seurat分群排序细胞
plot_cell_trajectory(mycds, color_by = "seurat_clusters", size=1,show_backbone=TRUE)
#时间分布图
library(ggpubr)
df <- pData(mycds)
ggplot(df, aes(Pseudotime, colour = EMX2OS.expression, fill = EMX2OS.expression)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
#拟时序展示单个基因表达量
colnames(pData(mycds))
pData(mycds)$EMX2OS = log2(exprs(mycds)['EMX2OS',]+1)
plot_cell_trajectory(mycds, color_by ="EMX2OS") + scale_color_gsea()
library(patchwork)
#pseudotime对比
library(ggpubr)
my_comparisons <- list(c("low","high"))
ggviolin(df, x="EMX2OS.expression", y="Pseudotime", color = "EMX2OS.expression", 
         ylab="Pseudotime",
         xlab="EMX2OS.expression",
         # legend.title=x,
         add.params = list(fill="white"),
         palette = c("#ECE700","#5bc0eb"),    
         width=1, add = "boxplot") +
  stat_compare_means(aes(group=group),comparisons = my_comparisons,
                     method="wilcox.test", # 可换其他统计方法
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif") + 
  ylim(0,45) 
#cytotrace
.libPaths(c("/home/pengdu/rpackage")) 
library(CytoTRACE)
data <- readRDS("~/gem/11singlecell/GSE130001/epithelial.rds")
expr <- as.matrix(data@assays$RNA@counts) 
# expr <- data@assays$RNA@counts
results <- CytoTRACE(mat = expr, 
                     #enableFast = F, 
                     ncores = 8) 
phenot <- data$EMX2OS.expression   
phenot <- as.character(phenot)
names(phenot) <- rownames(data@meta.data)
emb <- data@reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phenot, emb = emb)
#cytotrace联合monocle
monocle <- readRDS("mycds.rds")
monocle_meta <- data.frame(t(monocle@reducedDimS), 
                           monocle$Pseudotime, 
                           monocle$State, 
                           monocle$EMX2OS.expression)
colnames(monocle_meta) <- c("C1", "C2", "Pseudotime", "State", "EMX2OS.expression")
phenot1 <- monocle$EMX2OS.expression
phenot1 <- as.character(phenot1)
names(phenot1) <- rownames(monocle_meta)
emb_monocle <- monocle_meta[,1:2]
plotCytoTRACE(results, phenotype = phenot, emb = emb_monocle)
saveRDS(results, "results.rds")
#infercnvpy
library(AnnoProbe)
library(Seurat)
library(infercnv)
data <- readRDS("~/gem/11singlecell/GSE130001/gse130001.rds")
dat <- as.data.frame(GetAssayData(data))
Idents(data) <- data$Celltype..minor.lineage. 
groupinfo <- data.frame(v1=colnames(dat),v2=data@active.ident)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match(geneInfor[,1], rownames(dat) ),]
dim(dat)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
expFile='expFile_t.txt'## 注意这里的exp和R版不一样，做了转置
write.table(t(dat),file = expFile,sep = '\t',quote = F)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
table(data$Celltype..minor.lineage.)
#infercnvpy
setwd("/home/pengdu/gem/11singlecell/GSE130001")
infercnvpy <- read.table("infercnvpy.txt",header = T,sep = "\t")
my_comparisons <- list(c("EMX2OSlow","EMX2OShigh"))
ggviolin(infercnvpy, x="cell_groups", y="cnv_score", color = "cell_groups", 
         ylab="cnv_score",
         xlab="cell_groups",
         # legend.title=x,
         add.params = list(fill="white"),
         #palette = c("#ECE700","#5bc0eb"),    
         width=1, add = "boxplot") +
  stat_compare_means(#aes(group=group),
    comparisons = my_comparisons,
    #method="anova", # 可换其他统计方法
    symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
    label = "p.signif"
  ) +
  coord_flip()

##12pcr
setwd("/home/pengdu/gem/12pcr")
library(ggplot2)
library(ggpubr)
pcr <- read.csv("sample_pcr.csv",header = T,sep = ",")
ggplot(pcr, aes(x= group,y= expression,color= group)) +
  stat_boxplot(geom='errorbar',width=0.3,
               position=position_dodge(0.75))+
  geom_boxplot(position=position_dodge(0.75)) +
  # geom_jitter(mapping=aes(color = group),width = 0.05)+
  scale_color_manual(values = c("#87CEEB","#800000")) +
  geom_point(mapping=aes(color=group))+
  labs(x='Group',y='Relative expression')+
  guides(color='none')+
  stat_compare_means(
    comparisons=list(c("Normal","Tumor")),
    method="t.test",#t.test or wilcox.test
    paired = T,
    hide.ns = F,
    label = 'p.signif' #"p.signif" (shows the significance levels), "p.format" (shows the formatted p value)
  ) + theme_bw() + 
  geom_line(aes(group = paired),
                color='grey',size=0.5,
                position=position_dodge(0.01)
  ) + 
  facet_wrap(.~lncrna,nrow=1)
#生存曲线
bmt <-read.csv('survival_pcr.csv',header=T,sep = ",")
library(survminer)
#OS
fit1<-survfit(Surv(OS,status == 1)~LINC00930,bmt)
print(fit1)
ggsurvplot(fit1,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="BJCANCER cohort",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)
fit2<-survfit(Surv(OS,status == 1)~EMX2OS,bmt)
print(fit2)
ggsurvplot(fit2,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv",
           xlab="Time (days)",ylab="Overall survival",ncensor.plot = TRUE,
           legend=c(0.8,0.75),legend.title="BJCANCER cohort",ggtheme=theme_light(),xlim=c(0,1825),break.x.by=365)

