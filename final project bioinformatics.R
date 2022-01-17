#Zackery Wood
#Bioinformatics Final Project
#a gene expression analysis of the data set GSE43261
library(limma)
library(GEOquery)
library(ggplot2)
library(dplyr)
library(affy)

GSE43261 <- getGEO("GSE43261") 
View(GSE43261)

GSE43261


#38 samples
#45101 features


GSE43261.expr <- exprs(GSE43261[[1]])
GSE43261.p <- pData(GSE43261[[1]])
boxplot(GSE43261.expr, main = "processed data")

GSE43261.expr <- log2(GSE43261.expr) 
boxplot(GSE43261.expr, main = "log2 processed data")
GSE43261.p$`response:ch1`
status = as.character(GSE43261.p$`response:ch1`)

table(status)

#there's 16 control, 8 resistant, and 14 responder

design = model.matrix(~0+status)

#here i removed the control column along with the first 16 data points 
#because they were all control samples

colnames(design) <- c("Control","Resistant", "Responder")
design


dim(GSE43261.expr)
dim(design)

fit <- lmFit(GSE43261.expr, design)
head(fit$coefficients)
head(fit$sigma)

contrast.matrix <- makeContrasts(Resistant - Responder,levels=design)
fit1 <- contrasts.fit(fit, contrast.matrix)
head(fit1$coefficients)

fit1 <- eBayes(fit1)
tt <- topTable(fit1,sort.by = "p")
tt
fit1
probe <- rownames(tt)[1]
m <- match(probe, rownames(GSE43261.expr))
m

df <- data.frame(expr =GSE43261.expr[m,], status = status)

means <- df %>% group_by(status) %>% summarise(mean = mean(expr))
means
diff(means$mean) 
head(tt)
logFC <- tt[3,]$logFC
2**logFC

FC <- paste0("FC = ", round(2**logFC, 2))
main <- paste0("Expression of ", probe, ", ", FC)

ggplot(df, aes(x = status, y = expr, fill = status)) + geom_boxplot() +
  ylab("log2 expression") + ggtitle(main) +
  scale_fill_manual(values = c("pink", "lightblue", "orange")) +
  theme_classic() + theme(legend.position = "none")

tt.20 <- topTable(fit1,sort.by = "p", p.value = 0.20, number = nrow(GSE43261.expr))
nrow(tt.20)
top.3 = head(tt.20)
top.3 = top.3[-c(4,5,6),]
head(top.3, 3)
top.3 = top.3[,-c(2,3,4,6)]
top.3

top_probe = top.3[1,]
top_probe

probe2 = rownames(top.3)[1]
try1 = match(probe2, rownames(GSE43261.expr))
try1
df2 <- data.frame(expr =GSE43261.expr[try1,], status = status)
View(df2)
probe2
ggplot(df2, aes(x = status, y = expr, fill = status)) + geom_boxplot() +
  ylab("log2 expression") + ggtitle(main) +
  scale_fill_manual(values = c("pink", "lightblue", "orange")) +
  theme_classic() + theme(legend.position = "none")

#9.
platform <- annotation(GSE43261[[1 ]])  
pl <- getGEO(platform)
pl <- Table(pl)
probe3 = rownames(tt.20)
try2 = match(probe3, pl$ID)
df3 <- data.frame(pl[try2,])
View(df3)
tt.20 = tt.20[,-c(2,3,4,6)]

numbereleven = cbind(tt.20,df3$Gene.Title)
colnames(numbereleven)= c("logFC","adj.P.val","Gene Name")
View(numbereleven)
head(numbereleven, 5)
probes = numbereleven$`Gene Name`
keep <- probes!=""
genes <- probes[keep]
genes <- strsplit(genes, " ///" )
genes <- unlist(genes)

# get a unique set of genes 
genes <- unique(genes) 
write.txt(genes, row.names = FALSE, quote = FALSE, "filename.csv")


####
# Summary: In this analysis, GSE43261 is being analyzed. This data set is titled Fluoxetine
#resistance in mice is associated with attenuated
#progression of a stereotyped dentate gyrus gene expression program. This gene is  
#the most common treatment for major depression.
#The number of samples in this data set is 38 samples and 45101 features. 
#There are 3 groups 16 control, 8 resistant, 14 responder. There were 1886
#differentially expressed probes that corresponded to the FDR of .20.  
# The name of the probes are 1460092_at (no name), 1456121_at (centrosomal protein 97), and 1425608_at 
#(dual specificity phosphatase 3 (vaccinia virus phosphatase VH1-related)).
#
#
#
#
#
#
#
#
#
#
#



