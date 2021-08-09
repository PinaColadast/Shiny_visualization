library(inspectdf)
library(ggplot2)
library(dplyr)
library(colorspace)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(reshape2)
library(comprehenr)

Sys.setenv(language="en")
data.anno <-  read.table( "C:/Users/jtao/work_dir/RShiny/test_with_AML_data/data/AML_GPL570_annotation_subset.txt",  
                     header = TRUE, sep = "\t",check.names = FALSE, row.names = 1)
data <- read.table("C:/Users/jtao/work_dir/RShiny/test_with_AML_data/data/AML_GPL570_matrix_subset.txt", 
                        header = TRUE, sep = "\t", check.names = FALSE, row.names=1)
data.anno$sample_ID <- rownames(data.anno) 
data.anno <- data.anno[order(factor(data.anno$sample_ID, levels = colnames(data))), ]

table(data.anno$Indication)

mutation <- grep("mutation", colnames(data.anno), value = TRUE)

data_mutation <- data.anno[, mutation]
# data_plot <- melt(data_mutation, 
table(data_mutation[, c(3)])


ggboxplot(data[c("A1BG", "A1CF"), ])

data_gene <- data[data$"Gene Symbol" %in% c("A1BG", "A1CF"),]
dd <- data[data$"Gene Symbol" == "A1BG", ]
dataaa <- melt(data_gene, id.vars = "Gene Symbol")

dataaa$Sample_type <- data.anno$Sample_type
p1 <- ggboxplot(data_melt, x = "variable", y ="value", color = "sample type")
p1


# data_gene <- data_bind[rownames(data_bind) %in% c("A1BG", "A1CF", "Sample type"), ]
data_gene <- data[rownames(data) %in% c("CD33", "IL3RA", "CD70"), ]
a <- as.data.frame(t(data_gene))
str(a)
a$`sample type` <- c(data.anno$Sample_type)
# data_gene$"Gene Symbol" <- rownames(data_gene)
# a$id <- rownames(a)
# dataaa <- melt(data_gene)
data_melt <- melt(a, id.vars = c("sample type"))
data_melt <- data_melt %>% filter(variable %in% c("CD33", "IL3RA", "CD70")) %>% droplevels
all.equal(c(rownames(data.anno), rownames(data.anno)), c(data_melt$variable))
str(data_melt)


table(data.anno %>% filter(Indication_short == "AML") %>% select(Sample_type))
table(data.anno %>% filter(Sample_type == "normal control") %>% select(Indication_short))
table(data.anno$Indication_short)


data.anno[is.na(data.anno)] <- "no information" 
# data.anno$mutation_status <-

mutation_cate <- grep("mutation", colnames(data.anno), value =TRUE)



return_mutation<- function(x1, x2){
  # print(length(x1))
  if (is.na(x1) | is.na(x2)){
    mut = "no info"
  }
  else if (x1==x2){
    mut = x1}

  else if (x1!=x2){
    mut = paste(x1, ", ", x2)
  }
  return(mut)
}
  


data.anno$mutation_status1 <- apply(data.anno[, mutation_cate],
                                   1,
                            FUN = function(x) min(names(which(x=="positive",
                                                                          ))))

data.anno$mutation_status2 <- apply(data.anno[, mutation_cate],
                                    1,
                                    FUN = function(x) max(names(which(x=="positive",
                                    ))))

# data.anno %>% mutate(mutation_status== return_mutation(mutation_status1, mutation_status2))
  # apply(data.anno, 1, FUN=return_mutation, data.anno["mutation_status1"],
                                   # data.anno["mutation_status2"])

data.anno$mutation_status <- apply(data.anno[, c("mutation_status1",
                                               "mutation_status2")], 1,
                                   function(y) return_mutation(y[1], y[2]))
table(data.anno$mutation_status)
# mutation <- c()
# for i
                             
data.anno$mutation_status <- to_vec(for(i in (1:dim(data.anno)[1])) if (is.na(data.anno$mutation_status1[i])==FALSE) return_mutation(data.anno$mutation_status1[i],
                                                                                                                                     data.anno$mutation_status2[i]) else "no info")
data.anno
table(data.anno$mutation_status)
