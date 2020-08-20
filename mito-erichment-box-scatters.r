library(tidyverse)
library(ggplot2)
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("yaml")
install.packages("stringi")

version


setwd("C:/Users/joseph7e/Documents/mito-enrichment-data/")
results <- read.csv("mito-enrichment-results-with-gc-bins.csv", header = T, sep = ',')

head(results)




ggplot(data = results, mapping = aes(x = phylum2, y = CPM_original)) +
    geom_boxplot(alpha = 0) +
    geom_jitter(alpha = 0.3, color = "tomato") + ylab("CPM") + ggtitle("Mitochondrial Genomes Coverage Per Million Reads - Standard NGS Library") + xlab("Phylum")

ggplot(data = results, mapping = aes(x = phylum2, y = CPM_reverse)) +
    geom_boxplot(alpha = 0) +
    geom_jitter(alpha = 0.3, color = "tomato") + ylab("CPM") + ggtitle("Mitochondrial Genomes Coverage Per Million Reads - Enriched NGS Library") + xlab("Phylum")

ggplot(data = results, mapping = aes(x = phylum, y = CPM_enriched_reads)) +
    geom_boxplot(alpha = 0) +
    geom_jitter(alpha = 0.3, color = "tomato") + ylab("CPM") + ggtitle("Enriched Mitochondrial Genomes - Coverage Per Million Reads")


ggplot(data = results, mapping = aes(x = phylum2, y = factor)) +
    geom_boxplot(alpha = 0) +
    geom_jitter(alpha = 0.3, color = "tomato") + ylab("Enrichment factor") + ggtitle("Mitochondrial Genome Enrichment by Phylum") + xlab("Phylum")

results$Bin <- factor(results$Bin,levels = c("one", "two", "three", "four"))

ggplot(data = results, mapping = aes(x = GC, y = factor)) +
    geom_boxplot(alpha = 0) +
    geom_jitter(alpha = 0.3, color = "tomato") + ylab("factor") + ggtitle("Mitochondrial Genomes - Enrichment Factor (original CPM / enriched CPM)")
   



ggplot(data = results, mapping = aes(x = Bin, y = CPM_original)) +
    geom_boxplot(alpha = 0) +
    geom_jitter(alpha = 0.3, color = "tomato") + ylab("factor") + ggtitle("Mitochondrial Genomes - GC bin VS Enrichment Factor (original CPM / enriched CPM)")


ggplot(data = results, mapping = aes(x = phylum, y = GC)) +
    geom_boxplot(alpha = 0) +
    geom_jitter(alpha = 0.3, color = "tomato") + ylab("GC_percentage") + ggtitle("Mitochondrial Genomes - Enrichment Factor (original CPM / enriched CPM)")



sresults <- read.csv("mito-enrichment-results-with-gc-bins-no-all.csv", header = T, sep = ',')


ggplot(sresults, aes(x=GC, y=CPM_reverse)) + ylab("Average coverage per million reads") + ggtitle("Enriched Library CPM vs GC content") + xlab("GC content") +  geom_point() + 
  geom_smooth(method=lm) +geom_rug()

ggplot(sresults, aes(x=GC, y=CPM_original)) + ylab("Average coverage per million reads") + ggtitle("Standard Library CPM vs GC content") + xlab("GC content") +  geom_point() + 
  geom_smooth(method=lm) +geom_rug()

ggplot(sresults, aes(x=GC, y=CPM_original, color=phylum2)) + geom_point()  + ylab("Average coverage per million reads") + ggtitle("Mitochondrial Genome Coverage vs GC content") + xlab("GC content")



ggplot(sresults, aes(x=GC, y=CPM_original_reads, color=phylum2) + geom_point()  + ylab("Average coverage per million reads") + ggtitle("Mitochondrial Genome Coverage vs GC content") + xlab("GC content") + geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95)


library(ggplot2)
install.packages("stringi")
library(tidyverse)
setwd("C:/Users/joseph7e/Documents/mito-enrichment-data/")
results_basic <- read.csv("simplified-enrichment.csv", header = T, sep = ',')

head(results_basic)
#scale_x_discrete(limits = rev(levels(the_factor)))




ggplot(data = results_basic, mapping = aes(x = fct_reorder(Library,CPM), y = log(CPM))) +
    geom_boxplot(alpha = 0) +geom_jitter(alpha = 0.3, color = "red") + xlab('Library Prep Method') + ylab("log(CPM)") + ggtitle("Mitochondrial Genome Enrichment" ) 


ggplot(data = results_basic, mapping = aes(x = Library, y = log(CPM))) +
    geom_boxplot(alpha = 0) +
    geom_jitter(alpha = 0.3, color = "red") + xlab('Library Prep Method') + ylab("log(CPM)") + ggtitle("Mitochondrial Genome Enrichment" ) 



