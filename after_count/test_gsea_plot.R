library(fgsea)
library(ggplot2)
data(examplePathways)
data(exampleRanks)
p <- plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]], exampleRanks)
str(p$layers)
p$layers[[1]]$aes_params$linewidth <- 3
p$layers[[1]]$aes_params$size <- 3
ggsave("test.pdf", plot = p)
print("SUCCESS")
