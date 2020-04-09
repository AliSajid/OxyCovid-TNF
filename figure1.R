library(tidyverse)

df <- read_csv("results/tnf.csv") %>% 
  filter(cellline != "MCF7")

p <-  ggplot(df, 
             aes(x = perturbagen, 
                 y = similarity, 
                 group = cellline, 
                 color = cellline,
                 shape  = cellline,
                 label = similarity)
             )

p + geom_point() + geom_line() + 
  scale_x_discrete(limits = c("Carbetocin", "Chloroquine")) +
  xlab("Medication") +
  ylab("Similarity Score relative to TNF Knockout") +
  scale_color_discrete(name = "Cell Line") + 
  scale_shape_manual(
    name = element_blank(),
    values = c(16, 16)
  ) +
  scale_y_continuous(
    breaks = seq(0.2, 0.8, 0.1)
  ) +
  geom_label(hjust = -0.5, vjust = 0) +
  guides(size=FALSE,
         color = guide_legend(override.aes = list(shape=16)),
         shape = FALSE
  )

