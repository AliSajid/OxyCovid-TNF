library(tidyverse)

generate_figure <- function(gene, cell_lines, section = "inflammation") {
  selected <- cell_lines
  prefix <- "results/"
  file <- paste(gene, "csv", sep = ".")
  
  df <- read_csv(paste(prefix, file, sep = "/")) %>% 
    filter(cellline %in% selected)
  
  p <-  ggplot(df, 
               aes(x = perturbagen, 
                   y = similarity, 
                   group = cellline, 
                   color = cellline,
                   shape  = cellline,
                   label = similarity)
  )
  
  bp <- p + geom_point() + geom_line() + 
    scale_x_discrete(limits = c("Bupropion", "Chloroquine", "Carbetocin")) +
    xlab("Medication") +
    ylab(paste("Similarity Score relative to", gene, "Knockout")) +
    scale_color_discrete(name = "Cell Line") + 
    scale_shape_manual(
      name = element_blank(),
      values = c(16, 16)
    ) +
    geom_label(hjust = -0.5, vjust = 0) +
    guides(size=FALSE,
           color = FALSE,
           shape = FALSE
    )
  
  bp
  figure_file <- paste(paste(section, "comparison", gene, sep = "-"),
                       c("png", "pdf"), sep = ".")
  figure_prefix = "figures"
  fig_path <- paste(figure_prefix, figure_file, sep = "/")
  ggsave(fig_path[1], device = "png")
  ggsave(fig_path[2], device = "pdf")
  
}

generate_figure("TNF", c("HA1E"))
generate_figure("ARG1", c("HA1E"), section = "immune")
generate_figure("CD83", c("HA1E"), section = "immune")
generate_figure("TLR7", c("HA1E"), section = "immune")
generate_figure("TLR9", c("HA1E"), section = "immune")


# selected <- c("HA1E", "PC3")
# 
# df <- read_csv("results/tnf.csv") %>% 
#   filter(cellline %in% selected)
# 
# p <-  ggplot(df, 
#              aes(x = perturbagen, 
#                  y = similarity, 
#                  group = cellline, 
#                  color = cellline,
#                  shape  = cellline,
#                  label = similarity)
#              )
# 
# bp <- p + geom_point() + geom_line() + 
#   scale_x_discrete(limits = c("Bupropion", "Carbetocin", "Chloroquine")) +
#   xlab("Medication") +
#   ylab("Similarity Score relative to TNF Knockout") +
#   scale_color_discrete(name = "Cell Line") + 
#   scale_shape_manual(
#     name = element_blank(),
#     values = c(16, 16)
#   ) +
#   geom_label(hjust = -0.5, vjust = 0) +
#   guides(size=FALSE,
#          color = guide_legend(override.aes = list(shape=16)),
#          shape = FALSE
#   )
# 
# bp
# 
# ggsave("figures/inflammation-comparsion.png", device = "png")
# ggsave("figures/inflammation-comparsion.pdf", device = "pdf")
