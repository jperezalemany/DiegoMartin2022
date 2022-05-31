
setwd("~/../Dropbox/Bioinformatics/MNase-seq tpf chr/")

library(stringr)

library(ggplot2)

insert_sizes <- data.frame()
names <- c(
  "insert_size", "pairs_total", "inward_oriented_pairs", 
  "outward oriented pairs", "other pairs"
)

for (file in dir("samtools/", pattern = "is")) {
  sample <- file %>% str_remove(".clean.is.txt")
  condition <- sample %>% str_remove("_rep[0-9]")
  rep <- sample %>% str_remove(paste0(condition, "_"))
  inserts <- read.table(file.path("samtools/", file), sep = "\t", col.names = names)
  inserts$sample <- sample
  inserts$condition <- condition
  inserts$rep <- rep
  
  insert_sizes <- rbind(insert_sizes, inserts)
}

table(insert_sizes$sample)

for (c in unique(insert_sizes$condition)) {
  for (r in unique(insert_sizes$rep)) {
    max_point <- insert_sizes %>%
      filter(condition == c, rep == r) %>%
      select(pairs_total) %>%
      unlist() %>% max()
    max_is <- insert_sizes %>%
      filter(condition == c, rep == r, pairs_total == max_point) %>%
      select(insert_size)
    print(c)
    print(r)
    print(max_is)
  }
}
#color = brewer.pal(9, "Blues")[9], size = 0.75
insert_sizes %>%
  mutate(condition = factor(condition, levels = c("col0", "tpf", "chr"))) %>%
  mutate(rep = factor(rep, levels = c("rep1", "rep2"))) %>%
  ggplot(aes(insert_size, pairs_total / 10**6)) +
  geom_line(aes(color = rep), size = 0.75) +
  scale_color_manual(values = brewer.pal(3, "Dark2")[2:3]) +
  facet_grid(~condition) +
  xlab("Insert size (pb)") +
  ylab("Read count (M)") +
  theme_bw() +
  scale_x_continuous(limits = c(50, 200)) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))

dir.create("figures/supp")
ggsave(filename = "figures/supp/insert_sizes.pdf", height = 5, width = 15, units = "cm", dpi = 1200)
