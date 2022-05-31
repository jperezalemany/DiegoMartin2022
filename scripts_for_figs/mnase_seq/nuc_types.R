
library(ggplot2)
library(tidyr)
library(dplyr)

nuc_types <- rbind(
  read.table("processed_tables/fuzz.txt", header = T, sep = "\t") %>%
    mutate(change = "fuzziness"),
  read.table("processed_tables/occ.txt", header = T, sep = "\t") %>%
    mutate(change = "occupancy"),
  read.table("processed_tables/shift.txt", header = T, sep = "\t") %>%
    mutate(change = "shift")
)

nuc_types
nuc_types$change <- factor(
  c(rep("Fuzziness", 4), rep("Occupancy", 4), rep("Position shift", 4)),
  levels = c("Occupancy", "Fuzziness", "Position shift")
)
nuc_types[nuc_types$change == "Position shift" & nuc_types$mutant == "tpf", "from"] <- 465353
nuc_types[nuc_types$change == "Position shift" & nuc_types$mutant == "chr", "from"] <- 471570
nuc_types %>%
  mutate(group2 = paste(mutant, change)) %>%
  group_by(group2) %>%
  summarise(n = sum(number), group = change,
            from = max(from), mutant = mutant) %>%
  mutate(percent = n / from * 100) %>%
  ggplot(aes(group, percent, fill = mutant)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0, 0.05, 0.05), limits = c(0, 30)) +
  ylab("Percentage of dyads") +
  scale_fill_manual(values = brewer.pal(3, "Set2")[2:3]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))

ggsave(filename = "figures/categories/barplot.pdf", height = 6, width = 6, units = "cm", dpi = 1200)


