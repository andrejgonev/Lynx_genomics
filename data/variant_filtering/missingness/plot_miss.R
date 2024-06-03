library(ggplot2)
library(dplyr)

# read the missingness table
miss_table <- read.table("miss_table.txt",
                         header = TRUE)
miss_table$percent_out <- round(miss_table$filtered_snps/miss_table$total_snps*100, 2)
miss_table$percent_retained <- round((1 - miss_table$filtered_snps / miss_table$total_snps) * 100, 2)

# populations=("bal" "cau" "crp" "lva" "nor" "pol" "tva" "mng" "wru" "yak" "pyk")

miss_plot <- ggplot(miss_table) +
    geom_line(data = miss_table %>% filter(population == "bal"),
                        aes(x = f_miss, y = percent_retained, color = population), size = 0.5) +
    geom_line(data = miss_table %>% filter(population == "cau"),
                        aes(x = f_miss, y = percent_retained, color = population), size = 0.5) +
    geom_line(data = miss_table %>% filter(population == "crp"),
                        aes(x = f_miss, y = percent_retained, color = population), size = 0.5) +
    geom_line(data = miss_table %>% filter(population == "lva"),
                        aes(x = f_miss, y = percent_retained, color = population), size = 0.5) +
    geom_line(data = miss_table %>% filter(population == "nor"),
                        aes(x = f_miss, y = percent_retained, color = population), size = 0.5) +
    geom_line(data = miss_table %>% filter(population == "pol"),
                        aes(x = f_miss, y = percent_retained, color = population), size = 0.5) +
    geom_line(data = miss_table %>% filter(population == "tva"),
                        aes(x = f_miss, y = percent_retained, color = population), size = 0.5) +
    geom_line(data = miss_table %>% filter(population == "mng"),
                        aes(x = f_miss, y = percent_retained, color = population), size = 0.5) +
    geom_line(data = miss_table %>% filter(population == "wru"),
                        aes(x = f_miss, y = percent_retained, color = population), size = 0.5) +
    geom_line(data = miss_table %>% filter(population == "yak"),
                        aes(x = f_miss, y = percent_retained, color = population), size = 0.5) +
    geom_line(data = miss_table %>% filter(population == "pyk"),
                        aes(x = f_miss, y = percent_retained, color = population), size = 0.5) +
    geom_point(data = miss_table, 
                         aes(x = f_miss, y = percent_retained, fill = population),
                         shape = 21, size = 1.5) +
    theme(panel.grid.major.y = element_line(color = "black", size = 0.2, linetype = 2),
          panel.background = element_rect(fill = "white"),
          panel.grid.major.x = element_line(color = "lightgrey", size = 0.2, linetype = 2)) +
    scale_y_continuous(n.breaks = 12) +
    xlab(paste0("Maximum proportion of missing genotypes")) +
    ylab(paste0("Percentage of SNPs retained"))

ggsave(filename = "missing_plot.png",
       plot = miss_plot, width = 4, height = 3)
