library(tidyverse)

all_pairwise <- read_csv("/data/classes/2024/fall/biol343/course_files/all_pairwise_comparisons.csv")

all_counts <- read_csv("/data/classes/2024/fall/biol343/course_files/counts_clean.csv")

all_pairwise

all_counts

# tidy data

counts_long <- all_counts |>
    pivot_longer(-gene_id, names_to = 'sample', values_to = 'count')

counts_long |> 
    ggplot() +
    geom_col(aes(x = sample, y = count))

gene_of_interest <- counts_long |>
    filter(gene_id == 'Smp_315690') |>
    separate_wider_delim(cols = sample, delim = "_", names = c('tissue', 'misc')) |>
    separate_wider_position(cols = misc, widths = c('age' = 2, 'rep' = 1))

gene_of_interest |>
    ggplot(aes(x = tissue, y = count)) +
    # geom_boxplot() +
    geom_point(aes(shape = age, color = rep), size = 7)

ugly_plot <- counts_long |>
    filter(count > 0) |>
    separate_wider_delim(cols = sample, delim = "_", names = c('tissue', 'misc')) |>
    separate_wider_position(cols = misc, widths = c('age' = 2, 'rep' = 1)) |>
    ggplot() +
    geom_density(aes(x = count, color = tissue)) +
    lims(x = c(0, 100))
ugly_plot
