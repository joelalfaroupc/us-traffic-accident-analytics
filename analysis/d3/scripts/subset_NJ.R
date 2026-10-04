library(data.table)
library(dplyr)

input_path <- 'Datasets/US_Accidents_March23.csv'
output_path <- 'Datasets/US_Accidents_Final.csv'

df <- fread(input_path)

df_nj <- df %>% filter(State == 'NJ')
if (nrow(df_nj) > 9020) {
  set.seed(42)
  df_final <- df_nj %>% slice_sample(n = 9020)
} else {
  df_final <- df_nj
}

fwrite(df_final, output_path)

print(paste("Dimensiones finales:", nrow(df_final), ",", ncol(df_final)))