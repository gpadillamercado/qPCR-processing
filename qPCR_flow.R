library(tidyverse)
library(readxl)

# Loads in the Dataset
qPCRm <- read_excel("qPCR_testData.xls", sheet = "Results", skip = 7,
  .name_repair = "unique", na = c("Undetermined", "NTC"))
# Look into .name_repair

# Cleans up the subscripts for Ct values
colnames(qPCRm) <- colnames(qPCRm) %>%
  str_replace_all(pattern = "\u0442", replacement = "t") %>%
  str_replace_all(pattern = "\u0394", replacement = "d")

# Separates Biological Replicate number from Sample Name
qPCRm <- qPCRm %>%
  select(Well, `Sample Name`, `Target Name`, Ct) %>%
  separate(`Sample Name`, sep = "(\\sOXP\\s|\\s)", into = c("Sample", "Number")) %>%
  filter(!is.na(Number))

# Ensures each column has correct atomic vectors
qPCRm$Number <- qPCRm$Number %>% as.numeric()
qPCRm$Ct <- qPCRm$Ct %>% as.double()
names(qPCRm)[4] <- "Target"
qPCRm$Target <- qPCRm$Target %>% toupper()


# Test dataset for simple test
# Intend to upload a minimal but sufficiently complex test dataset
df <- tibble(
  Sample = c( rep("WT", 6), rep("Mut", 6) ),
  Number = rep( c("1" , "2"), 6),
  Target = rep( c("Ref", "Exp1", "Exp1", "Exp2", "Ref", "Exp2"), 2),
  Ct = c( 
    rnorm(2, 18, 0.2),
    rnorm(2, 25, 0.1),
    rnorm(2, 20, 0.1), 
    rnorm(2, 17, 0.2), 
    rnorm(2, 21, 0.1), 
    rnorm(2, 24, 0.2) )
)

# This defines which two samples to compare in t-test
ctrlGene <- "18S"
ctrlSample <- "WT"

# Averages technical replicates
df <- qPCRm %>% 
  select(-Well) %>% 
  group_by(Sample, Number, Target) %>%
  summarise(Ct = mean(Ct)) %>%
  as_tibble() %>% ungroup()

# This groups df by Target
# Then this mutates the grouped df by
# taking values in the Ct column and subtracting them by
# the Ct of Targets that are ctrlGene in that group
df <- df %>%
  group_by(Target) %>%
  mutate(dCt = Ct - filter(df, Target == ctrlGene)$Ct)

# Gets a dataframe with the average among control samples of dCt by Target genes
mean_dCt <- df %>% 
  filter(Sample == ctrlSample & !(Target == ctrlGene)) %>%
  group_by(Target) %>%
  summarize(Average = mean(dCt))

# Excludes the ctrlGene from data
# Adds a column that subtracts dCt's from mean_dCt in matching Target genes
# by biological replicate (Number)
df <- df %>% ungroup() %>%
  group_by(Number) %>%
  filter(!(Target == ctrlGene)) %>%
  mutate(ddCt = dCt - filter(mean_dCt, Target %in% df$Target)$Average) %>%
  ungroup() %>%
  mutate(FoldChange = 2 ^ ddCt)


# Performs a t.test on Sample1 and Sample2
# Returns pvalue, degrees of freedom, variances of both samples, and t.test method
# Updated from summarize_each, that function is deprecated
Sample1 <- "WT"
Sample2 <- "Mut"

df %>%
  group_by(Target) %>%
    summarise_at(
      vars(ddCt), list(
        pvalue = ~ t.test(.[Sample == Sample1], .[Sample == Sample2], var.equal = TRUE)$p.value,
        dfree  = ~ t.test(.[Sample == Sample1], .[Sample == Sample2], var.equal = TRUE)$parameter,
        vari1  = ~ sd(.[Sample == Sample1]) ^ 2,
        vari2  = ~ sd(.[Sample == Sample2]) ^ 2,
        test   = ~ t.test(.[Sample == Sample1], .[Sample == Sample2], var.equal = TRUE)$alternative)
    )
