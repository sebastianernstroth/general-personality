library(rprojroot)
library(tidyverse)
library(jmv)
library(rstatix)

# Establish project directory
root <- rprojroot::is_rstudio_project

# Set data directory
dataDir <- root$find_file("data/IPIP-FFM-data-8Nov2018")

# Grab data files
csvFiles <- list.files(dataDir, pattern = "*.csv", recursive = TRUE)

# Read data
data <- read.csv(paste(dataDir, csvFiles[1], sep = "/"), sep = "\t", na.strings = c("0", "NULL", "NONE"))

# Summary statistics of individual items
data %>%
  select(EXT1:OPN10) %>%
  pivot_longer(cols = c(EXT1:OPN10), names_to = "item", values_to = "value") %>%
  group_by(item) %>%
  get_summary_stats(type = "mean_sd") %>%
  print(n = 50)
  


df <- data %>%
  # filter_at(vars(EXT1:OPN10), all_vars(!is.na(.))) %>%
  mutate(across(c(EXT2,
                  EXT4,
                  EXT6,
                  EXT8,
                  EXT10,
                  EST1,
                  EST3,
                  EST5,
                  EST6,
                  EST7,
                  EST8,
                  EST9,
                  EST10,
                  AGR1,
                  AGR3,
                  AGR5,
                  AGR7,
                  CSN2,
                  CSN4,
                  CSN6,
                  CSN8,
                  OPN2,
                  OPN4,
                  OPN6),
                ~abs(as.numeric(.x)-6),
                .names="{.col}_R")) %>%
  mutate(EXT =  select(.,
                       EXT1,
                       EXT2_R,
                       EXT3,
                       EXT4_R,
                       EXT5,
                       EXT6_R,
                       EXT7,
                       EXT8_R,
                       EXT9,
                       EXT10_R) %>%
           rowSums(na.rm = FALSE),
         EST =  select(.,
                       EST1_R,
                       EST2,
                       EST3_R,
                       EST4,
                       EST5_R,
                       EST6_R,
                       EST7_R,
                       EST8_R,
                       EST9_R,
                       EST10_R) %>%
           rowSums(na.rm = FALSE),
         AGR =  select(.,
                       AGR1_R,
                       AGR2,
                       AGR3_R,
                       AGR4,
                       AGR5_R,
                       AGR6,
                       AGR7_R,
                       AGR8,
                       AGR9,
                       AGR10) %>%
           rowSums(na.rm = FALSE),
         CSN =  select(.,
                       CSN1,
                       CSN2_R,
                       CSN3,
                       CSN4_R,
                       CSN5,
                       CSN6_R,
                       CSN7,
                       CSN8_R,
                       CSN9,
                       CSN10) %>%
           rowSums(na.rm = FALSE),
         OPN =  select(.,
                       OPN1,
                       OPN2_R,
                       OPN3,
                       OPN4_R,
                       OPN5,
                       OPN6_R,
                       OPN7,
                       OPN8,
                       OPN9,
                       OPN10) %>%
           rowSums(na.rm = FALSE)) %>%
  mutate(GFP = select(.,
                      EXT,
                      EST,
                      AGR,
                      CSN,
                      OPN) %>%
           rowSums(na.rm = FALSE))

# Save dataset
# write.csv(df, "GFP_revised.csv", row.names = FALSE)

# Sample
df1 <- df %>%
  subset(!is.na(GFP))

# Subsample
df2 <- df1 %>%
  subset(IPC==1)

# Scaled
df_scaled <- df2 %>%
  select(EXT1,
         EXT2_R,
         EXT3,
         EXT4_R,
         EXT5,
         EXT6_R,
         EXT7,
         EXT8_R,
         EXT9,
         EXT10_R:EST1_R,
         EST2,
         EST3_R,
         EST4,
         EST5_R:AGR1_R,
         AGR2,
         AGR3_R,
         AGR4,
         AGR5_R,
         AGR6,
         AGR7_R,
         AGR8:CSN1,
         CSN2_R,
         CSN3,
         CSN4_R,
         CSN5,
         CSN6_R,
         CSN7,
         CSN8_R,
         CSN9:OPN1,
         OPN2_R,
         OPN3,
         OPN4_R,
         OPN5,
         OPN6_R,
         OPN7:OPN10) %>%
  scale() %>%
  as.data.frame()

# Correlations
df_cor <- df_scaled %>%
  cor()

ggcorrplot::ggcorrplot(df_cor)

#Eigen value
df_eigen <- eigen(df_cor)
df_eigen$values

#Proportion of variance explained
df_var_prop <- df_eigen$values/sum(df_eigen$values)
df_var_prop

#Eigen vectors
df_eigen$vectors

#Scree
plot(df_eigen$values, xlab = "Principal Component", ylab = "Eigen Values", type = "b")

plot(cumsum(df_var_prop), xlab = "Principal Component", ylab = "Cumulative Eigen Values", type = "b")


#Calculate initial factor loadings
pc1 <- psych::principal(df_scaled,
                        nfactors = length(df_scaled),
                        rotate = "none")
pc1

fa1 <- psych::fa(df_scaled,
                 rotate = "none",
                 fm = "ml")

fa1$values

fa2 <- factanal(df_scaled,
                factors = 1)

fa2

psych::scree(df_scaled, pc=FALSE)
psych::scree(df_scaled, pc=TRUE)


efa1 <- df_scaled %>%
  jmv::efa(
    nFactorMethod = "fixed",
    extraction = "ml",
    rotation = "none",
    hideLoadings = 0,
    screePlot = TRUE,
    eigen = TRUE,
    factorSummary = TRUE,
    modelFit = TRUE,
    kmo = TRUE)

efa1














# Correlations
df2_cor <- df2 %>%
  select(EXT1,
         EXT2_R,
         EXT3,
         EXT4_R,
         EXT5,
         EXT6_R,
         EXT7,
         EXT8_R,
         EXT9,
         EXT10_R:EST1_R,
         EST2,
         EST3_R,
         EST4,
         EST5_R:AGR1_R,
         AGR2,
         AGR3_R,
         AGR4,
         AGR5_R,
         AGR6,
         AGR7_R,
         AGR8:CSN1,
         CSN2_R,
         CSN3,
         CSN4_R,
         CSN5,
         CSN6_R,
         CSN7,
         CSN8_R,
         CSN9:OPN1,
         OPN2_R,
         OPN3,
         OPN4_R,
         OPN5,
         OPN6_R,
         OPN7:OPN10) %>%
  cor()

ggcorrplot::ggcorrplot(df2_cor)




# Archived
# Correlation matrix
df2 %>%
  select(.,
         EXT1:EXT10) %>%
  jmv::corrMatrix()

df2 %>%
  select(.,
         EST1_R,
         EST2,
         EST3_R,
         EST4,
         EST5_R:EST10_R) %>%
  jmv::corrMatrix()
         

df2 %>%
  select(.,
         EXT1,
         EXT2_R,
         EXT3,
         EXT4_R,
         EXT5,
         EXT6_R,
         EXT7,
         EXT8_R,
         EXT9,
         EXT10_R:EST1_R,
         EST2,
         EST3_R,
         EST4,
         EST5_R:AGR1_R,
         AGR2,
         AGR3_R,
         AGR4,
         AGR5_R,
         AGR6,
         AGR7_R,
         AGR8:CSN1,
         CSN2_R,
         CSN3,
         CSN4_R,
         CSN5,
         CSN6_R,
         CSN7,
         CSN8_R,
         CSN9:OPN1,
         OPN2_R,
         OPN3,
         OPN4_R,
         OPN5,
         OPN6_R,
         OPN7:OPN10) %>%
  jmv::corrMatrix()

# Exploratory Factor Analysis
efa1 <- df2 %>%
  select(.,
         EXT1,
         EXT2_R,
         EXT3,
         EXT4_R,
         EXT5,
         EXT6_R,
         EXT7,
         EXT8_R,
         EXT9,
         EXT10_R:EST1_R,
         EST2,
         EST3_R,
         EST4,
         EST5_R:AGR1_R,
         AGR2,
         AGR3_R,
         AGR4,
         AGR5_R,
         AGR6,
         AGR7_R,
         AGR8:CSN1,
         CSN2_R,
         CSN3,
         CSN4_R,
         CSN5,
         CSN6_R,
         CSN7,
         CSN8_R,
         CSN9:OPN1,
         OPN2_R,
         OPN3,
         OPN4_R,
         OPN5,
         OPN6_R,
         OPN7:OPN10) %>%
  jmv::efa(
    nFactorMethod = "fixed",
    extraction = "ml",
    rotation = "none",
    hideLoadings = 0,
    screePlot = TRUE,
    eigen = TRUE,
    factorSummary = TRUE,
    modelFit = TRUE,
    kmo = TRUE)

efa1

# Factor Analysis
df2 %>%
  select(.,
         EXT1,
         EXT2_R,
         EXT3,
         EXT4_R,
         EXT5,
         EXT6_R,
         EXT7,
         EXT8_R,
         EXT9,
         EXT10_R:EST1_R,
         EST2,
         EST3_R,
         EST4,
         EST5_R:AGR1_R,
         AGR2,
         AGR3_R,
         AGR4,
         AGR5_R,
         AGR6,
         AGR7_R,
         AGR8:CSN1,
         CSN2_R,
         CSN3,
         CSN4_R,
         CSN5,
         CSN6_R,
         CSN7,
         CSN8_R,
         CSN9:OPN1,
         OPN2_R,
         OPN3,
         OPN4_R,
         OPN5,
         OPN6_R,
         OPN7:OPN10) %>%
  jmv::efa(
    nFactorMethod = "fixed",
    nFactors = 5,
    rotation = "oblimin",
    extraction = "ml",
    hideLoadings = 0.2,
    factorCor = TRUE,
    factorSummary = TRUE)

# Correlation matrix
df2 %>%
  select(.,
         EXT,
         EST,
         AGR,
         CSN,
         OPN) %>%
  jmv::corrMatrix()

# Exploratory Factor Analysis
df2 %>%
  select(.,
         EXT,
         EST,
         AGR,
         CSN,
         OPN) %>%
  jmv::efa(
    nFactorMethod = "fixed",
    rotation = "none",
    extraction = "ml",
    hideLoadings = 0,
    screePlot = TRUE,
    eigen = TRUE,
    kmo = TRUE)

# Factor Analysis
df2 %>%
  select(.,
         EXT,
         EST,
         AGR,
         CSN,
         OPN) %>%
  jmv::efa(
    nFactorMethod = "fixed",
    nFactors = 1,
    rotation = "oblimin",
    extraction = "ml",
    hideLoadings = 0.0,
    factorCor = TRUE,
    factorSummary = TRUE)

# Factor Analysis
df2 %>%
  select(.,
         EXT,
         EST,
         AGR,
         CSN,
         OPN) %>%
  jmv::efa(
    nFactorMethod = "fixed",
    nFactors = 2,
    rotation = "oblimin",
    extraction = "ml",
    hideLoadings = 0.2,
    factorCor = TRUE,
    factorSummary = TRUE)

# Factor Analysis
df2 %>%
  select(.,
         EXT1,
         EXT2_R,
         EXT3,
         EXT4_R,
         EXT5,
         EXT6_R,
         EXT7,
         EXT8_R,
         EXT9,
         EXT10_R:EST1_R,
         EST2,
         EST3_R,
         EST4,
         EST5_R:AGR1_R,
         AGR2,
         AGR3_R,
         AGR4,
         AGR5_R,
         AGR6,
         AGR7_R,
         AGR8:CSN1,
         CSN2_R,
         CSN3,
         CSN4_R,
         CSN5,
         CSN6_R,
         CSN7,
         CSN8_R,
         CSN9:OPN1,
         OPN2_R,
         OPN3,
         OPN4_R,
         OPN5,
         OPN6_R,
         OPN7:OPN10) %>%
  jmv::efa(
    nFactorMethod = "fixed",
    nFactors = 1,
    rotation = "oblimin",
    extraction = "ml",
    hideLoadings = 0.2,
    factorCor = TRUE,
    factorSummary = TRUE)

