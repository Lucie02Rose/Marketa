#######################
### GMO7 ASSIGNMENT ###
#######################

## RQ ###
# Do promoter variants cause reduced PCYT2 gene expression in children with Epilepsy?

### Notes to datasets ###
# diagnoses, inherited, causal non-coding variants, gene expression regulators
# pcyt2 - chr17, candidate motifs, epilepsy - multifactorial, noncoding region upstream
# sirtuin7, 500 and 500 patients - matched, unrelated
# 3 technical replicates for the gene expression

################################
### clean the R environment ####
################################

rm(list = ls(all = T))

#########################################
### get and set the working directory ###
#########################################

# getting and setting working directory - replace for the actual directory with data
getwd()
setwd("/home/Desktop/genomic_medicine/gmo7/assignment")
getwd()

#######################################
### installing and loading packages ###
#######################################

# genetics for HW equilibrium
if (!require(genetics)) {
  install.packages("genetics")}
library(genetics)
# devtools for downloading from github
if (!require(devtools)) {
  install.packages("devtools") }
library(devtools)
# Check and install 'BiocManager'
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")}
library(BiocManager)
# install 'snpStats', as LDHeatmap dependency
if (!require("snpStats", quietly = TRUE)) {
  BiocManager::install("snpStats", force = TRUE)}
# the 'LDheatmap' package from GitHub
if (!require("LDheatmap", quietly = TRUE)) {
  devtools::install_github("SFUStatgen/LDheatmap")}
library(LDheatmap)
# scientific notation
if (!require(scales)) {
  install.packages("scales")}
library(scales)
# tidyverse, dplyr and tidyr for handling tables
if (!require(tidyverse)) {
  install.packages("tidyverse")}
library(tidyverse)

if (!require(dplyr)) {
  install.packages("dplyr")}
library(dplyr)

if (!require(tidyr)) {
  install.packages("tidyr")}
library(tidyr)
# purrr for iterations
if (!require(purrr)) {
  install.packages("purrr")}
library(purrr)
# ggplot2 for plotting
if (!require(ggplot2)) {
  install.packages("ggplot2")}
library(ggplot2)
# for numbers and plots
if (!require(ggsignif)) {
  install.packages("ggsignif")}
library(ggsignif)
# cochran-armitage trend test
if (!require(DescTools)) {
  install.packages("DescTools")}
library(DescTools)
# dunn test as post hoc for KW test
if (!require(dunn.test)) {
  install.packages("dunn.test")}
library(dunn.test)
# ancova
if (!require(sandwich)) {
  install.packages("sandwich")}
library(sandwich)
# GMM mclust package
if (!require(mclust)) {
  install.packages("mclust")}
library(mclust)
# power calculations
if (!require(pwr)) {
  install.packages("pwr")}
library(pwr)
# randomForest for random forest algorithm
if (!require(randomForest)) {
  install.packages("randomForest")}
library(randomForest)
# caret for random forest metrics
if (!require(caret)) {
  install.packages("caret")}
library(caret)


########################
### loading datasets ###
########################

# ensure the directory is set up so that the code is in the same folder as datasets

cohort <- read_tsv("cohort.txt") # tab-delimited
gene_expression <- read_csv("gene_expression.csv") #comma-separated
variants <- read_tsv("variants.vcf") # tab-delimited

########################
### data exploration ###
########################

# inspecting the datasets for column names (head), missing values (na)
head(cohort)
head(gene_expression) 
head(variants) # all datasets share columns for patient IDs, data formats

sum(is.na(cohort))
sum(is.na(gene_expression))
sum(is.na(variants)) # no missing values for any of the datasets (all 0)

summary(gene_expression) # sirt7 max is 250, pcyt2 max is 310, all whole numbers
summary(variants) # locus position treated as a number, all QUAL above 200

table(cohort$phenotype) # there are 500 epilepsy and 500 unaffected patients
table(variants$POS, variants$ALT) # 4 different loci affected, ALT mostly A
table(variants$POS, variants$REF) # each has 1 ref and 1 alt allele type
table(variants$ID) # 2 known and 2 unknown SNVs
# rs373468901 (81911685), rs377442650 (81911686)
table(variants$ID, variants$POS) # later on to tell which variant has which locus
table(variants$REF) # 3 SNVs and one deletion
table(variants$POS, variants$GT) # 8 people homozygous recessive, all others 0/1

variants %>% # plotting variant IDs and positions - complementary to table
  ggplot(aes(x = ID, fill= as.factor(POS))) + 
  geom_bar() + 
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position="right")

variants %>% # plotting variant IDs and reference - complementary to table
  ggplot(aes(x = GT, fill= REF)) + 
  geom_bar() + 
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position="right")

variants %>% # plotting variant IDs and alternative - complementary to table
  ggplot(aes(x = GT, fill= ALT)) + 
  geom_bar() + 
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position="right")

# pivotting gene expression for better plotting 
gene_expression_long <- gene_expression %>%
  pivot_longer(cols = -patient_id, names_to = "Variable", values_to = "Counts")

# plotting count distributions (look quite normal within each "peak", multimodal)
ggplot(gene_expression_long, aes(x = Counts, color = Variable)) +
  geom_density() +
  theme_minimal() +
  labs(title = "Density Plot for gene expressions for all patients",
       x = "gene count",
       y = "density") # the pcyt2 has 3 peaks, sirt7 2 peaks ("levels of expression")

######################
### data wrangling ###
######################

# for gene_expression data:
# create mean columns since technical replicates
# using the rowMeans function, selecting appropriate columns 
# and rounding to whole numbers to create pcyt2_mean and sirt7_mean
# scaling them by min-max for counts to be continuous

gene_expression <- gene_expression %>%
  mutate(
    pcyt2_mean = round(rowMeans(select(., 2:4)), 0),
    sirt7_mean = round(rowMeans(select(., 5:7)), 0)) %>%
  mutate(
    pcyt2_mm = (pcyt2_mean-min(pcyt2_mean))/(max(pcyt2_mean) - min(pcyt2_mean)),
    sirt7_mm = (sirt7_mean-min(sirt7_mean))/(max(sirt7_mean) - min(sirt7_mean))
  )

# pivotting to the long version of gene_expression for...
mmlong <- gene_expression %>%
  pivot_longer(cols = c(pcyt2_mm, sirt7_mm), 
               names_to = "Variable", 
               values_to = "Counts")
# ...density plotting of the min-maxed variants
ggplot(mmlong, aes(x = Counts, color = Variable)) +
  geom_density() +
  theme_minimal() +
  labs(title = "Density Plot for PCYT2 and SIRT7",
       x = "normalised count (min-max)",
       y = "density",
       color = "gene")

# variants handling
variants <- variants %>%
  rename(
    patient_id = SAMPLE) %>% # rename sample column to patient_id to match other datasets
  mutate( # rename deletion so it is not that long
    REF =ifelse(REF =="ACCCACACCTGGCCTCTCCGCACCG", "WT", REF)) %>% # explain then what WT is for
  mutate(
    combined_genotype = case_when( # creating a combined genotype depending on alleles
      GT == "0/1" ~ paste(REF, ALT, sep = "/"), # heterozygous
      TRUE ~ paste(ALT, ALT, sep = "/"))) %>% # homozygous recessive
  mutate(
    POS = paste0("loc", POS)) # adding loc in front of all positions to not be treated as numbers

# creating values to store homozygous dom. alleles to then fill all wild-type people with
# homozygous dominant means choosing both REF alleles and group by locus
wt_alleles <- variants %>%
  group_by(POS) %>%                                
  summarize(
    wt = ifelse(POS == "loc81911677", "WT/WT", paste(unique(REF), unique(REF), sep = "/")),
    .groups = "drop"
  ) %>%
  deframe()

variants_pivoted <- variants %>%
  pivot_wider( # pivoting variants to create columns for all 4 loci 
    # (those are deterministic of variants)
    names_from = POS, # column to be divided is POS
    values_from = combined_genotype # new columns filled up with combined_genotype
  ) %>%
  group_by(patient_id) %>% # group by patient_ID fort the following operations 
  mutate(
    # collapse loc columns into new columns
    across(starts_with("loc"), ~ { 
      # since several patients have more than one mutation,
      # their rows need to be collapsed into one
      collapsed <- paste(na.omit(.), collapse = ", ")
      if (collapsed == "") NA else collapsed # rest is filled out with NA
    })
  ) %>%
  ungroup() %>% # ungroup
  distinct(patient_id, .keep_all = TRUE) %>% # get rid of now duplicate rows 
  # for people who have more than 1 mutation
  select(patient_id, starts_with("loc")) # select only columns needed


# combining the datasets - cohort, gene expression with variants
all_patients <- cohort %>%
  left_join( # combine cohort and vars by patient_id
    variants_pivoted, by = "patient_id")  %>% 
  left_join( # combine selected mean columns from gene expression 
             # with all previous by patient_id
    gene_expression %>% 
              select(patient_id, pcyt2_mm, sirt7_mm, pcyt2_mean, sirt7_mean), 
            by = "patient_id") %>%
  mutate(
    across( # add to all current NA values in loc columns the 
            # appropriate wild type allele by loc column in the wt_alleles
      all_of(names(wt_alleles)),
      ~ replace(., is.na(.), wt_alleles[cur_column()])
    )
  ) %>%
  mutate(
    across(where(is.character), as.factor) # mutate every character to factor
  )

head(all_patients) # viewing if it worked - it did

######################
### further graphs ###
######################

#plotting density graph for pcyt2 only to see how that looks
ggplot(all_patients, 
       aes(x = pcyt2_mm, fill = phenotype)) +
  geom_density(alpha = 0.5) + # alpha means how see-through
  theme_minimal() +
  labs(title = "Density Plot by Group",
       x = "Gene Count",
       y = "Density")

##############################
### genotype and phenotype ###
##############################

### quantifying and identifying variants present in each cohort ###
###################################################################

# How will you identify and quantify what variants are present in the Epilepsy cohort and the control group.
# How would you test for an association between the frequencies of any identified variant(s), and the phenotype?

# creating a list of contingency tables for each locus
loc_tables <- list(
  loc77 = table(all_patients$phenotype, all_patients$loc81911677),
  loc85 = table(all_patients$phenotype, all_patients$loc81911685),
  loc86 = table(all_patients$phenotype, all_patients$loc81911686),
  loc97 = table(all_patients$phenotype, all_patients$loc81912397))

# print each table out to see the genotypes for each locus by phenotype
loc_tables$loc77
loc_tables$loc85 #rs373468901
loc_tables$loc86 #rs377442650
loc_tables$loc97

# fisher test and cochran armitage test on count data using a for loop
# creating a list of the tables
test_results <- list() # list to store results

for (loc in names(loc_tables)) {
  
  loc_table <- loc_tables[[loc]] # for each locus in the list of locus tables
  
  fisher_result <- fisher.test(loc_table) # perform fisher
  
  cochran_armitage_result <- CochranArmitageTest(loc_table) # perform CA test
  
  test_results[[loc]] <- list(
    fisher = fisher_result,
    cochran_armitage = cochran_armitage_result # store test results for both in a list
  )
}

test_results  # print the results of the tests (fo all loci)

# odds ratios - not present for loc85 and loc86 (some zero counts in table)
# meaning that the odds of having epilepsy with homozygous recessive are infinite
# there is no significant association of phenotype and genotype in loc81912397, 
# nor in loc81911677. There is an association in the two remaining loci, 
# although the problem is the zero counts in homozygous recessive individuals
# power of the test needed

### plotting count plots for each locus ###
###########################################

# doing the same as above but directly storing the p values in the plots
# note that frequency in this dataset would just be everything scaled down by 1000
# it is better to run the nonparametric tests on whole numbers not frequencies (continuous)

# extracting column names for loci
loc_cols <- grep("^loc", names(all_patients), value = TRUE)

# freq_tables - list for storing raw contingency tables
count_tables <- list()

for (loc in loc_cols) {
  # a raw contingency table for each locus
  loc_table <- table(all_patients$phenotype, all_patients[[loc]])
  
  # storing into count list
  count_tables[[loc]] <- loc_table
  
  # fisher test for all loci and storing the p value
  fisher_result <- fisher.test(loc_table)
  fisher_pvalue <- fisher_result$p.value  
  
  # CA trend test for all loci and storing the p values
  cochran_armitage_result <- CochranArmitageTest(loc_table)
  cochran_armitage_pvalue <- cochran_armitage_result$p.value  
  
  # converting to more meaningful notation with 2 sig figs
  fisher_pvalue_formatted <- format(fisher_pvalue, scientific = TRUE, digits = 2)
  cochran_armitage_pvalue_formatted <- format(cochran_armitage_pvalue, scientific = TRUE, digits = 2)
  
  # frequency tables for plotting - renaming some variables which were created by default
  count_df_tab <- as.data.frame(loc_table) %>%
    rename(phenotype = Var1, genotype = Var2, count = Freq)
  
  # adding significance stars for fisher's and CA p values
  fisher_star <- ifelse(fisher_pvalue < 0.05, "*", "")
  cochran_armitage_star <- ifelse(cochran_armitage_pvalue < 0.05, "*", "")
  
  # plotting a bar plot for each locus by phenotype
  plot <- ggplot(count_df_tab, 
                 aes(x = genotype, y = count, fill = phenotype)) + 
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_brewer(palette = "Set1") + # palette set1 used throughout for consistency
    theme(legend.position = "right") + # legend here is distracting
    labs(x = "allelic genotype", y = "count", title = paste("Locus", loc))
  
  # adding p values to the bar plots and position them meaningfully (x and y axis)
  plot <- plot +
    geom_text(aes(x = 1, y = max(count) * 0.9, label = paste("Fisher p =", fisher_pvalue_formatted, fisher_star)),
              size = 5, hjust = 0.5, color = "black") +
    geom_text(aes(x = 1, y = max(count) * 0.8, label = paste("Cochran p =", cochran_armitage_pvalue_formatted, cochran_armitage_star)),
              size = 5, hjust = 0.5, color = "black")
  
  # printing the plots
  print(plot)
}


### power of the Fisher tests ###
#################################

# power of the test - by simulation, number of simulations at 1000 and significance of 0.05
set.seed(123) # seed for consistency

fisher_power_simulation <- function(loc_table, n_simulations = 1000, alpha = 0.05) {
  significant_count <- 0  # counter for significant results
  
  # simulating n_simulations datasets with the same marginal totals as the observed locus tables
  for (i in 1:n_simulations) {
    # simulating a table using the r r2dtable
    simulated_table <- r2dtable(1, rowSums(loc_table), colSums(loc_table))[[1]]
    # applying fisher test to the simulated table
    fisher_result <- fisher.test(simulated_table)
    # is the p value less than significance threshold alpha - a significant count
    if (fisher_result$p.value < alpha) {
      significant_count <- significant_count + 1
    }
  }
  # power is the number of significant results in the simulations divided by no of simulations
  power <- significant_count / n_simulations
  return(power)
}

power_results <- list() # list for power results

# iterating through contingency tables
# doing the power calculation by applying the simulation function from above
for (loc in names(count_tables)) {
  loc_table <- count_tables[[loc]] 
  power_results[[loc]] <- fisher_power_simulation(loc_table)
}

print(power_results) # printing the results
# the results suggest that there is not sufficient power of the Fisher test
# locus 86 very around the threshold (marginally)
# Fisher test is likely not sufficiently powered because of the low counts for 
# the homozygous recessive (or heterozygous in loc 77) individuals
# study not sufficiently powered by the counts


# Locus 81911685 - rs373468901, Locus 81911686 - rs377442650
# ACCCACACCTGGCCTCTCCGCACCG = WT

#########################################################
### Hardy-Weinberg equilibrium for the loci genotypes ###
#########################################################

# Testing whether the identified variants are in Hardy-Weinberg Equilibrium (HWE) 
# in the sampled population.
# Perform Hardy-Weinberg equilibrium test on the genotypes - all genotypes are in the HWE

hwe_test_fcn <- function(genotypes, label) { # making a function to do all at once
  cat("Results for", label, "\n")
  print(HWE.test(genotypes)) # applying the hwe test on the genotypes
  print(HWE.chisq(genotypes)) # using the chi square test from the package
  cat("\n")
}

# HWE is independent of phenotypes so it should be for all patients together
# using each locus column from the all_patients dataset
for (locus in loc_cols) {
  genotypes <- genotype(all_patients[[locus]])
  hwe_test_fcn(genotypes, paste(locus, "in all patients"))
}

# some counts less than 5 so df does not work (Chisq test worse than Fisher), 
# but the p values all indicate equilibrium, another method would then be running a simulation

##############################################################
### Linkage disequilibrium and multiple testing adjustment ###
##############################################################

# a problem for a future regression model - LD high and problem with collinearity, same haplotype
# null hypothesis of independence
# adjustment for multiple testing (bonferroni)

# a lots of multiple testing performed in this study, always need for adjustment 
# of the p values - bonferroni is stricter than Benjamini Hochberg used elsewhere

# running genotype on all patients all loci - again LD is not about the phenotype
all77 <- genotype(all_patients$loc81911677)
all85 <- genotype(all_patients$loc81911685)
all86 <- genotype(all_patients$loc81911686)
all97 <- genotype(all_patients$loc81912397)

# putting all of them into a dataframe
LD_all_vars <- data.frame(all77, all85, all86, all97)
# extracting the D' values
LD.Mat.D <- LD(LD_all_vars)$"D'" 
LD.Mat.D
# extracting the p values
LD.Mat.p <- LD(LD_all_vars)$"P-value"
LD.Mat.p # print p values
# making a matrix from the p values
colnames(LD.Mat.p) <- rownames(LD.Mat.p) <- c("all77", "all85", "all86", "all97")
# extracting the upper part of the matrix (it is symmetrical)
p_values <- LD.Mat.p[upper.tri(LD.Mat.p, diag = FALSE)]
# removing any NA values
p_values <- p_values[!is.na(p_values)]

# adjusting for multiple testing using bonferroni - not too many comparisons here (6)
adjusted_p <- p.adjust(p_values, method = "bonferroni") 
# adjusted p values matrix (upper corner again)
LD.Mat.p.adj <- LD.Mat.p
LD.Mat.p.adj[upper.tri(LD.Mat.p.adj)][!is.na(LD.Mat.p[upper.tri(LD.Mat.p)])] <- adjusted_p
LD.Mat.p.adj # print the adjusted p values

# necessary steps for figure
# positions of the different snps
snp_positions <- c(81911677, 81911685, 81911686, 81912397)
# actual plot
LD_obj <- LDheatmap(LD.Mat.D, LDmeasure = "D'",
          genetic.distances = snp_positions, 
          SNP.name = c("all77", "all85", "all86", "all97"),
          color = colorRampPalette(c("red", "blue", "white"))(20),
          title = "Pairwise LD (D')",
          add.map = TRUE)
# all in LD apart from 77-97
# permutation test for 85-97 and 86-97 - the D' values are quite low
#####################################################################
set.seed(123) # setting a seed

# selecting the columns
loc85_86 <- list(loc85 = all_patients$loc81911685, loc86 = all_patients$loc81911686)
loc97 <- all_patients$loc81912397

results <- list() # list for results
n_perm <- 1000 # no of permutations

# 'computing the LD' function, using the genotype and returning D prime
compute_ld <- function(geno1, geno2) {
  geno1 <- genotype(geno1)
  geno2 <- genotype(geno2)
  ld_stats <- LD(geno1, geno2)
  return(ld_stats$"D'")         
}

# one of the loci is kept the same and for the other one the genotypes are
# shuffled around, in this case for loc97 in both cases they are shuffled
for (locus_name in names(loc85_86)) {
  locus <- loc85_86[[locus_name]]
  observed_ld <- compute_ld(locus, loc97) # applying the function for LD
  permuted_ld <- numeric(n_perm) #permuted - convert to numeric
  
  # this iterates through both the list of loci and the permutations (1000)
  for (i in 1:n_perm) {
    permuted_loc97 <- sample(loc97) 
    permuted_ld[i] <- compute_ld(locus, permuted_loc97)
  }
    empirical_p <- mean(permuted_ld >= observed_ld) # there is an empirical p value
  results[[locus_name]] <- list(observed_ld = observed_ld, p_value = empirical_p)
  # printing out the output by concatenating
  cat(locus_name, "vs loc97: observed D':", observed_ld, ", p-val from permutation=", empirical_p, "\n")
  
}

# not significantly different from what could be expected under the null hypothesis
# of no linkage disequilibrium - the p values are too high for both now (about 0.6)
# this permutation test is a post hoc, it can be said that there is some degree of LD
# just by the package but doing a permutation - not significantly different

########################################################
### Gene expression difference depending on genotype ###
########################################################

# making separate list of tibbles for the loci
loci_subsets <- list()

# in the tibbles per locus, including the min maxed genes and phenotype as well
for (locus in loc_cols) {
  loc_subset <- all_patients %>%
    select(phenotype, pcyt2_mm, sirt7_mm, !!sym(locus)) %>% # selecting columns
    rename(genotype = !!sym(locus)) %>% # renaming
    mutate(locus = locus)
  loci_subsets[[locus]] <- loc_subset
  assign(paste0(locus, "_data"), loc_subset) # naming each tibble
} 

plots <- list() # list for plots
# creating plots for all loci for genotypes separated by phenotype
for (locus in names(loci_subsets)) {
  # plots for pcyt2 - independent (genotype) on x axis, scatterplot - jittered to be clearer
  plot_pcyt2 <- ggplot(loci_subsets[[locus]], aes(x = genotype, y = pcyt2_mm, color = phenotype)) +
    geom_jitter(width = 0.3, height = 0, size = 1.2, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = paste("PCYT2 expression per genotype for", locus), # titles for all pcyt2 plots
      x = "genotype",
      y = "PCYT2 expression" # titles for axes
    )
  # plots for sirt7 - independent (genotype) on x axis, scatterplot - jittered to be clearer
  plot_sirt7 <- ggplot(loci_subsets[[locus]], aes(x = genotype, y = sirt7_mm, color = phenotype)) +
    geom_jitter(width = 0.3, height = 0, size = 1.2, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = paste("SIRT7 expression per genotype for", locus), # titles for all sirt7 plots
      x = "genotype",
      y = "SIRT7 expression" # titles for axes
    )
  # save both types of plots to the list
  plots[[paste0(locus, "_PCYT2")]] <- plot_pcyt2
  plots[[paste0(locus, "_SIRT7")]] <- plot_sirt7
}
# display plots for each locus - all quite informative
for (locus in names(plots)) {
  print(plots[[locus]])
}

# in all graphs the categorical independent variable (genotype) on the x axis
# continuous variable - dependent (gene expression) is on y axis throughout

# testing for normalty using the shapiro test
# again in a for loop for each locus and for both genes within the locus and 
# extracting the p values

for (locus_name in names(loci_subsets)) {
  
  locus <- loci_subsets[[locus_name]]

  pcyt2_shapiro <- shapiro.test(locus$pcyt2_mm)
  sirt7_shapiro <- shapiro.test(locus$sirt7_mm)
  
  print(paste(locus_name, "Pcyt2 p-value:", pcyt2_shapiro$p.value))
  print(paste(locus_name, "Sirt7 p-value:", sirt7_shapiro$p.value))
}

# as per the graphs and the Shapiro Wilk normalty test (p values below threshold)
# Shapiro null is that there is normalty (we can reject null due to p vals)
# none of the gene expressions exhibit normal distribution, 
# rather the expression seems to be linked to dosage of alleles
# this is no surprise as gene expression is multimodal

##############################
### graphing density plots ###
##############################

# only plotting the density plots for pcyt2-genotype relating to the research question 
# as sirt7 results not stat significant apart from locus 77'
# for each locus, again with a for loop by locus storing all in a loop

pcyt2_plots <- list()

for (locus in names(loci_subsets)) {
  subset_data <- loci_subsets[[locus]]  
  
  pcytplots <- ggplot(subset_data, aes(x = pcyt2_mm, fill = genotype)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    labs(title = paste("Density Plot for PCYT2 at", locus),
         x = "normalised count (min-max)",
         y = "density")
  
  pcyt2_plots[[paste0(locus, "_PCYT2")]] <- pcytplots
}

for (locus in names(pcyt2_plots)) {
  print(pcyt2_plots[[locus]])
} # then printing all in the list

############################
### non-parametric tests ###
############################

# Since this data is not normally distributed and there are discrete variables, 
# a non-parametric test should be applied, since there are more groups
# kruskal wallis test and a permutation simulation
# note that Kruskal Wallis tests significance of difference in group mean rank sums, is ranked

# a vector to store the p values of the KW test
kruskal_p_values <- c()

# running the KW test for each phenotype separately (all, epilepsy, unaffected)
# testing for how genotypes affect the gene expression
# iterating through all 4 loci
for (locus in names(loci_subsets)) {
  
  loc_all <- loci_subsets[[locus]]
  # pcyt2 and sirt7 for all patients - p values concatenated
  pcyt2_pvalue <- kruskal.test(pcyt2_mm ~ genotype, data = loc_all)$p.value 
  cat(locus, "| all | pcyt2 | p val:", pcyt2_pvalue, "\n")
  sirt7_pvalue <- kruskal.test(sirt7_mm ~ genotype, data = loc_all)$p.value
  cat(locus, "| all | sirt7 | p val:", sirt7_pvalue, "\n")
  kruskal_p_values <- c(kruskal_p_values, pcyt2_pvalue, sirt7_pvalue)
  
  loc_epilepsy <- loci_subsets[[locus]] %>% filter(phenotype == "Epilepsy")
  # pcyt2 and sirt7 for all epileptics - p values concatenated
  pcyt2_pvalue <- kruskal.test(pcyt2_mm ~ genotype, data = loc_epilepsy)$p.value
  cat(locus, "| epilepsy | pcyt2 | p val:", pcyt2_pvalue, "\n")
  sirt7_pvalue <- kruskal.test(sirt7_mm ~ genotype, data = loc_epilepsy)$p.value
  cat(locus, "| epilepsy | sirt7 | p val:", sirt7_pvalue, "\n")
  kruskal_p_values <- c(kruskal_p_values, pcyt2_pvalue, sirt7_pvalue)
  
  loc_unaffected <- loci_subsets[[locus]] %>% filter(phenotype == "unaffected")
  # pcyt2 and sirt7 for all unaffected - p values concatenated
  pcyt2_pvalue <- kruskal.test(pcyt2_mm ~ genotype, data = loc_unaffected)$p.value
  cat(locus, "| unaffected | pcyt2 | p val:", pcyt2_pvalue, "\n")
  sirt7_pvalue <- kruskal.test(sirt7_mm ~ genotype, data = loc_unaffected)$p.value
  cat(locus, "| unaffected | sirt7 | p val:", sirt7_pvalue, "\n")
  kruskal_p_values <- c(kruskal_p_values, pcyt2_pvalue, sirt7_pvalue)
}

# adjusting for multiple testing - both bonferroni (stricter) and BH
# with 0.05 alpha, the end conclusions do not differ between adjustment tests
adjusted_kruskal_p_values <- p.adjust(kruskal_p_values, method = "bonferroni")
bhadjusted_kruskal_p_values <- p.adjust(kruskal_p_values, method = "BH")
# printing p values and adjusted p values
kruskal_p_values
adjusted_kruskal_p_values
bhadjusted_kruskal_p_values

# performing the dunn test if the KW test is significant, again stratifying by phenotype
# all patients - for loop for loci subsets
# always printing p adjusted for the Dunn test

for (locus in names(loci_subsets)) {
  loc_all <- loci_subsets[[locus]]
  loc_all$genotype <- factor(loc_all$genotype)
  # pcyt2 gene,  if statement to do the dunn test if significant (below 0.05)
  pcyt2_kw_all <- kruskal.test(pcyt2_mm ~ genotype, data = loc_all)
  cat(locus, "| pcyt2 | KW p val:", pcyt2_kw_all$p.value, "\n")
  if (pcyt2_kw_all$p.value < 0.05) {
    dunn_pcyt2_all <- dunn.test(x = loc_all$pcyt2_mm, g = loc_all$genotype, kw = TRUE)
    print(dunn_pcyt2_all$P.adjusted)
  }
  # sirt7 gene, if statement to do the dunn test if significant (below 0.05)
  sirt7_kw_all <- kruskal.test(sirt7_mm ~ genotype, data = loc_all)
  cat(locus, "| sirt7 | KW p val:", sirt7_kw_all$p.value, "\n")
  if (sirt7_kw_all$p.value < 0.05) {
    dunn_sirt7_all <- dunn.test(x = loc_all$sirt7_mm, g = loc_all$genotype, kw = TRUE)
    print(dunn_sirt7_all$P.adjusted)
  }
}
# epileptics - for loop for loci subsets, subset by epileptics and convert to factors
for (locus in names(loci_subsets)) {
  loc_epilepsy <- loci_subsets[[locus]] %>% filter(phenotype == "Epilepsy")
  loc_epilepsy$genotype <- factor(loc_epilepsy$genotype)
  # pcyt2 gene, if statement to do the dunn test if significant (below 0.05)
  pcyt2_kw_epi <- kruskal.test(pcyt2_mm ~ genotype, data = loc_epilepsy)
  cat(locus, "| pcyt2 | KW p val:", pcyt2_kw_epi$p.value, "\n")
  if (pcyt2_kw_epi$p.value < 0.05) {
    dunn_pcyt2_epi <- dunn.test(x = loc_epilepsy$pcyt2_mm, g = loc_epilepsy$genotype, kw = TRUE)
    print(dunn_pcyt2_epi$P.adjusted)
  }
  # sirt7 gene, if statement to do the dunn test if significant (below 0.05)
  sirt7_kw_epi <- kruskal.test(sirt7_mm ~ genotype, data = loc_epilepsy)
  cat(locus, "| sirt7 | KW p val:", sirt7_kw_epi$p.value, "\n")
  if (sirt7_kw_epi$p.value < 0.05) {
    dunn_sirt7_epi <- dunn.test(x = loc_epilepsy$sirt7_mm, g = loc_epilepsy$genotype, kw = TRUE)
    print(dunn_sirt7_epi$P.adjusted)
  }
}
# unaffected - subset again for unaffected and convert to factor
for (locus in names(loci_subsets)) {
  loc_unaffected <- loci_subsets[[locus]] %>% filter(phenotype == "unaffected")
  loc_unaffected$genotype <- factor(loc_unaffected$genotype)
  # pcyt2 gene, if statement to do the dunn test if significant (below 0.05)
  pcyt2_kw_unaf <- kruskal.test(pcyt2_mm ~ genotype, data = loc_unaffected)
  cat(locus, "| pcyt2 | KW p val:", pcyt2_kw_unaf$p.value, "\n")
  if (pcyt2_kw_unaf$p.value < 0.05) {
    dunn_pcyt2_unaf <- dunn.test(x = loc_unaffected$pcyt2_mm, g = loc_unaffected$genotype, kw = TRUE)
    print(dunn_pcyt2_unaf$P.adjusted)
  }
  # sirt7 gene, if statement to do the dunn test if significant (below 0.05)
  sirt7_kw_unaf <- kruskal.test(sirt7_mm ~ genotype, data = loc_unaffected)
  cat(locus, "| sirt7 | KW p val:", sirt7_kw_unaf$p.value, "\n")
  if (sirt7_kw_unaf$p.value < 0.05) {
    dunn_sirt7_unaf <- dunn.test(x = loc_unaffected$sirt7_mm, g = loc_unaffected$genotype, kw = TRUE)
    print(dunn_sirt7_unaf$P.adjusted)
  }
}
# Dunn test results are the same regardless of stratification, the significant
# results are PCYT2 for loc 77', 85', 86' and SIRT7 for loc 77'
# for loci 85' and 86', genotype effects all important 
# (recessive have significantly lower expression)

########################
### permutation test ###
########################

# permutation can be run with the t statistic even though the t test itself is
# not applicable for the data, similarly to the lesson example
# the p value itself is computed on the permutations

# doing a similar test by permutation - genotype and gene expression
permutation_test <- function(genotype, gene_expression, n_permutations = 1000) {
  # a vector to store p values for each genotype level
  p_values <- numeric(length = length(unique(genotype)))
  names(p_values) <- unique(genotype)
  set.seed(123)
  # perform permutations similarly according to the lesson examples - t test
  for (factor_level in unique(genotype)) {
    # calculating the t statistics for all genotypes
    t_obs <- t.test(
      x = gene_expression[genotype == factor_level], 
      y = gene_expression[genotype != factor_level]
    )$statistic
    
    # vector for storing permuted t-statistics
    t_perm <- rep(NA, n_permutations)
    # doing 1000 permutations
    for (i in 1:n_permutations) {
      # and randomly shuffling the genotypes
      random_genotype <- sample(genotype)
      
      # calculating the permuted t-statistic for shuffled genotypes
      t_perm[i] <- t.test(
        x = gene_expression[random_genotype == factor_level], 
        y = gene_expression[random_genotype != factor_level]
      )$statistic
    }
    # permutation p-value (proportion of permuted t-statistics greater than or equal to observed t-statistic)
    p_values[factor_level] <- (sum(abs(t_perm) >= abs(t_obs)) + 1) / (n_permutations + 1)
  }
  return(p_values) # return the p values for each factor level
}

# using the above function by iterating over the four loci
results <- list() # list to store the results
 
# looping over each and subsetting
for (locus in names(loci_subsets)) {
  locus_data <- loci_subsets[[locus]]
  # extracting the genotypes and gene expressions for each
  genotypes <- locus_data$genotype
  gene_expression_pcyt2 <- locus_data$pcyt2_mm
  gene_expression_sirt7 <- locus_data$sirt7_mm
  
  # permutation test for the pcyt2
  p_values_pcyt2 <- permutation_test(genotypes, gene_expression_pcyt2, n_permutations = 1000)
  # permutation test for the sirt7
  p_values_sirt7 <- permutation_test(genotypes, gene_expression_sirt7, n_permutations = 1000)
  # storing the results in the list
  results[[locus]] <- list(
    p_values_pcyt2 = p_values_pcyt2,
    p_values_sirt7 = p_values_sirt7
  )
}
print(results) # print the results for each locus

# applying a correction test for multiple testing (BH)
all_p_values <- unlist(c(results$loc81911685$p_values_pcyt2, results$loc81911685$p_values_sirt7, 
                         results$loc81911686$p_values_pcyt2, results$loc81911686$p_values_sirt7,
                         results$loc81911677$p_values_pcyt2, results$loc81911677$p_values_sirt7,
                         results$loc81912397$p_values_pcyt2, results$loc81912397$p_values_sirt7))

adjusted_p_values <- p.adjust(all_p_values, method = "BH")
print(adjusted_p_values) # printing the adjusted values

# here, loc 85', 86' genotypes important for PCYT2
# loc 77' important for SIRT7

###################################
### gene expression correlation ###
###################################

genes <- cor.test(all_patients$pcyt2_mm, all_patients$sirt7_mm, method = "spearman")
genes # the two gene expressions are not correlated at all on each other

# in between genotypes it does not make much sense to discretize genotype into 
# numbers and use spearman (a bit redundant)

# problem with the multimodalty is that not all the modes have a representation of
# all genotypes and phenotypes - so it also doesn't make sense to divide them

#####################################
### gene expression and phenotype ###
#####################################

# here it would also be possible to run a simulation (in theory)

# non-parametric test for pcyt2 expression and phenotype
kruskal_pcyt2_patients <- kruskal.test(pcyt2_mm ~ phenotype, data = all_patients)
kruskal_pcyt2_patients

dunn_test_pcyt2 <- dunn.test(x = all_patients$pcyt2_mm, g = all_patients$phenotype, kw = TRUE)
dunn_test_pcyt2 # post-hoc test to tell the difference in medians and which group is lower

# plotting a violin plot of the phenotype and pcyt2 min-maxed expression
kw_pcyt2_p_value <- 0.01 
dunn_pcyt2_p_value <- 0.0052 # use the values from tests
dunn_pcyt2_z <- -2.56
kw_pcyt2_label <- paste0("KW p = ", format(kw_pcyt2_p_value, digits = 3), "*")
dunn_pcyt2_label <- paste0("Dunn p = ", format(dunn_pcyt2_p_value, digits = 3), "*")
dunn_pcyt2_z_label <- paste0("Dunn Z = ", format(dunn_pcyt2_z, digits = 3))
# format the values from tests
# plot the violin plot - phenotype on x axis and pcyt2mm on y axis
# overlay box and violin plots
ggplot(all_patients, aes(x = phenotype, y = pcyt2_mm, fill = phenotype)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.color = "red", outlier.shape = 16, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue") +
  labs(
    title = "Violin Plot of PCYT2 expression by phenotype",
    x = "phenotype",
    y = "pcyt2 expression (min-maxed)"
  ) +
  # Add the KW p value result to the graph
  geom_text(
    label = kw_pcyt2_label,
    x = 1.5, 
    y = 0.3, 
    size = 5, hjust = 0.5, color = "black"
  ) +
  # Add Dunn p value result to the graph
  geom_text(
    label = dunn_pcyt2_z_label,
    x = 1.5, 
    y = 0.2, 
    size = 5, hjust = 0.5, color = "black"
  ) +
  # add Dunn z result to the graph
  geom_text(
    label = dunn_pcyt2_label,
    x = 1.5, 
    y = 0.1, 
    size = 5, hjust = 0.5, color = "black"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# non-parametric test for sirt7 expression and phenotype
kruskal_sirt7_patients <- kruskal.test(sirt7_mm ~ phenotype, data = all_patients)
kruskal_sirt7_patients

dunn_test_sirt7 <- dunn.test(x = all_patients$sirt7_mm, g = all_patients$phenotype, kw = TRUE)
dunn_test_sirt7 # post-hoc test for sirt7 expression - difference in medians

# save and format the values from the tests
kw_sirt7_p_value <- 0.44
dunn_sirt7_p_value <- 0.22
dunn_sirt7_z <- -0.78
kw_sirt7_label <- paste0("KW p = ", format(kw_sirt7_p_value, digits = 3))
dunn_sirt7_label <- paste0("Dunn p = ", format(dunn_sirt7_p_value, digits = 3))
dunn_sirt7_z_label <- paste0("Dunn Z = ", format(dunn_sirt7_z, digits = 3))

# plot the graph (same as the previous one)
ggplot(all_patients, aes(x = phenotype, y = sirt7_mm, fill = phenotype)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.color = "red", outlier.shape = 16, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue") +
  labs(
    title = "Violin Plot of SIRT7 expression by phenotype",
    x = "phenotype",
    y = "sirt7 expression (min-maxed)"
  ) +
  # Add KW result annotation at a fixed x-y location
  geom_text(
    label = kw_sirt7_label,
    x = 1.5, 
    y = 0.3, 
    size = 5, hjust = 0.5, color = "black"
  ) +
  # Add Dunn p value result to the graph
  geom_text(
    label = dunn_pcyt2_z_label,
    x = 1.5, 
    y = 0.2, 
    size = 5, hjust = 0.5, color = "black"
  ) +
  # Add Dunn test result annotation
  geom_text(
    label = dunn_sirt7_label,
    x = 1.5, 
    y = 0.1, 
    size = 5, hjust = 0.5, color = "black"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# logistic regression instead for the loci 
# pcyt2 logistic regression
glm_pcyt2 <- glm(phenotype ~ pcyt2_mean, all_patients, family = "binomial")
levels(all_patients$phenotype) # find out how the encoding is
summary(glm_pcyt2)
exp(0.012783)
# plot(glm_pcyt2)
# the log odds of being unaffected increase with increasing pcyt2 expression
# the unaffected is 1, epilepsy is 0 factor level (reference)
# therefore when pcyt2 expression is 0, the log odds of epilepsy is -3.46
# then the pcyt2_mean for unaffected increases by log odds of 0.012783 with each increase of pcyt2 unit
# meaning a higher likelihood of unaffected with higher pcyt2 expression
# problem with this model is that gene expression is very much not normal, nonlinear and multimodal
# this is not the most appropriate test for this dataset - null deviance is only about 36 
# (still better than null model though)

glm_sirt7 <- glm(phenotype ~ sirt7_mean, all_patients, family = "binomial")
summary(glm_sirt7)
# no statistical difference here, again very cautious with the interpretation
# glm cannot capture the collinearity, multimodality and ab-normal distributions

# stratification by locus - the same results (would jeust need to be put in a for loop)
# per locus KW and Dunn test for phenotype and gene expression 
for (locus in names(loci_subsets)) {
  loc_all_pheno <- loci_subsets[[locus]]
  loc_all_pheno$phenotype <- factor(loc_all_pheno$phenotype)
  # pcyt2 gene,  if statement to do the dunn test if significant (below 0.05)
  pcyt2_k_all <- kruskal.test(pcyt2_mm ~ phenotype, data = loc_all_pheno)
  cat(locus, "| pcyt2 | KW p val:", pcyt2_k_all$p.value, "\n")
  if (pcyt2_k_all$p.value < 0.05) {
    d_pcyt2_all <- dunn.test(x = loc_all_pheno$pcyt2_mm, g = loc_all_pheno$phenotype, kw = TRUE)
    print(d_pcyt2_all$P.adjusted)
  }
  # sirt7 gene, if statement to do the dunn test if significant (below 0.05)
  sirt7_k_all <- kruskal.test(sirt7_mm ~ phenotype, data = loc_all_pheno)
  cat(locus, "| sirt7 | KW p val:", sirt7_k_all$p.value, "\n")
  if (sirt7_k_all$p.value < 0.05) {
    d_sirt7_all <- dunn.test(x = loc_all_pheno$sirt7_mm, g = loc_all_pheno$phenotype, kw = TRUE)
    print(d_sirt7_all$P.adjusted)
  }
}

################################################
### factors all together - ANCOVA lm and glm ###
################################################

# answering the question 
# Do promoter variants cause reduced PCYT2 gene expression in children with Epilepsy?
#####################################################################################

# since the only significant result for sirt7 is the genotype effect at locus 77'
# analysing the relationship of pcyt2 and genotype and phenotype (predictors)
# storing results in a list
robust_ancova_results <- list()

for (locus in names(loci_subsets)) {
  # for loop for iterating through the loci and storing each
  locus_data <- loci_subsets[[locus]]
  # converting all categorical to factors (not necessary because already should be)
  locus_data$phenotype <- factor(locus_data$phenotype, levels = c("unaffected", "Epilepsy"))
  locus_data$genotype <- factor(locus_data$genotype)
  # performing an equivalent to ancova - fitting a full model including interaction
  # of genotype and phenotype (a linear model here, which is not the best for collinearity, 
  # multimodality and non-normal data - caution with the interpretation)
  model_full <- lm(pcyt2_mm ~ phenotype + genotype + phenotype:genotype, data = locus_data)
  # calculation of thw R squared for the model - proportion of variance in pcyt2 explained
  # by genotype and phenotype per locus
  R2_full <- summary(model_full)$r.squared
  # here fitting a reduced model without the interaction term
  model_reduced <- lm(pcyt2_mm ~ phenotype + genotype, data = locus_data)
  # and calculating the r squared for the reduced model - a high r squared for full model means
  # that interaction term is important
  R2_reduced <- summary(model_reduced)$r.squared
  # calculation of the Cohen's f2 which is measuring the effect size
  f2 <- (R2_full - R2_reduced) / (1 - R2_full)
  # here the robust standard errors (heteroscedasticity consistent errors) are calculated for the full model
  # makes it a bit more reliable - if variance of pcyt2mm expression differs across the groups
  robust_se <- sqrt(diag(vcovHC(model_full, type = "HC3")))
  # get the coefficients for the full model
  coefficients <- summary(model_full)$coefficients
  # here, pulling out the t statistics and calculating p values based on them
  robust_t_stats <- coefficients[, "Estimate"] / robust_se
  robust_p_values <- 2 * (1 - pnorm(abs(robust_t_stats)))
  # here p values are adjusted using the Benjamini Hochberg test
  adjusted_p_values <- p.adjust(robust_p_values, method = "BH")
  # the results for each locust are saved
  robust_ancova_results[[locus]] <- list(
    adjusted_p_values = adjusted_p_values,
    R_squared = R2_full,
    Cohen_f2 = f2
  )
  # results for each locus printed
  cat("'Robust ANCOVA' for", locus, "\n")
  print(data.frame(Coefficients = coefficients[, "Estimate"], 
                   Robust_SE = robust_se, 
                   Robust_T = robust_t_stats, 
                   Robust_P = robust_p_values, 
                   Adjusted_P = adjusted_p_values, 
                   R_Squared = R2_full,
                   Cohen_f2 = f2))
} 

# it would not make sense to divide the gene expression of pcyt2 into three ranks 
# because not all groups are represented per the different expression level
# problem of study: sometimes no unaffected people for homozygous recessive phenotype

# glm is more robust in the sense that it does not assume normalty (homoskedasticity) 
# and works better on count data, allows for different error distributions
# here storing the results in a list, problem here again might still be the non-normalty of data
# after min-max transformation - the ancova though is actually testing the differences between models
# focuses more on the variances and not the relationships themselves (careful not to overinterpret)
robust_ancova_glm_results <- list()

# again doing a for loop for all the loci
for (locus in names(loci_subsets)) {
    locus_data <- loci_subsets[[locus]] # extracting the loci
  # converting all to factors (should in theory already be converted)
  locus_data$phenotype <- factor(locus_data$phenotype, levels = c("unaffected", "Epilepsy"))
  locus_data$genotype <- factor(locus_data$genotype)
  # again fitting a full model with the interaction term
  model_glm_full <- glm(pcyt2_mm ~ phenotype + genotype + phenotype:genotype, 
                    family = "gaussian", data = locus_data) # not really gaussian
  # fitting the reduced model without the interaction term
  model_glm_reduced <- glm(pcyt2_mm ~ phenotype + genotype, family = "gaussian", data = locus_data)
  # calculating the deviance of teh full model 
  deviance_glm_full <- deviance(model_glm_full)
  # calculating the deviance of the reduced model
  deviance_glm_reduced <- deviance(model_glm_reduced)
  # deviance for the null model (intercept-onlz) - all the deviances measure how well other models fit compared to null
  model_glm_null <- glm(pcyt2_mm ~ 1, family = "gaussian", data = locus_data)
  deviance_glm_null <- deviance(model_glm_null)
  # here, calculating the McFadden pseudo R^2 for the models
  # pseudo R^2 tells how much of the deviance is explained
  pseudo_r2_full <- 1 - (deviance_glm_full / deviance_glm_null)
  pseudo_r2_reduced <- 1 - (deviance_glm_reduced / deviance_glm_null)
  # the f2 is assessing the effect size 
  f2 <- (pseudo_r2_full - pseudo_r2_reduced) / (1 - pseudo_r2_full)
  # here, we are calculating the robust standard error again similarly to the lm ancova
  robust_se <- sqrt(diag(vcovHC(model_glm_full, type = "HC3")))
  # extracting the model coefficients for each locus
  coefficients <- summary(model_glm_full)$coefficients
  # again here computing the t statistic and p values
  robust_t_stats_glm <- coefficients[, "Estimate"] / robust_se
  robust_p_values_glm <- 2 * (1 - pnorm(abs(robust_t_stats_glm)))
  # here adjusting the p values for multiple testing using BH
  adjusted_p_values_glm <- p.adjust(robust_p_values_glm, method = "BH")
  # saving teh results per locus
  robust_ancova_results[[locus]] <- list(
    adjusted_p_values_glm = adjusted_p_values_glm,
    Pseudo_R_squared_full = pseudo_r2_full,
    Cohen_f2 = f2
  )
  # printing the results per locus
  cat("\nrobust ANCOVA glm for", locus, "\n")
  print(data.frame(Coefficients = coefficients[, "Estimate"], 
                   Robust_SE = robust_se, 
                   Robust_T = robust_t_stats_glm, 
                   Robust_P = robust_p_values_glm, 
                   Adjusted_P = adjusted_p_values_glm, 
                   Pseudo_R_Squared = pseudo_r2_full,
                   Cohen_f2 = f2))
}
# results analogous to the lm robust ANCOVA

# some diagnostic plots and why lm and glm alone are not the best models for this dataset
for (locus in names(loci_subsets)) {
  locus_data <- loci_subsets[[locus]] # again for loop for the loci extracting each
  # make sure everything categorical is a factor
  locus_data$phenotype <- factor(locus_data$phenotype, levels = c("unaffected", "Epilepsy"))
  locus_data$genotype <- factor(locus_data$genotype)
  # the LM
  lm_model <- lm(pcyt2_mm ~ phenotype + genotype + phenotype:genotype, data = locus_data)
  # diagnostic plots for the lm full model
  par(mfrow=c(2, 2))  # fit all 4 on one page
  plot(lm_model)  # plot them
  # diagnostic plots for the glm full model
  glm_model <- glm(pcyt2_mm ~ phenotype + genotype + phenotype:genotype, family = "gaussian", data = locus_data)
  # diagnostic plots for the glm model 
  plot(glm_model$residuals, main = paste("Residuals for GLM - Locus", locus), ylab = "Residuals", xlab = "Index")
  # using the predicted values, model fit plot for the glm
  glm_preds <- predict(glm_model, type = "response") 
  plot(locus_data$pcyt2_mm, glm_preds, main = paste("Fit for GLM - Locus", locus), xlab = "Observed Values", ylab = "Predicted Values")
  # for the linear model, also plotting the predicted values
  lm_preds <- predict(lm_model)
  plot(locus_data$pcyt2_mm, lm_preds, main = paste("Fit for LM - Locus", locus), xlab = "Observed Values", ylab = "Predicted Values")
}
# as can be seen, the data is not normally distributed and the expression data, although min-maxed is multimodal
# there are more complex tests that can be performed as well - below

###############################################
### power analysis for ancova (modified lm) ###
###############################################

#How would we test if the study is sufficiently powered to draw conclusions?
#What would be the sample size to guarantee a 90% power of the test?

# here, testing the power for the four loci in the ancova lm
# it is quite difficult to do power analyses for non-parametric tests or glm etc. 
# because these often need simulations (quite complex)

# Parameters for power analysis - using the pwr package and examples from the lessons
# Define effect sizes and group sizes in a list - using the Cohen f2 from above for the 4 loci
effect_sizes <- c(0.097, 0.239, 0.037, 0.0003) # loc77, loc85, loc86, loc97
group_sizes <- c(4, 6, 6, 4) # group sizes per locus, for 3 possible genotypes - 6, for 2 there are 4
sample_sizes <- seq(10, 300, by = 10) # simulation of sample sizes - for graphing

# Define alpha and power
alpha <- 0.05        # significance level at 0.05
power <- 0.9         # power of 90% (0.9)

# list to store the sample sizes
required_sample_sizes <- list()

# loop through the effect size list
for (i in 1:length(effect_sizes)) {
  # calculating the power result for each effect size
  power_result <- pwr.f2.test(u = group_sizes[i] - 1,  # numerator degrees of freedom
                              v = NULL,                # denominator degrees of freedom (sample size)
                              f2 = effect_sizes[i],    # Effect size (Cohen's f2)
                              sig.level = alpha,       # significance
                              power = power)           # power
  
  # calculating and storing the values according to the lesson example
  required_sample_sizes[[i]] <- ceiling(power_result$v / group_sizes[i])
}

# printing the required sample sizes for each effect size
required_sample_sizes 
# problem with the power values here are that the lm and 
# ancova are not the best ways to calculate the power
# also because the loc97 is non-significant, the Cohen f2 is very small 
# and therefore sample required unrealistically large - this locus should not 
# be used for investigation because it is not significant in most of the performed tests

# here, for graphing, creating a dataframe in for loops
power_data <- data.frame()

for (i in 1:length(effect_sizes)) { # for loop for all four effect sizes
  for (n in sample_sizes) { # applying all the sample size sequences to all effect sizes
    # calculating the power result - same as above but the v which here is the sample size n
    power_result <- pwr.f2.test(u = group_sizes[i] - 1, # numerator df group sizes
                                v = n - group_sizes[i], # sample sizes - denominator df
                                f2 = effect_sizes[i],   # effect sizes - cohen f2
                                sig.level = alpha)      # significance 0.05
    # here, storing the results in the dataframe using rbind
    power_data <- rbind(
      power_data, 
      data.frame(
        Locus = paste0("loc", c("77", "85", "86", "97")[i]), #pasting for concatenation
        Effect_Size = effect_sizes[i], # effect size
        Sample_Size = n, # sample size
        Power = power_result$power # power
      )
    )
  }
}

# the using ggplot  - labels, title etc
ggplot(power_data, aes(x = Sample_Size, y = Power, color = as.factor(Locus))) +
  geom_line(linewidth = 1) +  # width of the lines
  labs(
    title = "Power vs. Sample Size for lm robust ANCOVA",
    x = "sample size [n]",
    y = "power",
    color = "locus"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") + # set1 for plot consistency
  theme(
    plot.title = element_text(size = 18),  # text size adjustment (for poster purposes)
    axis.title.x = element_text(size = 14),             
    axis.title.y = element_text(size = 14),              
    axis.text.x = element_text(size = 12),               
    axis.text.y = element_text(size = 12),                
    legend.title = element_text(size = 14),               
    legend.text = element_text(size = 12)                 
  )
# again set1 for consistency of plot colours

########################################################
### more advanced techniques - GMM and random forest ###
########################################################

# problem with ML methods usually also is a low sample size even though they are robust

### GMM ###
###########

# Gaussian mixture model is an unsupervised clustering method similar to K means
# it determines the probability of a datapoint belonging somewhere - ML method
# it can model gene expression as a mixture of distributions 
# better fit for multimodalty than glm and lm

# storing the results for each locus
mixture_model_results <- list()

for (locus in names(loci_subsets)) {
  locus_data <- loci_subsets[[locus]] # looping through all the loci and seeing 
  # again making sure everything is a factor
  locus_data$phenotype <- factor(locus_data$phenotype, levels = c("unaffected", "Epilepsy"))
  locus_data$genotype <- factor(locus_data$genotype)
  # everything but the pcyt2 as a predictor - combining here
  locus_data$group <- interaction(locus_data$phenotype, locus_data$genotype)
  response <- locus_data$pcyt2_mm # extract the response variable - here the pcyt2
  # here, fitting the GMM to the response variables
  gmm_fit <- Mclust(response, G = 2:4)  # 3 components (2 to 4 components in terms of genotype)
  # summary of the mixture model
  cat("\nMixture Model Results for", locus, "\n")
  print(summary(gmm_fit))
  mixture_model_results[[locus]] <- gmm_fit # store the results somewhere
  # extracting the classification of the datapoint and add to the locus data
  locus_data$component <- gmm_fit$classification
  # printing the summary tables of components (vs genotype and phenotype)
  cat("\ncomponent distribution across genotype", locus, "\n")
  print(table(locus_data$component, locus_data$genotype)) # for genotype
  cat("\ncomponent distribution across phenotype", locus, "\n")
  print(table(locus_data$component, locus_data$phenotype)) # for phenotype
  # visualization of PCYT2 expression colored by GMM components
  # the components predicted here are the phenotypes because genotypes are graphed
  plot <- ggplot(locus_data, aes(x = genotype, y = pcyt2_mm, color = factor(component))) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    theme_minimal() +
    labs(title = paste("PCYT2 expression with GMM components for", locus),
         y = "pcyt2 expression (min-max)",
         color = "component")
  # print all the plots for all the loci
  print(plot)
}

### problem here: GMM puts all lower expression data to Epilepsy patients
# does not capture the intricacies of the dataset e.g. people who have epilepsy 
# and do not need to have the genotype or people who are unaffected but have variants
# one important insight is that epilepsy are labelled as the lower expressed pcyt2

### random forest ###
#####################

# random forest is a robust ML method that is good for any complex categorical data
# it is robust for heteroskedasticity and multimodality - here with the gene expression

rf_results <- list() # again a list for the results

for (locus in names(loci_subsets)) {
  locus_data <- loci_subsets[[locus]] # again, looping through each locus and extracting it
  
  set.seed(123)  # setting the seed as always
  # splitting the data into a train (0.7 - 70%) and test set (0.3 - 30%)
  train_index <- createDataPartition(locus_data$pcyt2_mm, p = 0.7, list = FALSE)
  train_data <- locus_data[train_index, ] # train subset
  test_data  <- locus_data[-train_index, ] # test subset - anything but the train subset
  # training the RF model on the train data, importance as true
  rf_model <- randomForest(pcyt2_mm ~ genotype + phenotype, data = train_data, importance = TRUE)
  # then prediction on the test set
  test_predictions <- predict(rf_model, newdata = test_data)
  # calculating the R squared on the test set
  rss <- sum((test_data$pcyt2_mm - test_predictions)^2)  # residual sum of squares
  tss <- sum((test_data$pcyt2_mm - mean(test_data$pcyt2_mm))^2)  # total sum of squares
  R2 <- 1 - (rss / tss) # the R2 is 1 minus the RSS divided by TSS
  # mean squared error clculation for the test set (squared - ^2)
  rmse <- sqrt(mean((test_data$pcyt2_mm - test_predictions)^2))
  # storing all these metrics per locus in a list
  rf_results[[locus]] <- list(
    model = rf_model,
    R_squared = R2,
    RMSE = rmse,
    variable_importance = rf_model$importance
  )
  # print all of the metrics per locus
  cat("\nRandom forest for", locus, "\n")
  cat("R squared:", R2, "\n")
  cat("root mean squared error:", rmse, "\n")
  print(rf_model$importance)
  # importance tells which of the two - genotype or phenotype is
  # more important per locus - can plot them by putting the two variables and
  # importance in a dataframe
  importance_df <- data.frame(
    Variable = rownames(rf_model$importance),
    Importance = rf_model$importance[, 1]
  )
  # then plot each plot per locus as bar graphs
  importance_plot <- ggplot(importance_df, aes(x = reorder(Variable, Importance), y = Importance, fill = Variable)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Set1") +
    theme(legend.position="none") +
    labs(title = paste("Random forest feature importance for", locus), 
         x = "importance", 
         y = "discrete variables") +
    # adding the R squared for each model into the graph, with rounding to 3 sig figs and to somewhere not to overlay the colors
    annotate("text", x = 1, y = max(importance_df$Importance) * 0.9, label = paste("R^2 = ", round(R2, 3)), size = 5, color = "black")
  # display all the plots
  print(importance_plot)
}

# importance tells which one - genotype or phenotype - is more important
# R^2 tells how much variation in the data is explained by the model

# final important note to code: there is likely a way more elegant method like purr
# for per locus extraction than with for loops but that's what works well in e.g. Python 
# and what I know how to do

