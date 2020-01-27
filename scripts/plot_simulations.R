# Load the functions that will be used
source('./scripts/functions/plot_functions.R')
# For reproducibility
set.seed(3030)

# Generate a number of samples in which the first half of the sample will be considered a group
# and the second half of the sample will be considered another group to then calculate t-tests
# of all of the parameters estimated to generate a H0 distribution of p-values

# What is the sample size
sample_size <- 100
# What is the number of samples that will be run
n_samples <- 1000000
# Alpha level of the tests employed
alpha_level <- 0.05

# Read the data that will be used
my_data <- read_full_data(file_fast = ('./output/data/fast_dm_results.csv.gz'),
                          file_EZ_means = ('./output/data/means_EZ_results.csv.gz'))

# Find the number of subject ids in the data tibble
n_ids <- my_data$subject_id %>% unique() %>% length()

# Find the number of simulations per participant
n_sims <- my_data$simulation %>% max()

# Preallocate the tibble to store the p values
p_values <- create_p_values_tibble(n_samples)

# Preallocate the tibble to store the number of significant p values per simulation
significant_simulations <- create_significant_simulations_tibble(n_samples)

for (cur_sample in 1:n_samples){

  if (cur_sample %% 1000 == 0) {print(cur_sample)}
  # Get the rows of this sample
  this_sample <- my_data %>% 
    random_sample(sample_size, n_ids, n_sims)
  
  # Run the t tests for all of the parameters and stats, and calculate how many of them
  # are significant
  p_values[cur_sample, ] <- full_test(this_sample, sample_size)
  significant_simulations[cur_sample, ] <- sum_significant(p_values[cur_sample, ], alpha_level)

}

# Save the p values tibble
write_csv(as.data.table(p_values), './output/data/p_values.csv.gz')

# Show the frequency of different levels of false positive for each group of parameters
for group_of_parameters in c("full", "fast", "robust"){
  significant_simulations %>% 
    group_by({{group_of_parameters}}) %>%
    count() %>%
    mutate (prop = n / n_samples)
  }

# Plot the distribution of p values for each stat or parameter
p_values %>%
  mutate(sample = row_number()) %>%
  pivot_longer(cols = -sample, names_to = "parameter", values_to = "p_val") %>%
  ggplot() +
  geom_histogram(aes(p_val, y = stat(count)/n_samples), bins = 100) +
  geom_vline(xintercept = alpha_level, color = "red") +
  geom_hline(yintercept = alpha_level / 5, color = "green") + #This works only for alpha_level = 0.05 and bins = 100, check generalization
  facet_grid(parameter ~ .)



# Check what would happen if the p values were independent and randomly distributed
base_pval <- p_values
base_significant <- significant_simulations %>%
  select(-robust) %>%
  rename (n_param_9 = full, n_param_3 = fast)

for (cur_sample in 1:n_samples){
  base_pval[cur_sample,] <- runif(9)
  base_significant[cur_sample, ] <- c(
    number_of_significant_tests(base_pval[cur_sample, ], alpha_level),
    number_of_significant_tests(base_pval[cur_sample, 1:3], alpha_level))
  }

for group_of_parameters in c("n_param_9", "n_param_3"){
base_significant %>% 
  group_by({{group_of_parameters}}) %>%
  count() %>%
  mutate (prop = n / n_samples)
  }
