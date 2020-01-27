library (data.table)
library (tidyverse)

read_full_data <- function(file_fast, file_EZ_means){
  # Merge the data generated in the fast-dm modelling, the Robust-EZ modelling and the
  # standard stats for the trials, returning a tidied tibble 

  print("FIX THIS FUNCTION")
  return(as_tibble(fread("results/full_results.csv.gz")))
  
  fast_data <- as_tibble(fread('./results/fast_dm_results.csv.gz')) %>% arrange(id, simulation)
  EZ_means_data <- as_tibble(fread('./results/means_EZ_results.csv.gz')) %>% arrange(id, simulation)

  merge(fast_data, EZ_means_data, by = c("dataset", "id", "simulation")) %>%
    arrange(subject_id, simulation)
  
}

locate_row <- function(cur_subject_index, cur_sim_index, n_sims){
  # Function that given a subject index and a simulation index returns the
  # row number of that specific simulation.
  # It works as long as the data frame or tibble is adequately sorted
  # (first by subject ids, then simulation indexes) and all of the subjects have
  # the same number of simulation per subject (n_sims)
  return(((cur_subject_index - 1)*n_sims) + cur_sim_index)
}

my_t.test <- function(selected_var, t_data, sample_size){
  # Compare the two halves of a data frame's selected variable  
  selected_data <- t_data %>% pull({{selected_var}})
  return(t.test (x = selected_data[1:floor(sample_size/2)],
                 y = selected_data[ceiling(sample_size/2):sample_size])$p.value)
  
}

full_test <- function(cur_sample, sample_size){
  # Run t test for independent samples on all of the fast and Robust parameters
  # and the mean and trimmed mean, and a Wilcox test in the median
  map(.x = c("fast_a","fast_v","fast_t0",
             "Robust_EZ_a", "Robust_EZ_v", "Robust_EZ_Ter",
             "stats_mean", "stats_trimmed_mean"),
      .f = my_t.test,
      t_data = cur_sample, sample_size = sample_size) %>%
    as.double() %>%
    c(wilcox.test(
      x = cur_sample$stats_median[1:floor(sample_size/2)],
      y = cur_sample$stats_median[ceiling(sample_size/2):sample_size]
    )$p.value)
}

random_sample <- function(my_data, sample_size, n_ids, n_sims){
  # Select a random sample of IDs, and for each IDs
  # select a random simulation
  my_data %>%
    slice(as.integer(map2(sample(n_ids, sample_size), # select sample_size number of subject ids
                          sample(1:n_sims, sample_size, replace = TRUE), # for each subject id select a simulation
                          locate_row, n_sims = n_sims))) # and find their row numbers
}

create_p_values_tibble <- function(n_samples){
  # Preallocation of a tibble of 9 columns and n_sampls rows to store the p values of
  # all the tests that will be carried out
  p_values <- matrix(NA, nrow = n_samples, ncol = 9)
  colnames(p_values) <- c("fast_a","fast_v","fast_t0",
                          "Robust_EZ_a", "Robust_EZ_v", "Robust_EZ_t0",
                          "stats_mean", "stats_trimmed_mean", "stats_median")
  return (as_tibble(p_values))
}

create_significant_simulations_tibble <- function(n_samples) {
  # Preallocation of NA values tibble that will be filled
  # with the number of significant tests later
  tibble(full = rep(NA, n_samples),
         fast = rep(NA, n_samples),
         robust  = rep(NA, n_samples))
}

is_significant <- function(p_value, alpha_level) {
  # Compare the value passed with a given alpha level to determine if it is significant
  return(p_value <= alpha_level)
}

number_of_significant_tests <- function(p_value_vector, alpha_level){
  # For a vector that contains a series of p values, calculate how many of them are
  # significant
  sum(as.logical(map(p_value_vector, is_significant, alpha_level = alpha_level)))
}

sum_significant <- function(cur_p_values, alpha_level){
  # Calculate number of significant p values for all the comparisons,
  # for the fast-dm parameters comparisons and
  # for the Robust-EZ parameters comparisons
  c(number_of_significant_tests(cur_p_values, alpha_level),
    number_of_significant_tests(cur_p_values[1:3], alpha_level),
    number_of_significant_tests(cur_p_values[4:6], alpha_level))
}