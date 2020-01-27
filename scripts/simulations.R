# Load the functions that will be used
source('./scripts/functions/simulation_functions.R')
source('./scripts/functions/RobustEZ.R')

# For reproducibiility
set.seed(2020)

# Choose how many samples will be extracted for each participant
num_of_simulations <- 10000
# Choose how many trials will be extracted per participant on each simulation
trials_per_participant <- 200

# Preprocess the rt data that will be used in the simulations
# (see the function for more info) - ADD MORE INFO HERE
rt_data <- preprocess_rt_data('input/results.utf8.txt')

# Preallocate a tibble in which the results of Robust EZ modelling
# and traditionals stats will be saved
means_EZ <- create_EZ_tibble()
# Initialize the file in which those results will be saved
initialize_file(means_EZ, './output/data/means_EZ_results.csv.gz')

# This might take a while (from seconds to weeks depending on your
# computer and the parameters selected
print(paste("Launching simulations:", Sys.time()))

# Main loop
for (cur_id in sort(unique(rt_data$id))){
  
  cur_id_formatted <- sprintf("%03d",cur_id)
  check_directory('./fast_dm/S', cur_id_formatted)
  
  # Keep only the data for our current participant and the condition 'word'
  rt_id <- rt_data %>%
    filter(condition == 'word') %>%
    filter (id == cur_id) %>%
    select(-c(id,condition))
  
  # Calculate the proportion of correct trials. If it is over 0.99, choose 0.99).
  # Then, multiply if by the number of trials per participant and simulation, and round it.
  trials_correct_response <- min(trials_per_participant*0.99,
                                 round(mean(rt_id$RESPONSE)*trials_per_participant))
  # The rest of the trials have to be incorrect responses								 
  trials_incorrect_response <- trials_per_participant - trials_correct_response
  
  # Extract that number of correct and incorrect trials num_of_simulations times,
  # and each time save those values for fast-dm modelling,
  # and calculate stats and Robust-EZ parameters
  for (cur_sim in 1:num_of_simulations){
    
    cur_sim_formatted <- sprintf("%06d",cur_sim)
      
    cur_sample <- rt_id %>%
        sample_stratified(by = RESPONSE,
                          size = c(trials_incorrect_response, trials_correct_response))
  
    means_EZ[cur_sim, ] <- calculate_means_EZ(cur_sample, cur_id, cur_sim)
    save_fast_dm_data(cur_sample, path = './fast_dm/', cur_id_formatted, cur_sim_formatted)
  
  }
  
  # Save the fast-dm config file to run later fast-dm modelling
  create_fast_dm_config_file(path = './fast_dm/ctl/', cur_id_formatted)
  # Append the traditional stats and Robust-EZ parameters results file
  save_means_EZ (means_EZ, './results/means_EZ_results.csv.gz')
  
  # Message indicating when this participant has been completed
  print(paste("Completed:", cur_id, Sys.time()))
}

# Check if fast-dm is already downloaded and available, otherwise get it
check_fast_dm()
# Model all of the simulations stored for each participant using fast-dm
# Again, this might take a while (from seconds to weeks depending on your
# computer and the parameters selected
execute_parallel_fast_dm(sort(unique(rt_data$id)))
# Merge all of the results of fast-dm modelling
fast_dm_results <- join_fast_dm('./fast_dm/logs/')
# Save the fast-dm results in a single file
fwrite(as.data.table(fast_dm_results), './output/data/fast_dm_results.csv.gz')
