library (data.table)
library (tidyverse)
library (parallel)
library (pbapply)

check_fast_dm <- function(skip = FALSE){
  # We will use the fast-dm 30_2 version to adjust the diffusion model.
  # This function checks if it is already downloaded in the correct directory.
  # If it is not, it downloads and unzips the fast-dm file
  # It lacks UNIX support yet, it should download the fast-dm source and try to compile it.
  
  if (!skip) {
      
    print("ADD UNIX SUPPORT")
    
    if (!file.exists('fast_dm/fast-dm.exe')){
      fast_dm_url <- 'https://www.psychologie.uni-heidelberg.de/ae/meth/fast-dm/Win_Binaries_30_2.zip'
      utils::download.file(fast_dm_url, './fast_dm/Win_Binaries_30_2.zip')
      utils::unzip('./fast_dm/Win_Binaries_30_2.zip', files = 'Win_Binaries_30_2/fast-dm.exe',
                 exdir = './fast_dm', junkpaths = TRUE)
    }
  }

}

create_EZ_tibble <- function(){
  # Function to preallocate a tibble that will store the values of the Robust EZ modelling
  # and the standard stats (mean, trimmed mean and median) calculated
  EZ_tibble <- matrix(NA, nrow = num_of_simulations, ncol = 11)
  colnames(EZ_tibble) <- c("dataset", "id", "simulation", "accuracy",
                               "mean", "trimmed_mean", "median", "Robust_EZ_v",
                               "Robust_EZ_a", "Robust_EZ_Ter", "Robust_EZ_p_EG")
  EZ_tibble <- as_tibble(EZ_tibble)
}

check_directory <- function(path, id_formatted){
  # Function to test if a directory exists, and if this is not the case,
  # the function will create it
  
  if (!dir.exists(paste0(path, id_formatted))){
    dir.create(paste0(path, id_formatted))
  }
  
}
  
preprocess_rt_data <- function(data_filename = './input/results.utf8.txt'){
  
  # Origin of the data: http://www.lexique.org/databases/Chronolex/results.utf8.txt
  # Ferrand, L., New, B., Brysbaert, M., Keuleers, E., Bonin, P., Méot, A., … Pallier, C. (2010). The French Lexicon
  # Project: Lexical decision data for 38,840 French words and 38,840 pseudowords. Behavior Research Methods,
  # 42(2), 488–496. https://doi.org/10.3758/BRM.42.2.488
  
  if (!file.exists(data_filename)){
    data_url <- 'http://www.lexique.org/databases/Chronolex/results.utf8.txt'
    utils::download.file(data_url, data_filename)
  }
  
  # The file is separated by tabs
  # The colums available in the file are the following:
    # subject id (id)
    # identifier of the word tested in that trial (word_id)
    # word tested in the trial (word)
    # reaction time in milliseconds (rt)
    # accuracy of the response: 1 correct response, 0 incorrect response (accuracy)
    # type of stimulus presented: mot -word-, nonmot -nonword- (condition)
    # trial index for the complete task (trial)
    # trial index of that specific category (trial_condition)
  rt_data <- read_delim(data_filename, delim = '\t',
                        col_names = c('id', 'word_id', 'word', 'rt', 'accuracy',
                                      'condition', 'trial', 'trial_condition'))
  
  # Several participants have very low accuracy levels, indicating that they probably were not paying
  # attention to the task or they had poor knowledge of the vocabulary. Those with accuracy under 70%
  # in the task will not be used for the simulations
  selected_participants <- rt_data %>% 
    group_by(id) %>%
    summarise(mean_acc = mean(accuracy)) %>%
    filter (mean_acc > 0.7) %>%
    filter (mean_acc < 0.99) %>%
    select(id) %>%
    unique()
  
  rt_data <- rt_data %>%
    filter(id %in% selected_participants$id)
  
  # Preproccesing of the rt data
  rt_data <- rt_data %>%
    # remove several trials that were incorrectly codified
    filter (condition %in% c('mot', 'nonmot')) %>%
    # change the names of the conditions to English
    mutate(condition = ifelse(condition == 'mot', 'word', 'nonword')) %>%
    # create a new variable that instead of accuracy codifies the type of response
    # as needed by fast-dm: 0 (lower threshold, nonword) or 1 (upper threshold, word)
    mutate(RESPONSE = ifelse(condition == 'nonword', 1-accuracy, accuracy)) %>%
    # remove trials without response
    filter(rt < 2000) %>%
    # Change the reaction times from seconds to milliseconds
    mutate(TIME = rt/1000) %>%
    # remove the variables that will not be used for the analysis
    select(id, condition, TIME, RESPONSE)
  
  
  return(rt_data)
}

sample_stratified <- function(cur_tibble, by, size){
  # Function for stratified sampling. Choose a tibble, a variable to stratify
  # and the size of the total sample. The total sample will contain the same
  # proportion of each value of the stratifying variable as the tibble
  
  return(cur_tibble %>%
    group_nest({{by}}, keep = TRUE) %>%
    add_column(to_sample = size) %>%
    mutate(sampled  = map2(data , to_sample, ~ sample_n(.x, .y))) %>%
    .$sampled %>%
    bind_rows())
  
}

calculate_means_EZ <- function(cur_tibble, id, sim){
  # Function that will calculate the mean, trimmed mean and median of the data passed
  # andu will model the data using the Robust EZ model (insert reference)
  # It will also add additional information to match these results with other (dataset,
  # id and simulation indexes)
    
  id_formatted <- sprintf("%03d",id)
  sim_formatted <- sprintf("%06d",sim)
  
  c(
    paste0("S", id_formatted, "_Sim", sim_formatted), #dataset
      cur_id,                     # id
      cur_sim,                    # simulation
      cur_tibble$RESPONSE %>%     # accuracy rounded to 4 decimals
        mean() %>% round(4), 
      cur_tibble %>%              # mean of correct responses rounded to 4 decimals
        filter(RESPONSE == 1) %>%
        pull(TIME) %>%
        mean() %>% round(4),
      cur_tibble %>%              # trimmed mean of correct responses rounded to 4 decimals
        filter(RESPONSE == 1) %>%
        pull(TIME) %>%
        mean(trim = 0.05) %>% round(4),
      cur_tibble %>%              # median RT of correct responses
        filter(RESPONSE == 1) %>%
        pull(TIME) %>% median(),
      cur_tibble %>%              # Robust EZ parameters rounded to 4 decimals
        filter(RESPONSE == 1) %>%
        pull(TIME) %>%
        Get.Robust.vaTer(Pc = mean(cur_tibble$RESPONSE)) %>% round(4)
  )
  
}

create_fast_dm_config_file <- function(path, id) {
  # Function that creates the config file that we will use later
  # to launch the fast-dm program created by Voss & Voss (insert link)
  # This version of the function is used to calculate only a, v and t0,
  # while the rest of the parameters have fixed values.
  # Only the subject id and the path to the dat files is changed.
  # Here the use of Chi Square -method cs- of Kolmogorov Smirnoff -method ks-
  # reduces the running time 6 fold
  
  write_lines(paste0('method ks
precision 3
set st0 0
set zr 0.5
set d 0
set szr 0
set sv 0
set p 0
#depends v condition
#format condition TIME RESPONSE
format TIME RESPONSE
load .\\fast_dm\\S', id, '\\*.dat
#save .\\fast_dm\\*.par
log .\\fast_dm\\logs\\S', id, '.log'), paste0(path, 'S', id, '.ctl'))
  
}

save_fast_dm_data <- function(cur_data, path = './fast_dm/', id_formatted, sim_formatted) {
  # Save the sample generated (cur_data) in the selected path. In that path it will look for a
  # subfolder called SNNN, where NNN is the subject id formatted.
  # The file will have the name SNNN_SimZZZZZZ, where ZZZZZZ is the simulation id formatted
  write_delim(cur_data,
              paste0(path, 'S', id_formatted,
                     '/S', id_formatted,
                     '_Sim', sim_formatted,
                     '.dat'),
              col_names = FALSE)
}

save_means_EZ <- function (cur_data, filename, append){
  # Save the file making use of the fast data.table fwrite function
  # appending the data at the end of the selected file
  fwrite(as.data.table(cur_data), filename, append = TRUE)
}
  
run_fast_dm <- function(cur_id){
  # Launch the fast-dm application for the current participant using the OS
  cur_id_formatted <- sprintf("%03d",cur_id)
  system(paste0('./fast_dm/fast-dm.exe .\\fast_dm\\ctl\\S', cur_id_formatted, '.ctl'),
         show.output.on.console = FALSE)
  return(cur_id_formatted)
}


execute_parallel_fast_dm <- function(id_vector, n_cores = 4L){
  # Launch fast-dm application in parallel for all the participant ids
  # in the id_vector after creating a cluster that has n_cores cores
  # Make sure that n_cores is an integer
  
  local_cluster <- makeCluster(n_cores)
  clusterExport(local_cluster, c("run_fast_dm"))
  completed_ids <- pblapply(id_vector, run_fast_dm, cl = local_cluster)
  stopCluster(local_cluster)
}



join_fast_dm <- function(log_dir = './fast_dm/logs/'){
  # This function goes across all of the files that are stored in the directory selected
  # and joins them in a common tibble (final_data) that is returned
  
  # Create the basic tibble that will be updated
  fast_data <- tibble()
  
  # Find all the log files available
  completed_logs <- list.files(log_dir, pattern = '*.log', full.names = TRUE)
  
  for (cur_log in completed_logs){
    # The log files have a peculiar format, they have  variable nomber of spaces in several places
    # Read the file as text, remove the spaces and then convert it to a tibble
    log_text <- read_file(cur_log)
    log_text <- str_replace_all(string = log_text, pattern = '  ', replacement = ' ')
    log_text <- str_replace_all(string = log_text, pattern = '^ ', replacement = '')
    log_text <- str_replace_all(string = log_text, pattern = '\n ', replacement = '\n')
    log_text <- str_replace_all(string = log_text, pattern = '\r ', replacement = '')
    log_tibble <- read_delim(log_text, delim = ' ', skip = 1, col_names = FALSE)

    # Join the tibbles
    fast_data <- rbind(fast_data, log_tibble)
  }
  
  # Rename the tibble as appropriate
  names(fast_data) <- c("dataset",
                        "fast_a", "fast_v", "fast_t0",
                        "fast_penalty", "fast_fit", "fast_time", "fast_method")
  
  # Change the variable types as appropiate, extract the subject id
  # and number of simulation from the dataset name and reorder the variables
  fast_data <- fast_data %>%
    mutate_at(vars("fast_a", "fast_v", "fast_t0",
                   "fast_penalty", "fast_fit", "fast_time"), as.double) %>%
    mutate(subject_id = as.integer(substring(dataset, first = 2, last = 4))) %>%
    mutate(simulation = as.integer(substring(dataset, first = 9))) %>%
    select(dataset, subject_id, simulation, fast_method, fast_time,
           fast_a, fast_v, fast_t0, fast_penalty, fast_fit)
  
}


initialize_file <- function(df_names, filename = './output/data/means_EZ_results.csv.gz'){
  print("FIX THIS, IT ADDS THE EQUIVALENT OF AN ADDITIONAL PARTICIPANT INSTEAD OF WRITING ONLY THE HEADERS")
  fwrite(as.data.table(df_names), filename, append = FALSE) 
}