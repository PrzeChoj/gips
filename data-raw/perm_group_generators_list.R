## code to prepare `perm_group_generators_list` dataset
# See OEIS sequence A051625
# See ISSUE#21 **"BF" optimization**

get_perm_group_generators <- function(p) {
  stopifnot(2 < p, p < 11)
  message("Getting perm_group_generators for p = ", p, ". There will be 2 progress bars:")
  start_time <- Sys.time()
  
  S <- matrix(abs(rnorm(p*p)), nrow = p) |> (\(x){x %*% t(x)})() + diag(0.1, p)
  
  p_factorial <- factorial(p)
  perm_group_generators <- rep(NA, p_factorial)
  
  all_perms_list <- permutations::allperms(p)
  unclass_all_perms_list <- unclass(all_perms_list)
  projected_matrices_list <- vector("list", p_factorial)
  projected_matrices_hash_list <- numeric(p_factorial)
  
  progressBar <- utils::txtProgressBar(min = 0, max = p_factorial, initial = 0)
  for (i in seq_along(all_perms_list)) {
    utils::setTxtProgressBar(progressBar, i)
    projected_matrices_list[[i]] <- project_matrices_cpp_impl_(list(S), unclass_all_perms_list[i,])[[1]]
    projected_matrices_hash_list[i] <- round(sum(log(projected_matrices_list[[i]])), 6)
  }
  close(progressBar)
  
  
  hash_order <- order(
    projected_matrices_hash_list,
    method = "radix"
  )
  
  sorted_hashes <- projected_matrices_hash_list[hash_order]
  
  group_starts <- c(
    1L,
    which(sorted_hashes[-1L] != sorted_hashes[-length(sorted_hashes)]) + 1L
  )
  
  group_ends <- c(
    group_starts[-1L] - 1L,
    length(sorted_hashes)
  )
  
  progressBar <- utils::txtProgressBar(
    min = 0,
    max = p_factorial,
    initial = 0
  )
  
  for (group_index in seq_along(group_starts)) {
    group_start <- group_starts[group_index]
    group_end <- group_ends[group_index]
    
    indexes_with_same_hash <- hash_order[group_start:group_end]
    
    for (position in seq_along(indexes_with_same_hash)) {
      i <- indexes_with_same_hash[position]
      
      if (!is.na(perm_group_generators[i])) {
        next
      }
      
      perm_group_generators[i] <- TRUE
      
      if (position == length(indexes_with_same_hash)) {
        next
      }
      
      candidate_indexes <-
        indexes_with_same_hash[(position + 1L):length(indexes_with_same_hash)]
      
      candidate_indexes <-
        candidate_indexes[is.na(perm_group_generators[candidate_indexes])]
      
      for (j in candidate_indexes) {
        if (
          max(abs(
            projected_matrices_list[[j]] -
            projected_matrices_list[[i]]
          )) < 1e-8
        ) {
          perm_group_generators[j] <- FALSE
        }
      }
    }
    
    utils::setTxtProgressBar(progressBar, group_end)
  }
  
  close(progressBar)
  
  message("Time elapsed: ", format(Sys.time() - start_time))
  
  perm_group_generators
}


# Time it took on Apple M2, single core

perm_group_generators_3 <- get_perm_group_generators(3) # 0.0007 secs
perm_group_generators_4 <- get_perm_group_generators(4) # 0.003 secs
perm_group_generators_5 <- get_perm_group_generators(5) # 0.01 secs
perm_group_generators_6 <- get_perm_group_generators(6) # 0.04 secs
perm_group_generators_7 <- get_perm_group_generators(7) # 0.2 secs
perm_group_generators_8 <- get_perm_group_generators(8) # 1.2 secs
perm_group_generators_9 <- get_perm_group_generators(9) # 9.6 secs
perm_group_generators_10 <- get_perm_group_generators(10) # 1.77 mins
# perm_group_generators_11 <- get_perm_group_generators(11) # Takes > 40 GB of RAM; see ISSUE#21 **"BF" optimization**

perm_group_generators_list <- list(
  perm_group_generators_3,
  perm_group_generators_4,
  perm_group_generators_5,
  perm_group_generators_6,
  perm_group_generators_7,
  perm_group_generators_8,
  perm_group_generators_9,
  perm_group_generators_10
)

stopifnot(all(c(1, 2, sapply(perm_group_generators_list, sum)) == OEIS_A051625[1:(2+length(perm_group_generators_list))]))

usethis::use_data(
  perm_group_generators_list,
  internal = TRUE,
  overwrite = TRUE,
  compress = "bzip2" # the smallest of "gzip", "bzip2", "xz".
)
