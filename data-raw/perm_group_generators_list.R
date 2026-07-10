## code to prepare `perm_group_generators_list` dataset
# See OEIS sequence A051625
# See ISSUE#21 **"BF" optimization**

get_perm_group_generators <- function(p) {
  stopifnot(1 < p, p < 19)
  message("Getting perm_group_generators for p = ", p, ". There will be 2 progress bars:")
  start_time <- Sys.time()
  
  p_factorial <- factorial(p)
  perm_group_generators <- rep(NA, p_factorial)
  
  all_perms_list <- permutations::as.cycle(permutations::allperms(p))
  all_perms_hash <- hash::hash()
  
  progressBar <- utils::txtProgressBar(min = 0, max = p_factorial, initial = 0)
  for (i in seq_along(all_perms_list)) {
    utils::setTxtProgressBar(progressBar, i)
    all_perms_hash[[as.character(all_perms_list[i])]] <- i
  }
  close(progressBar)
  
  get_perms_in_group <- function(i) {
    this_perm <- all_perms_list[i]
    this_perm_cumulated <- this_perm
    perms_in_group <- integer(0)
    
    repeat {
      this_perm_cumulated <- permutations::as.cycle(this_perm_cumulated * this_perm)
      char <- as.character(this_perm_cumulated)
      if (char == "()") break
      perms_in_group <- c(perms_in_group, all_perms_hash[[char]])
    }
    
    perms_in_group
  }
  
  perm_group_generators[1] <- TRUE
  progressBar <- utils::txtProgressBar(min = 0, max = p_factorial, initial = 0)
  
  for (i in 2:p_factorial) {
    utils::setTxtProgressBar(progressBar, i)
    if (!is.na(perm_group_generators[i])) next
    
    perm_group_generators[i] <- TRUE
    
    perms_in_group <- get_perms_in_group(i)
    if (length(perms_in_group) == 0) next
    
    # group order = (non-identity powers) + generator + identity
    group_order <- length(perms_in_group) + 2L
    
    for (j in seq_along(perms_in_group)) {
      perm_j <- perms_in_group[j] # perm_j = (all_perms_list[i]) ^ (j+1)
      # sigma^(j+1) generates the same full cyclic group as sigma iff gcd(j+1, order) = 1.
      # Mark it FALSE so it isn't treated as a new canonical generator.
      if (is.na(perm_group_generators[perm_j]) && numbers::coprime(j + 1L, group_order)) {
        perm_group_generators[perm_j] <- FALSE
      }
    }
  }
  close(progressBar)
  
  message("Time elapsed: ", format(Sys.time() - start_time))
  
  perm_group_generators
}


# Time it took on Apple M2, single core

perm_group_generators_3 <- get_perm_group_generators(3) # 0.03 secs
perm_group_generators_4 <- get_perm_group_generators(4) # 0.04 secs
perm_group_generators_5 <- get_perm_group_generators(5) # 0.2 secs
perm_group_generators_6 <- get_perm_group_generators(6) # 0.7 secs
perm_group_generators_7 <- get_perm_group_generators(7) # 5.4 secs
perm_group_generators_8 <- get_perm_group_generators(8) # 1 min
perm_group_generators_9 <- get_perm_group_generators(9) # 1 hour 27 mins
#perm_group_generators_10 <- get_perm_group_generators(10) # see ISSUE#21 **"BF" optimization**

perm_group_generators_list <- list(
  perm_group_generators_3,
  perm_group_generators_4,
  perm_group_generators_5,
  perm_group_generators_6,
  perm_group_generators_7,
  perm_group_generators_8,
  perm_group_generators_9
)

usethis::use_data(
  perm_group_generators_list,
  internal = TRUE,
  overwrite = TRUE,
  compress = "bzip2" # the smallest of "gzip", "bzip2", "xz".
)
