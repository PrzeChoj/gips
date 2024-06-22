## code to prepare `perm_group_generators_list` dataset
# See OEIS sequence A051625
# See ISSUE#21 **"BF" optimization**

get_perm_group_generators <- function(p) {
  stopifnot(1 < p, p < 19)
  cat(paste0("Get perm_group_generators for p = ", p, ". There will be 2 progress bars:\n"))
  start_time <- Sys.time()

  p_factorial <- prod(1:p)
  perm_group_generators <- rep(NA, p_factorial)
  
  progressBar <- utils::txtProgressBar(min = 0, max = prod(1:p), initial = 1)

  all_perms_list <- permutations::as.cycle(permutations::allperms(p))
  all_perms_hash <- hash::hash()

  for (i in 1:length(all_perms_list)) {
    utils::setTxtProgressBar(progressBar, i)
    all_perms_hash[[as.character(all_perms_list[i])]] <- i
  }
  close(progressBar)

  get_perms_in_group <- function(i) {
    perms_in_group <- numeric(0)

    this_perm <- all_perms_list[i]
    this_perm_cumulated <- all_perms_list[i]
    while (TRUE) {
      this_perm_cumulated <- permutations::as.cycle(this_perm_cumulated * this_perm)
      this_perm_cumulated_char <- as.character(this_perm_cumulated)

      if (this_perm_cumulated_char == "()") {
        return(perms_in_group)
      }

      perms_in_group[length(perms_in_group) + 1] <- all_perms_hash[[this_perm_cumulated_char]]
    }
  }

  progressBar <- utils::txtProgressBar(min = 0, max = prod(1:p), initial = 1)

  perm_group_generators[1] <- TRUE
  perm_group_generators_where <- c(1)

  for (i in 2:length(perm_group_generators)) {
    utils::setTxtProgressBar(progressBar, i)
    if (!is.na(perm_group_generators[i])) {
      next
    }

    perm_group_generators[i] <- TRUE

    perms_in_group <- get_perms_in_group(i)
    num_perms_in_group_i <- length(perms_in_group) + 2 # + the generator + id
    if (length(perms_in_group) == 0) {
      next
    }
    for (j in 1:length(perms_in_group)) {
      perm_j <- perms_in_group[j] # perm_j = (all_perms_list[i]) ^ (j+1)
      if (is.na(perm_group_generators[perm_j]) && numbers::coprime(j + 1, num_perms_in_group_i)) {
        perm_group_generators[perm_j] <- FALSE
      }
    }
  }
  close(progressBar)

  end_time <- Sys.time()

  print(end_time - start_time)

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
