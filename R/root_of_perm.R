# https://math.stackexchange.com/questions/266569/how-to-find-the-square-root-of-a-permutation
# Therefore: a permutation is a square if and only if the number of cycles of any even length in its disjoint cycle decomposition is even
#' @noRd
#' @examples
#' p_perm <- permutations::as.cycle("(1,3,5,2,4,6)(7,8)(9,10,11)")
#' g_perm <- gips_perm(p_perm, 16)
#' is_perm_square(g_perm) # FALSE
#' 
#' p_perm_sq <- p_perm^2
#' g_perm_sq <- gips_perm(p_perm_sq, 12)
#' is_perm_square(g_perm_sq) # TRUE
#' 
#' p_perm <- permutations::as.cycle("(2,3,4,5)(6,7,8,9)(10,11,12,13)(14,15,16,17)(18,19)(20,21)")
#' g_perm <- gips_perm(p_perm, 22)
is_perm_square <- function(g_perm) {
  cycle_lengths <- sapply(g_perm, length)
  even_cycle_lengths <- cycle_lengths[(cycle_lengths %% 2) == 0]
  all((table(even_cycle_lengths) %% 2) == 0)
}

#' https://math.stackexchange.com/questions/266569/how-to-find-the-square-root-of-a-permutation
#' The odd cycles can be sqrt in a pair and without a pair
#' @noRd
#' @examples
#' g_perm <- gips_perm("(1,2,3)(4,5,6)(7,8,9,10,11)(12,13,14,15,16)(17,18,19)(20,21,22)", 22)
#' random_root_of_perm(g_perm) # 20 different possible sqrts.
#' g_perm <- gips_perm("(1,2)(3,4,5,6)(7,8,9,10,11)(12,13,14,15,16)(17,18,19)(20,21,22)(23,24)(25,26,27,28)", 28)
#' random_root_of_perm(g_perm)
#' g_perm <- gips_perm("(1,2)(3,4,5,6)(7,8)(9,10,11,12)", 12)
#' random_root_of_perm(g_perm)
random_root_of_perm <- function(g_perm) {
  num_of_cycles <- length(g_perm) # including all cycles of length 1
  cycle_lengths <- sapply(g_perm, length)
  
  # odd cycles:
  odd_cycles <- ((cycle_lengths %% 2) == 1)
  sq_odd_cycles <- list()
  if (any(odd_cycles)) {
    which_odd_cycles <- which(odd_cycles)
    odd_cycles_lengths <- sapply(
      which_odd_cycles,
      function(which_odd_cycle) {
        length(g_perm[[which_odd_cycle]])
      }
    )
    odd_cycles_lengths_unique <- unique(odd_cycles_lengths)
    
    sq_odd_cycles <- lapply(odd_cycles_lengths_unique, function(length_of_odd_cycle) {
      sqrt_for_odd_cycle_length(g_perm[cycle_lengths == length_of_odd_cycle])
    })
    sq_odd_cycles <- unlist(sq_odd_cycles, recursive = FALSE)
    # sq_odd_cycles is unsorted:
    sq_odd_cycles <- sort_cicles(sq_odd_cycles)
  }
  
  # even cycles:
  even_cycles <- !odd_cycles
  sq_even_cycles <- list()
  if (!all(odd_cycles)) {
    table_of_even_cycles <- table(cycle_lengths[even_cycles])
    for (k in as.integer(names(table_of_even_cycles))) {
      all_even_of_length_k <- which(cycle_lengths == k)
      
      while (length(all_even_of_length_k) > 0) {
        first_cycle_index <- all_even_of_length_k[1]
        second_cycle_index <- if (length(all_even_of_length_k) > 2) {
          sample(all_even_of_length_k[-1], 1) # sample has different behaviour for single number entry
        } else {
          all_even_of_length_k[-1]
        }
        
        sq_even_cycles[[length(sq_even_cycles) + 1]] <- even_cycle_sqrt(
          g_perm[[first_cycle_index]],
          g_perm[[second_cycle_index]]
        )
        
        all_even_of_length_k <- setdiff(
          all_even_of_length_k,
          c(first_cycle_index, second_cycle_index)
        )
      }
    }
    # sq_even_cycles is unsorted:
    sq_even_cycles <- sort_cicles(sq_even_cycles)
  }
  
  # combine even and odd cycle:
  i <- 1 # index for odd
  j <- 1 # index for even
  length_combined <- length(sq_odd_cycles) + length(sq_even_cycles)
  sq_combined <- list()
  while(length(sq_combined) != length_combined) {
    next_even <- if (i > length(sq_odd_cycles)) {
      TRUE
    } else if (j > length(sq_even_cycles)) {
      FALSE
    } else {
      sq_odd_cycles[[i]][1] > sq_even_cycles[[j]][1] # smaller will be next
    }
    
    if (next_even) {
      sq_combined[[length(sq_combined) + 1]] <- sq_even_cycles[[j]]
      j <- j + 1
    } else {
      sq_combined[[length(sq_combined) + 1]] <- sq_odd_cycles[[i]]
      i <- i + 1
    }
  }
  
  attr(sq_combined, "size") <- attr(g_perm, "size")
  attr(sq_combined, "class") <- attr(g_perm, "class")
  
  sq_combined
}

sort_cicles <- function(unsorted_cicles) {
  if (length(unsorted_cicles) < 2) {
    return(unsorted_cicles)
  }
  
  cycles_unordered_first_elements <- sapply(
    unsorted_cicles,
    function(single_cycle) {
      single_cycle[1]
    }
  )
  order_of_cycles <- order(cycles_unordered_first_elements)
  
  sorted_cicles <- vector('list', length(unsorted_cicles))
  for (i in 1:length(order_of_cycles)) {
    sorted_cicles[[i]] <- unsorted_cicles[[order_of_cycles[i]]]
  }
  
  sorted_cicles
}

# We want to draw a sqrt uniformly.
# Let f(n) be a number of possible sqrts for length(odd_cycles_list) == n.
# Then, f(0) = 1; f(1) = 1; f(n) = f(n-1) + (n-1)*f(n-2)
# Moreover, the probability of a specific cycle to be without a pair is f(n-1) / f(n).
#  So, the probability of a specific cycle to be with a pair is 1 - f(n-1) / f(n) = (n-1)*f(n-2) / f(n).
#  Moreover, when a specific cycle is supposed to be in a pair,
#  it has to ba in a pair with other cycle - this other cycle can be drawn uniformly.

# odd_cycles_list <- g_perm[cycle_lengths == length_of_odd_cycle]
sqrt_for_odd_cycle_length <- function(odd_cycles_list) {
  n <- length(odd_cycles_list)
  
  if (n == 0) { # no cycle to sqrt
    return(odd_cycles_list)
  } else if (n == 1) { # single cycle to sqrt
    return(list(odd_cycle_sqrt(odd_cycles_list[[1]])))
  }
  
  # this could be not calculated in recurrence, but R would have to copy it, which is also O(n).
  number_of_sqrts <- numeric(n)
  number_of_sqrts[1] <- 1
  number_of_sqrts[2] <- 2
  if (n >= 3){
    for (i in 3:n) {
      number_of_sqrts[i] <- number_of_sqrts[i-1] + (i-1) * number_of_sqrts[i-2]
    }
  }
  
  first_will_be_without_pair <- if (n == 2) {
    sample(c(TRUE, FALSE), size = 1) # Probability equal, 1/2
  } else {
    sample(
      c(TRUE, FALSE), size = 1,
      prob = c(number_of_sqrts[n-1], (n-1) * number_of_sqrts[n-2])
    )
  }
  
  if (first_will_be_without_pair) {
    first_sqrt <- odd_cycle_sqrt(odd_cycles_list[[1]])
  } else {
    i <- if (n > 2) {
      sample(2:n, size = 1)
    } else {
      2
    }
    first_sqrt <- even_cycle_sqrt(odd_cycles_list[[1]], odd_cycles_list[[i]])
    odd_cycles_list[[i]] <- NULL # this will shorten a list
  }
  odd_cycles_list[[1]] <- NULL # this will shorten a list
  
  reccurent_list_sqrt <- sqrt_for_odd_cycle_length(odd_cycles_list)
  
  reccurent_list_sqrt[[length(reccurent_list_sqrt) + 1]] <- first_sqrt
  
  return(reccurent_list_sqrt)
}

# Computes sqrt for a single odd cycle:
# odd_cycle_sqrt(c(1, 2, 3)) == c(1, 3, 2)
# odd_cycle_sqrt(c(1, 2, 3, 4, 5)) == c(1, 4, 2, 5, 3)
odd_cycle_sqrt <- function(odd_cycle) {
  l_odd_cycle <- length(odd_cycle)
  if (l_odd_cycle == 1) {
    return(odd_cycle)
  }
  
  k <- (l_odd_cycle + 1) / 2
  sqrt_odd_cycle <- numeric(k * 2 - 1)
  
  for (i in 1:k) {
    sqrt_odd_cycle[i * 2 - 1] <- odd_cycle[i]
  }
  for (i in 1:(k-1)) {
    sqrt_odd_cycle[2*i] <- odd_cycle[k + i]
  }
  
  sqrt_odd_cycle
}

even_cycle_sqrt <- function(even_cycle1, even_cycle2) {
  # even_cycle2 can start at any place:
  n <- length(even_cycle2)
  i <- sample(n, size = 1)
  even_cycle2_shuffle <- numeric(n)
  even_cycle2_shuffle[1:(n-i+1)] <- even_cycle2[i:n]
  if (i != 1){
    even_cycle2_shuffle[(n-i+2):n] <- even_cycle2[1:(i-1)]
  }
  
  combined_matrix <- rbind(even_cycle1, even_cycle2)
  as.vector(combined_matrix)
}
