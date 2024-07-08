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

# https://math.stackexchange.com/questions/266569/how-to-find-the-square-root-of-a-permutation
random_root_of_perm <- function(g_perm) {
  num_of_cycles <- length(g_perm) # including all cycles of length 1
  cycle_lengths <- sapply(g_perm, length)
  
  # odd cycles:
  odd_cycles <- ((cycle_lengths %% 2) == 1)
  which_odd_cycles <- which(odd_cycles)
  sq_odd_cycles <- lapply(
    which_odd_cycles,
    function(which_odd_cycle) {
      odd_cycle_sqrt(g_perm[[which_odd_cycle]])
    }
  )
  
  # even cycles:
  sq_even_cycles_unordered <- list()
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
        
        sq_even_cycles_unordered[[length(sq_even_cycles_unordered) + 1]] <- even_cycle_sqrt(
          g_perm[[first_cycle_index]],
          g_perm[[second_cycle_index]]
        )
        
        all_even_of_length_k <- setdiff(
          all_even_of_length_k,
          c(first_cycle_index, second_cycle_index)
        )
      }
    }
  }
  # sq_even_cycles_unordered is unsorted:
  sq_even_cycles_unordered_first_elements <- sapply(
    sq_even_cycles_unordered,
    function(single_cycle) {
      single_cycle[1]
    }
  )
  order_of_sq_even_cycles <- order(sq_even_cycles_unordered_first_elements)
  
  sq_even_cycles <- vector('list', length(sq_even_cycles_unordered))
  for (i in 1:length(sq_even_cycles_unordered)) {
    sq_even_cycles[[i]] <- sq_even_cycles_unordered[[order_of_sq_even_cycles[i]]]
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
  combined_matrix <- rbind(even_cycle1, even_cycle2)
  as.vector(combined_matrix)
}
