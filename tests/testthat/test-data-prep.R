test_that("raw neighbour records are sorted and filtered", {
  raw <- data.frame(
    final_id_1 = c(rep("a", 5), rep("b", 5), rep("c", 4)),
    final_id_2 = paste0("n", seq_len(14)),
    distance_m = c(5, 1, 4, 2, 3, 10, 8, 9, 6, 7, 1, 2, 3, 4)
  )

  out <- reorder_distance_data(raw, min_neighbors = 5L)
  expect_equal(out$distances$final_id_1, c("a", "b"))
  expect_equal(
    unname(as.numeric(unlist(out$distances[1, paste0("distance_", 1:5)]))),
    1:5
  )
  expect_equal(
    unname(as.numeric(unlist(out$distances[2, paste0("distance_", 1:5)]))),
    6:10
  )
  expect_equal(nrow(out$metadata), 2L)
})

test_that("complete distance extraction retains row indices", {
  x <- data.frame(
    distance_1 = c(1, 2),
    distance_2 = c(2, NA),
    distance_3 = c(3, 4),
    distance_4 = c(4, 5),
    distance_5 = c(5, 6)
  )
  profiles <- extract_distances(x, 5L)
  expect_equal(nrow(profiles), 1L)
  expect_equal(attr(profiles, "row_index"), 1L)
})
