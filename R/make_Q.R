#' Calculate ancestry coefficients
#'
#' This function creates individual ancestry coefficients (implemented through LEA) of subsetted .geno files. The subsetted .geno files are located in 'subsets/'.
#' @param K The K value (or range) used to calculate ancestry matrix
#' @param CPU A number of CPUs to run the parallel version of the algorithm (LEA)
#' @param alpha A numeric value corresponding to the snmf regularization parameter (LEA)
#' @param tolerance A numeric value for the tolerance error (LEA)
#' @param iterations An integer for the maximum number of iterations in SNMF algorithm (LEA)
#' @param ploidy 1 if haploid, 2 if diploid, n if n-ploid
#' @param plot TRUE/FALSE: Plots relationships between isolates visually
#' @return .Q file
#' @export
#' @examples
#' make_Q(K = 2:20, CPU = 1, alpha = 10, tolerance = 1e-05, iterations = 200, ploidy = 2)

make_Q <- function(K, CPU = 1, alpha = 10, tolerance = 1e-05, iterations = 200, ploidy = 2) {
  #create 'q_file' directory to store Q matrices
  if (dir.exists("q_files")){
    ans = "None"
    while (ans != "y" && ans != "n") {
      ans <- readline("'q_files' directory already exists. Do you wish to overwrite the directory? (y/n)")
      if (ans == "y") {
        print("Overwriting existing 'subset' directory'")
        unlink("q_files", recursive = TRUE, force = TRUE)
        dir.create("q_files") }
      else if (ans == "n") {
        print("Keeping existing 'subset' directory") }
      else {
        print("You need to enter y/n") }
    }
  }
  else {
    dir.create("q_files") }
  # make Q matrix for each K value using snmf() from LEA package
  file_list <- dir(path = "subsets", full.names = FALSE, recursive = TRUE)
  tryCatch(
    for (f in paste("subsets/", file_list, sep = "")) {
    LEA::snmf(f, K, project = 'continue', repetitions = 1,
         CPU = CPU, alpha = alpha, tolerance = tolerance, entropy = FALSE,
         percentage = 0.05, iterations = iterations, ploidy = ploidy, seed = -1)
    },
    error = function(e) {"An error has occurred. The most likely culprit is that 'subset' directory contains other file types or has been renamed. Please delete those and run again. On the command line this can be done simply with the command $rm !(*geno)"})
  #make file_list without extension
  file_list <- tools::file_path_sans_ext(file_list)
  # copy Q files to new diectory ('q_files') and remove snmf directories and classes from 'subsets'
  for (f in file_list) {
    file.copy(from = paste("subsets/", f, ".snmf/K", K, "/run1/", f, "_r1.", K, ".Q", sep = ""), to = "q_files/")
    unlink(paste("subsets/", f, ".snmf", sep = ""), recursive = TRUE)
    unlink(paste("subsets/", f, ".snmfProject", sep = ""))
  }
}

