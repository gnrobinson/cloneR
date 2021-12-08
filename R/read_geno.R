#' Reads .geno file for make_subsets()
#'
#' @param input.file .geno file created by LEA::vcf2geno()
#' @return .geno data frame
#' @export

read.geno <- function(input.file) {
  x = scan(file = input.file, what = "character", skip = 0, sep ="")
  return(as.data.frame(x))
}
