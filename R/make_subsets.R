#' Make subsets of .geno file for euclidean distance calculations
#'
#' This function creates multiple smaller .geno files from larger .geno file.
#' @param input.file .geno file created by LEA::vcf2geno()
#' @param snps Number of SNPs sub-sampled for ancestry calculations
#' @param subsets Number of ancestry calculations at each K value
#' @return .geno file
#' @export
#' @examples
#' make_subsets("example.geno", snps = 1000, subsets = 100)

make_subsets <- function(input.file, subsets = 100, snps = 1000) {
  #input file
  input.file = read.geno(input.file)
  #filter out multi-allelic sites
  if(ploidy == 1){
  input.file <- input.file |> dplyr::filter(!grepl('2|3|4', x))
  }
  #check for "subsets" directory
  if (dir.exists("subsets")){
    ans = "None"
    while (ans != "y" && ans != "n") {
    ans <- readline("'subset' directory already exists. Do you wish to overwrite the directory? (y/n)")
      if (ans == "y") {
        print("Overwriting existing 'subset' directory'")
        unlink("subsets", recursive = TRUE, force = TRUE)
        dir.create("subsets") }
      else if (ans == "n") {
       print("Keeping existing 'subset' directory") 
        stop() }
      else {
        print("You need to enter y/n") }
    }
  }
  else {
    dir.create("subsets")
  }
  #randonly subsampling snps from geno file
  for (r in 1:subsets) {
      outfile <- input.file[sample(nrow(input.file), size = snps), ]
      write.table(outfile, file = paste("subsets/rep_", snps, "_", r, ".geno", sep = ""), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "")
  }
}




