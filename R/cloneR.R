#' Entire cloneR functionality
#'
#' This function calculates the membership groups from isolates contained within a VCF file
#' @param input.file Unzipped VCF file
#' @param snps Number of SNPs sub-sampled for ancestry calculations
#' @param subsets Number of ancestry calculations at each K value
#' @param K The K value (or range) used to calculate ancestry matrix
#' @param CPU A number of CPUs to run the parallel version of the algorithm (LEA)
#' @param alpha A numeric value corresponding to the snmf regularization parameter (LEA)
#' @param tolerance A numeric value for the tolerance error (LEA)
#' @param iterations An integer for the maximum number of iterations in SNMF algorithm (LEA)
#' @param ploidy 1 if haploid, 2 if diploid, n if n-ploid
#' @param plot TRUE/FALSE: Plots relationships between isolates visually
#' @return membership.txt File containing membership assignments at each K value
#' @return Kplot.tiff TIFF file containing network plot at specified K value
#' @export
#' @examples
#' cloneR("example.vcf", snps = 1000, subsets = 100, K = 2:20)

cloneR <- function (input.file, snps = 1000, subsets = 100, K,
                                CPU = 1, alpha = 10,
                                percentage = 0.05, iterations = 200, ploidy = 2, plot = TRUE) {
  #creates labels for ancestry
  extract_labels(input.file)
  
  #remove erroneous SNP metadata from VCF file
  x = paste0(input.file)
  exe = paste0("sed 's#:[0-9][^[:space:]]*##g;s#:[.][^[:space:]]*##g' ", x, " > ", x, ".concise.vcf")
  system(exe)

  #convert VCF file to geno file
  if(ploidy == 2){
  inter.file <- LEA::vcf2geno(paste0(input.file,".concise.vcf"), force = TRUE) #from LEA
  }
  if(ploidy == 1){
  prep_haploid2.sh <- paste('sed "s#[.]/[.]#9#g" $1 |', " awk '($0 !~ ", '"#")', "' | cut -f10- | sed 's#[.]#9#g'")
  write.table(prep_haploid2.sh, file = "prep_haploid2.sh", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "")
  exe2 <- "chmod +x ./prep_haploid2.sh"
  system(exe2)
  y <- paste0(input.file,".concise.vcf")
  exe3 <- paste0("./prep_haploid2.sh ", y, " | sed 's/\t//g' > ", y, ".geno")
  system(exe3)
  inter.file <- paste0(y, ".geno")
  }
  else{
    stop("Data need to either be haploid or diploid")}

  #convert VCF file to geno file
  inter.file <- LEA::vcf2geno(input.file, force = TRUE) #from LEA

  #make subsets of vcf file
  make_subsets(inter.file, snps = snps, subsets = subsets, ploidy = ploidy)

  #calculate Q matrices for all subsets
  make_Q(K, CPU = CPU, alpha = alpha, iterations = iterations, ploidy = ploidy)

  #find euclidean distances for all subsets Q matrices
  calc_dist(K, plot = plot)
}



