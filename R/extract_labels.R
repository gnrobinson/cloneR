#' Extracts isolate names from VCF file
#'
#' This function calculates the euclidean distance from ancestry matrices (Q). It assigns isolates to groups using a pairwise distance of 0 (+/- 1 SD). These membership groups are then plotted (if plot = TRUE) as a network and output as a membership table at each specified K value.
#' @param file Unzipped VCF file
#' @return List of isolate names
#' @export
#' @examples
#' extract_labels("example.vcf")

extract_labels <- function(file) {
  con = file(file, "r")
  line = readLines(con, n = 1)
  while ( substring(line[1], 1, 1) == "#" ) {
    line = readLines(con, n = 1)
    if (any(grep("#CHR", line))){
      out <- strsplit(line,"\t")
      break
    }
  }
  close(con)
  out <- unlist(out)
  out <- as.data.frame(out[c(10:length(out))])
  write.table(out, file = "label.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

