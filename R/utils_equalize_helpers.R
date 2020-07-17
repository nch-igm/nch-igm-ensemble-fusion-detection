#' Get the ordered fusions given two equalized genes.
#'
#' @param Gene1 A character vector of the equalized first fusion gene partners.
#' @param Gene2 A character vector of the equalized second fusion gene partners.
#' @return A character vector with the ordered fusions 'Gene1->Gene2'.
#' @examples
#' get_ordered(star_fusion$Gene1, star_fusion$Gene2)
get_ordered <- function(Gene1, Gene2) {
  ordered_fusion <- paste(Gene1, Gene2, sep = ">>")
  return(ordered_fusion)
}

#' Get the unordered fusions given two equalized genes.
#'
#' @param Gene1 A character vector of the equalized first fusion gene partners.
#' @param Gene2 A character vector of the equalized second fusion gene partners.
#' @return A character vector with the unordered fusions 'GeneA+GeneB'.
#'   'GeneA' < 'GeneB' evaluates to TRUE so that 'GeneA->GeneB' and
#'   'GeneB->GeneA' are converted to the same unordered fusion.
#' @examples
#' get_unordered(Gene1, Gene2)
get_unordered <- function(Gene1, Gene2) {
  unordered_fusion <- purrr::map2_chr(
    Gene1,
    Gene2,
    function(g1, g2) {
      if (g1 < g2) {
        result <- paste(g1, g2, sep = "+")
      } else {
        result <- paste(g2, g1, sep = "+")
      }
      return(result)
    }
  )
  return(unordered_fusion)
}

#' Convert chromosome records to a standard format.
#'
#' Convert chromosome records to a standard format. Common chromosome formats
#' include 'chrX' and 'X'. This function will use the format 'X', so that a
#' record of 'chrX' -> 'X'.
#'
#' @param chrom A character vector of the chromosomes.
#' @return A character vector with the chromosomes in a standard format.
#' @examples
#' equalize_chrom("chr3")
#' equalize_chrom("3")
#' equalize_chrom(c("chr3", "3"))
equalize_chrom <- function(chrom) {
  chr_opt <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
                   "12", "13", "14", "15", "16", "17", "18", "19", "20", "21",
                   "22", "X", "Y", "M")
  pattern <- paste0("(", stringr::str_c(chr_opt, collapse = "|"), ")(?:_.+)?")
  chrom <- purrr::map_chr(
    chrom,
    function(chr) {
      if ((nchar(chr) >= 3) & (tolower(substr(chr, 1, 3)) == "chr")) {
        chr <- substr(chr, 4, nchar(chr))
      }
      if (stringr::str_detect(chr, pattern) == FALSE) {
        print(paste0("Chr '", chr, "' is not a valid chr (equalize_chrom)."))
        stop()
      }
      return(chr)
    }
  )
  return(chrom)
}
