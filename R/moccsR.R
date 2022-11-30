## usethis namespace: start
#' @useDynLib moccsR, .registration = TRUE
## usethis namespace: end
## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
getvalue <- Vectorize(
    function(result, kmer, num) {
        return(result[[kmer]][num])
    },
    vectorize.args = "kmer"
)

#' Calculate MOCCS2score
#' @param fasta Input fasta file
#' @param k Length of k-mer
#' @param ignoreLowerCase Logical parameter. If TRUE, the function ignores
#' lower case character in fasta file.
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr %>%
#' @importFrom dplyr sym
#' @importFrom stats p.adjust
#' @importFrom stats pnorm
#' @return data.frame that consists of kmer, auc, MOCCS2score, p-value, and
#' q-value.
#' @examples 
#' fastaPath <- system.file("count.fasta", package = "moccsR")
#' calcMOCCS2score(fastaPath, k = 6, ignoreLowerCase = TRUE)
#' @export
calcMOCCS2score <- function(fasta, k, ignoreLowerCase = TRUE) {
    auc <- sym("AUC")
    pValue <- sym("pValue")
    qValue <- sym("qValue")
    count <- sym("count")
    kmer <- sym("kmer")
    aucList <- cppCalcMOCCS2score(fasta, k, ignoreLowerCase)
    rettb <- tibble(kmer = names(aucList)) %>%
        mutate(AUC = getvalue(aucList, !!kmer, 1)) %>%
        mutate(count = (getvalue(aucList, !!kmer, 2))) %>%
        mutate(MOCCS2score = getvalue(aucList, !!kmer, 3)) %>%
        mutate(
            pValue = 1 - pnorm(
                !!auc,
                mean = 0,
                sd = sqrt(
                    getvalue(aucList, !!kmer, 4)^2 / 12 / !!count
                )
            )
        ) %>%
        mutate(qValue = p.adjust(!!pValue)) %>%
        mutate("p-value" = !!pValue) %>%
        mutate("q-value" = !!qValue) %>%
        select(-pValue, -!!qValue)
    return(rettb)
}