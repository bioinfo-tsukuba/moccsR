# moccsR: Caluculate MOCCS2score in R

## Introduction


MOCCS2 is a method for detecting transcription factor recognition sequences ([Yoshitane et al., 2019](https://doi.org/10.1038/s42003-019-0522-3)) and it was implemented in Perl ([Ozaki et al., 2016](https://doi.org/10.1016/j.compbiolchem.2016.01.014); https://github.com/yuifu/moccs). This package allows MOCCS2 to be run in R.

## Installation and loading

You can install moccsR package with this command:

```{r}
## from github
remotes::install_github("https://github.com/bioinfo-tsukuba/moccsR.git")
```
You can load moccsR package with this command:
```{r setup, message = FALSE}
library(moccsR)
```

## Data preparation

You can use ChIP-seq data as fasta file.
moccsR package has test data and You can obtain its
path with this command:
```{r}
fastaPath <- system.file("count.fasta", package = "moccsR")
```

## Caluculate MOCCS2score

You can calculate MOCCS2score with this command:
```{r}
calcMOCCS2score(fastaPath, k = 6, ignoreLowerCase = TRUE)
```
You will obtain the result like this:

```
# A tibble: 12 × 6
##    kmer     AUC count MOCCS2score `p-value` `q-value`
##    <chr>  <dbl> <dbl>       <dbl>     <dbl>     <dbl>
##  1 AAAAAA  -3   17800       -7.56      1.00      1.00
##  2 CCCAAA  88.5    50       11.8       0         0   
##  3 AAAAAC  86.5    50       11.5       0         0   
##  4 AAAACC  87.5    50       11.7       0         0   
##  5 CCCCAA  89.5    50       11.9       0         0   
##  6 CCAAAA  87.5    50       11.7       0         0   
##  7 AAACCC  88.5    50       11.8       0         0   
##  8 AACCCC  89.5    50       11.9       0         0   
##  9 ACCCCC  90.5    50       12.1       0         0   
## 10 CCCCCC  91.5    50       12.2       0         0   
## 11 CCCCCA  90.5    50       12.1       0         0   
## 12 CAAAAA  86.5    50       11.5       0         0
```
The parameter `fastaPath` represents the path to fasta file you want to run with. `k` is length of mer (if you want to caluculate MOCCS2score of 6-mer, set `k` to 6).  If parameter `ignoreLowerCase` is TRUE (default), moccsR ignores lower-case characters in fasta file. This behavior assumes that repeat sequences are represented in lower case in fasta files. You can include lower-case characters setting `ignoreLowerCase` to `FALSE`.

## Reference
Yoshitane, H., Asano, Y., Sagami, A. et al. Functional D-box sequences reset the circadian clock and drive mRNA rhythms. Commun Biol 2, 300 (2019). 
DOI: https://doi.org/10.1038/s42003-019-0522-3

Ozaki H, Iwasaki W. MOCCS: Clarifying DNA-binding motif ambiguity using ChIP-Seq data. Comput Biol Chem 63 (2016).
DOI: https://doi.org/10.1016/j.compbiolchem.2016.01.014