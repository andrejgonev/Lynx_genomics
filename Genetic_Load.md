# Genetic load
For the genetic load study, I will subsample the lynx populations. For this study, I will include the balkan (x), carpathian (x), caucasian (x), kirov (x) and poland and norway




Summary of the first round of filtering:


| Filter                               | N of variants  | Filtered variants |
|:------------------------------------:|:--------------:|:-----------------:|
| 0. GLNexus merging (Qual >10)        |     8,331,076  |   0               |
| 1. Low complexity and repeats        |     4,165,581  |   4,165,495       |
| 2. Non-biallelic sites and INDELs    |     3,149,975  |   1,015,606       |
| 3. Invariant sites                   |     3,118,927  |   31,048          |
| 4. Quality filter (QUAL>=30)         |     2,661,404  |   457,523         |




# depth filtering

[1] "bal fail: 710"
[1] "cau fail: 652"
[1] "crp fail: 594"
[1] "kir fail: 536"
[1] "nor fail: 607"
[1] "pol fail: 615"

[1] "all fail: 449"
[1] "bal fail and others pass: 0"

after applying the filter on the filter4.vcf I went from 2,661,404 to 2,643,151 (lost 18,253). Sounds a bit fishy, to check with Enrico after the vacation!