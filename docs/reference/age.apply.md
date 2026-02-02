# Age Category Specific apply Function

This function performs the split-apply-combine methodology on
Phylostrata or Divergence Strata stored within the input
PhyloExpressionSet.

This function is very useful to perform any phylostratum or
divergence-stratum specific analysis.

## Usage

``` r
age.apply(phyex_set, FUN, ..., as.list = FALSE)
```

## Arguments

- phyex_set:

  a standard PhyloExpressionSet object.

- FUN:

  a function to be performed on the corresponding expression matrix of
  each phylostratum or divergence-stratum.

- ...:

  additional arguments of FUN.

- as.list:

  a boolean value specifying whether the output format shall be a matrix
  or a list object.

## Value

Either a numeric matrix storing the return values of the applied
function for each age class or a numeric list storing the return values
of the applied function for each age class in a list.

## Details

This function uses the [`split`](https://rdrr.io/r/base/split.html)
function to subset the expression matrix into phylostratum specific
sub-matrices. Internally using
[`lapply`](https://rdrr.io/r/base/lapply.html), any function can be
performed to the sub-matrices. The return value of this function is a
numeric matrix storing the return values by `FUN` for each phylostratum
and each developmental stage s. Note that the input `FUN` must be an
function that can be applied to a matrix (e.g.,
[`colMeans`](https://rdrr.io/r/base/colSums.html)). In case you use an
anonymous function you could use `function(x) apply(x , 2 , var)` as an
example to compute the variance of each phylostratum and each
developmental stage s.

## See also

[`split`](https://rdrr.io/r/base/split.html),
[`tapply`](https://rdrr.io/r/base/tapply.html),
[`lapply`](https://rdrr.io/r/base/lapply.html)

## Author

Hajk-Georg Drost

## Examples

``` r
 
# source the example dataset
data(example_phyex_set)
 
# Example 1
# get the relative expression profiles for each phylostratum
age.apply(example_phyex_set, relative_expression)
#>                       Preglobular   Globular Early Heart  Late Heart
#> Cellular Organisms   2.605911e-01 0.12793981 0.000000000 0.026881705
#> Eukaryota            1.000000e+00 0.55293348 0.189559124 0.196201915
#> Viridiplantae        9.487069e-01 0.64837547 0.383978993 0.429882908
#> Streptophyta         4.841418e-01 0.32916458 0.135164508 0.268663521
#> Streptophytina       4.658941e-01 0.23957101 0.000000000 0.121559172
#> Embryophyta          1.000000e+00 0.66912196 0.128931740 0.000000000
#> Traceophyta          9.801590e-02 0.06042883 0.135543409 0.378021395
#> Euphyllophyta        0.000000e+00 0.04677271 0.087985328 0.168200972
#> Spermatophyta        7.810519e-01 0.58144016 0.223285879 0.337486796
#> Magnoliopsida        1.629285e-01 0.09061350 0.000000000 0.018285123
#> Mesangiospermae      4.038082e-02 0.01121083 0.000000000 0.005737602
#> Pentapetalae         3.246299e-01 0.21700156 0.144509595 0.097944252
#> Rosids               3.326239e-01 0.11927725 0.032900891 0.027334952
#> Malvids              9.283219e-05 0.00000000 0.002602139 0.010052870
#> Brassicales          5.345106e-02 0.01215026 0.000000000 0.007310340
#> Brassicaceae         1.000000e+00 0.60384609 0.241587219 0.110540287
#> Camelineae           1.000000e+00 0.33199892 0.344002183 0.000000000
#> Arabidopsis          9.861585e-01 0.16623807 0.055818498 0.000000000
#> Arabidopsis thaliana 4.013569e-01 0.15321763 0.000000000 0.116830157
#>                      Early Torpedo Late Torpedo Bent Cotyledon Mature Green
#> Cellular Organisms      0.20663298   0.21843800    0.097548168    1.0000000
#> Eukaryota               0.69151453   0.57874362    0.000000000    0.4467661
#> Viridiplantae           1.00000000   0.97389494    0.179728213    0.0000000
#> Streptophyta            0.94067019   1.00000000    0.000000000    0.6574747
#> Streptophytina          1.00000000   0.92977692    0.022839379    0.6142885
#> Embryophyta             0.32543544   0.12254727    0.132335580    0.5435090
#> Traceophyta             0.95513332   1.00000000    0.530095060    0.0000000
#> Euphyllophyta           0.41611064   0.64059690    0.608733049    1.0000000
#> Spermatophyta           1.00000000   0.93870781    0.000000000    0.8373587
#> Magnoliopsida           0.24116992   0.24928704    0.206324828    1.0000000
#> Mesangiospermae         0.03784404   0.03927685    0.005984007    1.0000000
#> Pentapetalae            0.11205103   0.06150410    0.000000000    1.0000000
#> Rosids                  0.06003207   0.03575030    0.000000000    1.0000000
#> Malvids                 0.11014167   0.08366978    0.198889083    1.0000000
#> Brassicales             0.04177697   0.05724916    0.087448340    1.0000000
#> Brassicaceae            0.21256956   0.09035905    0.000000000    0.3169490
#> Camelineae              0.61326644   0.28343817    0.061167981    0.6789578
#> Arabidopsis             0.28692584   0.14326831    0.156220962    1.0000000
#> Arabidopsis thaliana    0.39498925   0.36292280    0.314706027    1.0000000

# this is analogous to 
rel_exp_matrix(example_phyex_set)
#>     Preglobular   Globular Early Heart  Late Heart Early Torpedo Late Torpedo
#> 1  2.605911e-01 0.12793981 0.000000000 0.026881705    0.20663298   0.21843800
#> 2  1.000000e+00 0.55293348 0.189559124 0.196201915    0.69151453   0.57874362
#> 3  9.487069e-01 0.64837547 0.383978993 0.429882908    1.00000000   0.97389494
#> 4  4.841418e-01 0.32916458 0.135164508 0.268663521    0.94067019   1.00000000
#> 5  4.658941e-01 0.23957101 0.000000000 0.121559172    1.00000000   0.92977692
#> 6  1.000000e+00 0.66912196 0.128931740 0.000000000    0.32543544   0.12254727
#> 7  9.801590e-02 0.06042883 0.135543409 0.378021395    0.95513332   1.00000000
#> 8  0.000000e+00 0.04677271 0.087985328 0.168200972    0.41611064   0.64059690
#> 9  7.810519e-01 0.58144016 0.223285879 0.337486796    1.00000000   0.93870781
#> 10 1.629285e-01 0.09061350 0.000000000 0.018285123    0.24116992   0.24928704
#> 11 4.038082e-02 0.01121083 0.000000000 0.005737602    0.03784404   0.03927685
#> 12 3.246299e-01 0.21700156 0.144509595 0.097944252    0.11205103   0.06150410
#> 13 3.326239e-01 0.11927725 0.032900891 0.027334952    0.06003207   0.03575030
#> 14 9.283219e-05 0.00000000 0.002602139 0.010052870    0.11014167   0.08366978
#> 15 5.345106e-02 0.01215026 0.000000000 0.007310340    0.04177697   0.05724916
#> 16 1.000000e+00 0.60384609 0.241587219 0.110540287    0.21256956   0.09035905
#> 17 1.000000e+00 0.33199892 0.344002183 0.000000000    0.61326644   0.28343817
#> 18 9.861585e-01 0.16623807 0.055818498 0.000000000    0.28692584   0.14326831
#> 19 4.013569e-01 0.15321763 0.000000000 0.116830157    0.39498925   0.36292280
#>    Bent Cotyledon Mature Green
#> 1     0.097548168    1.0000000
#> 2     0.000000000    0.4467661
#> 3     0.179728213    0.0000000
#> 4     0.000000000    0.6574747
#> 5     0.022839379    0.6142885
#> 6     0.132335580    0.5435090
#> 7     0.530095060    0.0000000
#> 8     0.608733049    1.0000000
#> 9     0.000000000    0.8373587
#> 10    0.206324828    1.0000000
#> 11    0.005984007    1.0000000
#> 12    0.000000000    1.0000000
#> 13    0.000000000    1.0000000
#> 14    0.198889083    1.0000000
#> 15    0.087448340    1.0000000
#> 16    0.000000000    0.3169490
#> 17    0.061167981    0.6789578
#> 18    0.156220962    1.0000000
#> 19    0.314706027    1.0000000
# Example 2
# compute the mean expression profiles for each phylostratum
age.apply(example_phyex_set, colMeans)
#>                       Preglobular     Globular Early Heart Late Heart
#> Cellular Organisms    319.2256868 257.52020702  198.006382  210.51096
#> Eukaryota             140.9937811 114.14220354   92.317321   92.71630
#> Viridiplantae         167.2580207 137.35667928  111.033065  115.60331
#> Streptophyta           91.2879207  84.80962271   76.700108   82.28058
#> Streptophytina        121.2804565 108.09604554   94.139882  101.22129
#> Embryophyta           345.0889115 296.06388693  216.025821  196.92246
#> Traceophyta           111.0771983  92.59719465  129.527919  248.74430
#> Euphyllophyta         352.3710864 386.91405362  417.350739  476.59226
#> Spermatophyta         110.5848040 102.51114553   88.024950   92.64401
#> Magnoliopsida          85.7850927  77.24344506   66.540425   68.70021
#> Mesangiospermae       103.8641965  78.17351420   68.299874   73.35311
#> Pentapetalae          151.8902778 117.40807204   94.182927   79.26421
#> Rosids                138.7250813  62.58309299   31.755946   29.76950
#> Malvids                 0.1975973   0.06982345    3.651396   13.90655
#> Brassicales           350.8078278 203.10827402  159.656636  185.79981
#> Brassicaceae         1240.8340834 830.33617517  454.960581  319.16868
#> Camelineae            501.9910183 287.23392235  291.092877  180.49888
#> Arabidopsis           236.1796142 126.12962280  111.309071  103.81709
#> Arabidopsis thaliana   79.2626235  60.87854425   49.526997   58.18268
#>                      Early Torpedo Late Torpedo Bent Cotyledon Mature Green
#> Cellular Organisms       294.12595    299.61730      243.38291    663.17689
#> Eukaryota                122.46562    115.69240       80.93208    107.76561
#> Viridiplantae            172.36482    169.76577       90.69763     72.80368
#> Streptophyta             110.37154    112.85162       71.05001     98.53351
#> Streptophytina           152.39469    148.30386       95.47039    129.92514
#> Embryophyta              245.14108    215.07986      216.53016    277.45226
#> Traceophyta              532.48631    554.54540      323.51261     62.88684
#> Euphyllophyta            659.68056    825.47001      801.93765   1090.89936
#> Spermatophyta            119.44055    116.96148       78.99375    112.86223
#> Magnoliopsida             95.02676     95.98553       90.91095    184.65770
#> Mesangiospermae          101.63000    102.89190       73.57013    949.02303
#> Pentapetalae              83.78378     67.58943       47.88460    368.26696
#> Rosids                    41.43888     32.77288       20.01384    376.90719
#> Malvids                  151.66834    115.23254      273.81990   1376.46550
#> Brassicales              309.05905    364.39056      472.38859   3735.84753
#> Brassicaceae             424.89225    298.25672      204.62594    533.05112
#> Camelineae               377.65922    271.62202      200.16390    398.77848
#> Arabidopsis              142.32838    123.04661      124.78512    238.03743
#> Arabidopsis thaliana      78.79086     76.41513       72.84285    123.61474

# Example 3
# compute the variance profiles for each phylostratum
age.apply(example_phyex_set, function(x) apply(x , 2 , var))
#>                       Preglobular    Globular  Early Heart  Late Heart
#> Cellular Organisms     375769.502  213691.418   93456.6263  111376.489
#> Eukaryota               18957.667    9246.500    4881.1656    5275.564
#> Viridiplantae           19765.672    9955.686    5088.6714    5898.013
#> Streptophyta             8988.070    5289.605    2911.2031    4200.383
#> Streptophytina          14051.461    9938.314    5196.7120    5942.141
#> Embryophyta           1682940.617  865145.071  456717.0145  340113.846
#> Traceophyta             11741.175    6082.060   27027.1280  232242.024
#> Euphyllophyta          305880.389  397342.608  519645.4590  678605.730
#> Spermatophyta            7915.181    6101.064    3870.7790    3901.142
#> Magnoliopsida            5305.453    4352.801    3774.8536    4394.796
#> Mesangiospermae          7492.982    5565.843    6237.1524   12496.953
#> Pentapetalae            41129.591   23603.672   16333.7718   10152.252
#> Rosids                  53859.055    5195.222    1230.0765    1446.926
#> Malvids                        NA          NA           NA          NA
#> Brassicales            316537.394   88058.415   47820.5856   53928.058
#> Brassicaceae         14885810.615 9941830.931 2669040.4130 1089242.891
#> Camelineae             307467.949   53236.390  115819.8364   16224.737
#> Arabidopsis             76646.093   16854.057   12918.5441   10053.409
#> Arabidopsis thaliana     1126.399    1146.464     999.1358    1331.448
#>                      Early Torpedo Late Torpedo Bent Cotyledon Mature Green
#> Cellular Organisms      388410.149  472014.5965    397835.8529 63911898.395
#> Eukaryota                13465.408   15107.5394     11822.8106   152722.168
#> Viridiplantae            16364.920   19284.1470      5212.9363     3250.544
#> Streptophyta             12769.289   15878.6702      2984.0422    27031.826
#> Streptophytina           16402.149   13306.9785      5833.6976    29315.188
#> Embryophyta             412912.826  295668.0230    790522.8496   414144.210
#> Traceophyta            1364256.658 1559803.5400    506744.6395     4122.238
#> Euphyllophyta          1257597.946 2173425.2889   2266699.0105  4675378.101
#> Spermatophyta             6605.489    5854.6060      4676.2139    56952.986
#> Magnoliopsida             8062.390    8888.6519      5562.2734    51672.097
#> Mesangiospermae          30301.931   42528.5521      9445.6795  8819608.027
#> Pentapetalae              9174.203    4831.4157       867.6210   649519.677
#> Rosids                    3161.268    1742.4209       261.2453   193249.936
#> Malvids                         NA           NA             NA           NA
#> Brassicales              76336.441  145638.0550    275410.3640 41379294.447
#> Brassicaceae           1578566.326  369065.8774    100803.7193  1284948.737
#> Camelineae              144921.918   63440.1248     96206.7819   338162.950
#> Arabidopsis              20162.613   12372.7956     20081.3084   100005.445
#> Arabidopsis thaliana      1192.349     580.5684      2527.3283    21297.735

# Example 4
# compute the range for each phylostratum
# Note: in this case, the range() function returns 2 values for each phylostratum
# and each developmental stage, hence one should use the argument 'as.list = TRUE'
# to make sure that the results are returned properly 
age.apply(example_phyex_set, function(x) apply(x , 2 , range), as.list = TRUE)
#> $`Cellular Organisms`
#>      Preglobular Globular Early Heart Late Heart Early Torpedo Late Torpedo
#> [1,]       0.000    0.000       0.000      0.000          0.00 3.544692e-02
#> [2,]    5564.301 4519.245    2919.754   4556.479      14832.43 1.875515e+04
#>      Bent Cotyledon Mature Green
#> [1,]           0.00          0.0
#> [2,]       15615.46     232540.6
#> 
#> $Eukaryota
#>      Preglobular     Globular  Early Heart Late Heart Early Torpedo
#> [1,]      0.0000   0.00668537   0.02359617     0.0000        0.0000
#> [2,]    892.8625 719.63352333 500.92871828   502.5445      836.5676
#>      Late Torpedo Bent Cotyledon Mature Green
#> [1,] 1.913363e-02          0.000        0.000
#> [2,] 1.066288e+03       1076.485     5007.653
#> 
#> $Viridiplantae
#>      Preglobular  Globular Early Heart Late Heart Early Torpedo Late Torpedo
#> [1,]    5.170611  10.58502     18.8333   26.27958      52.45156     48.44215
#> [2,]  543.351669 416.02187    290.4093  325.57187     508.29134    671.60642
#>      Bent Cotyledon Mature Green
#> [1,]       19.67485     11.01641
#> [2,]      339.55379    230.36384
#> 
#> $Streptophyta
#>      Preglobular     Globular Early Heart  Late Heart Early Torpedo
#> [1,]   0.1112127   0.02726316   0.2400369   0.4643987      1.437311
#> [2,] 429.2365901 298.91980215 249.1829693 368.5563000    669.422226
#>      Late Torpedo Bent Cotyledon Mature Green
#> [1,]     3.658204       20.90519    0.6574457
#> [2,]   766.647581      338.46027  742.8222730
#> 
#> $Streptophytina
#>       Preglobular     Globular Early Heart Late Heart Early Torpedo
#> [1,]   0.06808085   0.01162535    8.607448   21.13225      24.16181
#> [2,] 592.75789921 455.82413441  281.509058  280.86175     455.68970
#>      Late Torpedo Bent Cotyledon Mature Green
#> [1,]     35.45795        27.6913     5.139044
#> [2,]    444.17275       315.3078   848.215171
#> 
#> $Embryophyta
#>       Preglobular     Globular  Early Heart   Late Heart Early Torpedo
#> [1,] 7.526608e-02 3.427276e-02 1.299371e-02 6.769752e-02     0.2958471
#> [2,] 9.981932e+03 6.858560e+03 5.181343e+03 4.485187e+03  4853.5604685
#>      Late Torpedo Bent Cotyledon Mature Green
#> [1,]    0.1505351       3.776421     0.191978
#> [2,] 4100.7403444    6952.029473  3766.979231
#> 
#> $Traceophyta
#>      Preglobular  Globular Early Heart Late Heart Early Torpedo Late Torpedo
#> [1,]     28.6467  33.54971    38.44432   40.84347      43.16473      48.1074
#> [2,]    342.2688 259.05460   489.39943 1338.92886    3176.31691    3384.2641
#>      Bent Cotyledon Mature Green
#> [1,]        33.7693    0.6672493
#> [2,]      1937.1539  189.5959594
#> 
#> $Euphyllophyta
#>       Preglobular     Globular Early Heart Late Heart Early Torpedo
#> [1,]    0.1068995 2.429332e-02    1.783572   25.50669      54.57407
#> [2,] 1323.0953155 1.491979e+03 1689.893350 1933.14094    2652.15309
#>      Late Torpedo Bent Cotyledon Mature Green
#> [1,]     55.45201       26.89962     3.153947
#> [2,]   3453.34947     3489.60651  4949.764295
#> 
#> $Spermatophyta
#>       Preglobular     Globular Early Heart   Late Heart Early Torpedo
#> [1,]   0.01479277   0.05749413   0.2470909   0.08216436     0.4607543
#> [2,] 396.58922233 330.25820591 295.1633296 316.74962222   353.8924251
#>      Late Torpedo Bent Cotyledon Mature Green
#> [1,]     1.425323        4.34375     1.667889
#> [2,]   361.639571      308.00956  1321.199400
#> 
#> $Magnoliopsida
#>      Preglobular Globular Early Heart Late Heart Early Torpedo Late Torpedo
#> [1,]      0.0000   0.0000       0.000     0.0000    0.01505505    0.1552744
#> [2,]    227.9104 208.5809     261.676   301.9089  416.19400915  398.7333107
#>      Bent Cotyledon Mature Green
#> [1,]       2.935162     1.141558
#> [2,]     288.410357   923.185945
#> 
#> $Mesangiospermae
#>      Preglobular Globular Early Heart   Late Heart Early Torpedo Late Torpedo
#> [1,]   0.1096234   0.0000      0.0000   0.09332517     0.5496438     1.451886
#> [2,] 292.5547037 267.7522    361.8584 542.73093229   821.1595351   988.818290
#>      Bent Cotyledon Mature Green
#> [1,]       2.128736     1.144999
#> [2,]     484.976400 12960.780277
#> 
#> $Pentapetalae
#>       Preglobular     Globular  Early Heart   Late Heart Early Torpedo
#> [1,]   0.05261297   0.01197911   0.08091946   0.08694481     0.0982695
#> [2,] 635.42989937 495.34699334 417.37091062 317.80365833   272.6756148
#>      Late Torpedo Bent Cotyledon Mature Green
#> [1,]    0.1884891     0.06334834 7.368801e-02
#> [2,]  209.1104428    87.99405005 2.482625e+03
#> 
#> $Rosids
#>      Preglobular Globular Early Heart Late Heart Early Torpedo Late Torpedo
#> [1,]      0.0000   0.0000     0.00000    0.00000     0.1432673    0.1435256
#> [2,]    642.4961 185.3317    91.70859   98.72456   147.6228087  103.7562025
#>      Bent Cotyledon Mature Green
#> [1,]       2.773028    0.2618479
#> [2,]      41.354124 1199.1175298
#> 
#> $Malvids
#>      Preglobular   Globular Early Heart Late Heart Early Torpedo Late Torpedo
#> [1,]   0.1975973 0.06982345    3.651396   13.90655      151.6683     115.2325
#> [2,]   0.1975973 0.06982345    3.651396   13.90655      151.6683     115.2325
#>      Bent Cotyledon Mature Green
#> [1,]       273.8199     1376.466
#> [2,]       273.8199     1376.466
#> 
#> $Brassicales
#>      Preglobular     Globular Early Heart Late Heart Early Torpedo Late Torpedo
#> [1,]    0.108494   0.01235082   0.4582448   2.966973       23.2142      44.6563
#> [2,] 1185.025172 632.57525785 475.6846959 524.167552      550.5373     878.1621
#>      Bent Cotyledon Mature Green
#> [1,]       48.21886     9.177722
#> [2,]     1211.64171 13326.382738
#> 
#> $Brassicaceae
#>       Preglobular Globular Early Heart   Late Heart Early Torpedo Late Torpedo
#> [1,] 2.576976e-02     0.00       0.000 9.689722e-03  1.217018e-02 7.805728e-02
#> [2,] 1.749544e+04 15175.64    7916.273 5.056297e+03  6.073790e+03 2.698556e+03
#>      Bent Cotyledon Mature Green
#> [1,]      0.8951774    0.1766512
#> [2,]   1066.8117527 5174.1424821
#> 
#> $Camelineae
#>      Preglobular Globular Early Heart Late Heart Early Torpedo Late Torpedo
#> [1,]    103.8834 105.0861    74.49384   77.86892      123.7724     74.73676
#> [2,]   1322.0438 620.8688   797.86782  359.76975      932.8575    626.21049
#>      Bent Cotyledon Mature Green
#> [1,]       28.42762     19.43867
#> [2,]      664.64225   1256.45665
#> 
#> $Arabidopsis
#>       Preglobular Globular Early Heart  Late Heart Early Torpedo Late Torpedo
#> [1,]   0.06434192   0.0000   0.1530719   0.2052625      1.515649      5.95713
#> [2,] 841.28115092 414.7169 412.3216774 351.0940496    571.246363    445.54896
#>      Bent Cotyledon Mature Green
#> [1,]       21.84301     30.71553
#> [2,]      455.03646   1117.21434
#> 
#> $`Arabidopsis thaliana`
#>      Preglobular  Globular Early Heart Late Heart Early Torpedo Late Torpedo
#> [1,]    27.47878  24.27659    4.510751   5.873973      42.98079     37.08141
#> [2,]   118.24917 118.91433   97.911597 115.962891     127.05882    108.71476
#>      Bent Cotyledon Mature Green
#> [1,]       22.26581      10.4845
#> [2,]      162.11333     398.9621
#> 
```
