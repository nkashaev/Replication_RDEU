README Document for Reproducing Results in
==========================================
"Random rank-dependent expected utility"
=============================================
Nail Kashaev
nkashaev@uwo.ca

Victor H. Aguiar
vaguiar@uwo.ca

Software
========

A version 1.6.4 of the `Julia` programming language was used in coding the analysis files. For details about how to install `Julia` on different platforms and make `Julia` programs executable from the command line see <https://julialang.org/downloads/platform/>. After installation of `Julia 1.6.4.` run `using Pkg` and `Pkg.instantiate()` in the `Julia` terminal after setting the replication folder as the main one.

Some simulations and estimations use `KNITRO 12.3`.  
In this version if you do not have `KNITRO 12.3` then you can use `Ipopt`. Results are close and qualitative results are the same. Use Knitro when possible. 

Hardware
========

The code was run on Mac mini (M1, 2020) with 16 Gb of RAM

Content
=======

-   `data`  -- the folder contains the data used in the application.

-    `results` -- the testing results for EU and RDEU models.

-   `testing_RDEU`  -- the folder contains the analysis files to replicate all testing results.

-   `Manifest.toml` and `Project.toml`  -- toml files with all necessary `Julia` packages.



Below we describe the content of every folder.


`data`
============

-    `csv` file that contain the data low cost frame.

`results`
===========

-    `Tn_pval_MODEL.csv` -- the value of the test statistic and its p-value for MODEL in {EU, RDEU}.



`testing_RDEU`
============

-    `functions_common_RDEU.jl`-- the functions used in `testing_RDEU.jl`.


-    `testing_RDEU.jl` -- the code generates the results in Section 4.
