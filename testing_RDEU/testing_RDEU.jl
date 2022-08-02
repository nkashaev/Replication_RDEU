using Pkg
Pkg.activate(".")
using Distributed
using Statistics
using DataFrames, CSV
addprocs(7)

@everywhere begin
  using Random
  using Combinatorics
  using LinearAlgebra
  using JuMP
  using KNITRO
end

@everywhere model="RDEU"   # put "EU" or "RDEU"
println(model)


## Defining the file directories
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("Replication_AK_RDEU",tempdir1)[end]]

## Functions
@everywhere include($(rootdir)*"/testing_RDEU/functions_common_RDEU.jl")

## Common parameters
dYm=5                               # Number of varying options in menus
dYu=6                               # Number of all options
Menus=collect(powerset(vec(1:dYm))) # Menus
gindex=gindexes(Menus,dYu)          # Indices that correspond to nonzero linearly independent frequencies
U=preferencesRDEU(dYu, model)       # All preference orders consistent with RD expected utility
G=matrixcons(gindex, Menus, U)      # Marix of 0 and 1
Omegadiag=ones(size(G,1))           # Vector of weigts

## Data for all 3 frames
X=Matrix(CSV.read(rootdir*"/data/menu_choice_low.csv", DataFrame))
println("Data is ready!")
# Sample sizes
N=size(X,1)
# Smallest sample per menu
Ntau=nintaun(X)
# Tuning parameter as suggested in KS
taun=sqrt(log(Ntau)/Ntau)

## Testing
println("Testing...")
# Estimates of g for all 3 frames
ghat=estimateg(X,gindex)
etahat=kstesstat(ghat,G,Omegadiag,taun,true)
@everywhere begin
  X=$X; N=$N; gindex=$gindex; ghat=$ghat; G=$G;
  Omegadiag=$Omegadiag; taun=$taun; etahat=$etahat; Menus=$Menus
end
# Test statistic
Tn=N*kstesstat(ghat,G,Omegadiag,0.0,false)
# Bootstrap statistics
@time Boot=pmap(ksbootseed,1:1000) # use ksbootseedstable function
# Pvalue. If it is 0.0, then pvalue<0.001
pvalue=mean(Tn.<=N.*collect(Boot))

## Saving the output
CSV.write(rootdir*"/results/Tn_pval_$(model).csv", DataFrame(Tn=Tn, pvalue=pvalue))
println("Done!")
