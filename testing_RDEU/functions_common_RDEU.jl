#This function computes the smallest number of observations per menu (needed to τ_n)
## X=is a dataset
## There are 32 menus
function nintaun(X)
    dM=maximum(X[:,1])
    Ntau=10000.0
    for i=2:dM # First menu is empty
      Ntau=minimum([Ntau,sum(X[:,1].==i)])
    end
    return Ntau
end

# Finding nonredundant indices
function gindexes(Menus,dYu)
    dM=length(Menus)
    MM=zeros(dM,dYu)
    if dYu==6
        #For RUM we add the default to every menu
        Menus=[vcat(1,Menus[i].+1) for i in eachindex(Menus)]
        Menus[1]=[]
    end
    for i in 1:dM
        MM[i,Menus[i][1:end-1]].=1.0 # The last option in the menu is dropped
    end
    # Nonzero and linearly independent frequencies in calibrated probabilities
    gindex=findall(MM.==1.0)
    return gindex
end

# Given data computes empirical frequencies for all menus and options
function frequencies(data)
    dM=maximum(data[:,1]); dY=maximum(data[:,2])-minimum(data[:,2])+1;
    F=zeros(Float64,dM,dY)
    for i in 1:size(data,1)
        F[data[i,1],data[i,2]+1]= F[data[i,1],data[i,2]+1]+1.0
    end
    # Computing sample frequencies
    P=zeros(Float64,dM,dY)
    P[1,1]=1.0
    P[2:end,:]=F[2:end,:]./sum(F[2:end,:],dims=2)
    # Throwing away P(default) since it is equal to 1-sum(P[1:5])
    P=P[:,2:end]
    return P
end

# This function translates integer menu identifiers to actual menus
function int2menu(number, menus=Menus)
    return menus[number]
end

# This function translates actual menus to integer menu identifiers
function menu2int(vector::Array{Int64,1})
    p=true; i=1;
    while p
        if int2menu(i)==vector
            p=false
        end
        i=i+1
    end
    return i-1
end

# This function computes all subsets of a given set
# Uses Combinatorics
function subsetsint(intset,menus=Menus)
    Temp=collect(powerset(int2menu(intset, menus)))
    return [menu2int(Temp[i]) for i in 1:length(Temp)]
end


# Computing G
# U is the set of preferences
# M is the set of menus
# gindexsh are coordinates of nonzero linearly independent p_pi
function matrixcons(gindexsh, M, U)
    dYu=length(U[1])
    dYm=length(M)
    if dYu==6
        M=[vcat(1,M[i].+1) for i in eachindex(M)]
        M[1]=[]
    end
    d=length(U)
    d2=length(gindexsh)
    B=zeros(d2,d)
    m1=1
    for j in 1:d
        pp=zeros(dYm,dYu)
        for i in eachindex(M)
            if length(M[i])>0 # Skipping empty menu
                for k in 1:dYu
                    if k==M[i][argmax(U[j][M[i]])]
                        pp[i,k]=1.0
                    end
                end
            end
        end
        B[:,m1]=pp[gindexsh] #picking only relevant elements, indexes that are not always zero
        m1=m1+1
    end
    return B
end

# This function computes the vector if linearly independent nozero elements of
# p_\pi and m (if needed)
function estimateg(X,gindex)
  Y=frequencies(X)
  return [1.0 .- sum(Y,dims=2) Y][gindex]
end

# This function computes the test statistic
function kstesstat(ghat,G,Omegadiag,taun,solution)
    if sum(isnan.(ghat))>0
        return -100
    end
    dr,dg=size(G)
    KS=Model(KNITRO.Optimizer)
    set_optimizer_attribute(KS,"outlev",0)
    @variable(KS,etavar[1:dg]>=taun/dg) #taun is a tuning parameter
    @objective(KS,Min,sum((sqrt(Omegadiag[r])*(ghat[r]-sum(G[r,l]*etavar[l] for l in 1:dg)))^2 for r in 1:dr))
    JuMP.optimize!(KS)
    if solution==true
        return G*value.(etavar)
    else
        return objective_value(KS)
    end
end

# This function computes the bootstrap statistic given the seed for a given frame
function ksbootseed(seed)
  ghatb=estimateg(genbootsample(X,seed),gindex)
  return kstesstat(ghatb-ghat+etahat,G,Omegadiag,taun,false)
end

#This function generates a bootstrap sample that has positive probability of the outside option for every menu
function genbootsample(Xt,seed)
    rng1=MersenneTwister(seed)
    dd=false
    Xtb=zeros(size(Xt))
    while dd==false
        Xtb=Xt[rand(rng1,1:size(Xt,1),size(Xt,1)),:]
        dd=minimum(1.0 .- sum(frequencies(Xtb),dims=2))>0.0
    end
    return Xtb
end 



#This function computes the set of RDEU and EU orders
function preferencesRDEU(dYu, model)
    U=collect(permutations(vec(1:dYu))) # All preference orders
    Lot=[[] for i=1:dYu]
    Lot[1]=[0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5]
    Lot[2]=[0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0]
    Lot[3]=[0.25, 0.0, 0.25, 0.0, 0.0, 0.25, 0.25]
    Lot[4]=[0.25, 0.2, 0.0, 0.15, 0.0, 0.0, 0.4]
    Lot[5]=[0.0, 0.2, 0.25, 0.15, 0.0, 0.25, 0.15]
    if dYu==6
        Lot[6]=[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    end
    dPi=length(Lot)
    K=length(Lot[1])
    Cump=zeros(K,dPi)
    for k in 1:K, l in 1:dPi
        Cump[k,l]=sum(Lot[l][1:k])
    end
    UCump=sort(unique(Cump))
    Index=[[] for i=1:length(UCump)]
    F1base=zeros(K,dPi)
    for i in 1:length(UCump)
        F1base[Cump.==UCump[i]].=i
    end
    B=[[] for i=1:length(U)]
    t=1
    for j in 1:length(U)
        KS=Model(KNITRO.Optimizer)
        set_optimizer_attribute(KS,"outlev",0)
        @variable(KS,uvar[1:K]>=0)
        @objective(KS,Min,uvar[1]^2);
        if model=="EU"
            M=Lot[U[j]]
            @constraint(KS,cons[l=1:length(Lot)-1],sum((M[l][k]-M[l+1][k])*uvar[k] for k in 1:length(Lot[1]))>=0.01 )
        else
            Cump1=Cump[:,U[j]]
            F1=F1base[:,U[j]]
            for i in 1:length(Index)
                Index[i]=findall(F1.==i)
            end
            @variable(KS,0.0<=fvar[1:K,1:dPi]<=1.0)
            
            #RDEU constraint
            @constraint(KS,cons1[l=1:dPi-1],(fvar[1,l]-fvar[1,l+1])*uvar[1]+sum((fvar[k,l]-fvar[k-1,l]-fvar[k,l+1]+fvar[k-1,l+1])*uvar[k] for k in 2:K)>=0.001 )
            
            #Equality constraint
            for l in 1:length(Index)
                @constraint(KS,[i=2:length(Index[l])],fvar[Index[l][1]] - fvar[Index[l][i]]==0.0)
            end
            #0 constraint
            if UCump[1]==0.0
                @constraint(KS,fvar[Index[1][1]]==0.0)
            end
            #1 constraint 
            @constraint(KS,fvar[Index[end][1]]==1.0)
            
        end
        JuMP.optimize!(KS)
        if termination_status(KS) in [MOI.OPTIMAL,MOI.LOCALLY_SOLVED]
            B[t]=U[j]
            t=t+1
        end
    end
    return B[[length(B[i])>0 for i in 1:length(B)]]
end
