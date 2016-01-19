module NeutralEvolution

################################################################################
################################################################################

#initialize required packages

using Distributions
using Gadfly
using DataFrames
using GLM
using PyCall
using Colors

################################################################################
################################################################################

export new_pop,setupprob,newGenerationGrow2,tumourgrow,tumourgrow_nonzero,setID,setcelltype,allele_frequency,analysedata,pop_expand,idealbranching,idealsampling,fullmodel_sampling,multsample,tumourgrow_sel,idealsampling_mutation,idealsampling_selection,func,mixing,prob,popstruct,simstat

################################################################################
################################################################################

#type definitions
type prob
    # type for offspring distributions
    pu::Array{Float64,1}
    pv::Array{Float64,1}
end

type popstruct
    #Type for a population, includes an array of "alive" cells, the number of mutations each cell has acquired
    # the celltype (u(neutral) vs v (fitter)) and the parent id of the cell.
    alive::Array{Int,1}
    mutations::Array{Int,1}
    celltype::Array{Int,1}
    parentid::Array{Int,1}
end

type simstat
    #struct for the statistics that may be useful from a single simulation
    numU::Int64
    popsize::Int64
    N_generation::Int64
    popsize_event::Int64
    selection_gen::Int64
    mutantid::Int64
end

################################################################################
################################################################################

#function definitions
function new_pop(N,lambda) #new population size given by sampling from Poisson distribution
    w_fit=1
    rand(Poisson(N*(1+w_fit*lambda)))
end

function setupprob(s)
    #set as our baseline that E[offspring prob dist] = 1.2
    #p0=0, p1=0.8, p2= 0.2

    #u is our normal population, v is the population that has fitness advantage
    pu=[0,0.8,0.2]

    #We define our selective advantage as 1+S=E[fit]/E[normal] to get the required fitness advantage s,
    # we calculate rho which we use to modify the original offpsring probability distribution

    p1 = pu[2]
    p2 = pu[3]

    ρ = (s*p1 + 2*s*p2)/p1

    pv=[0,pu[2]-pu[2]*ρ,pu[3]+pu[2]*ρ]

    sum(abs(pv)) == 1.0 || error("Offspring probability distribution must sum to 1 and have no negative elements, s cannot be greater than 0.25")

    return prob(pu,pv)
end

function newGenerationGrow2(N,celltype,p)
    # takes a population of pop individuals which grows according to their offspring probability distributions
    # and samples from these distributions to get daughter cells for each cell.

    #celltype refers to the type, type 0 indicating normal, type 1 have a fittness advantage.

    #denote the composition of generation i+1 with vi
    newpopulation=zeros(Int,N)
    for i in 1:N
        newpopulation[i]=wsample([0,1,2],p.pu*(1-celltype[i])+p.pv*celltype[i])
    end

    return newpopulation
end

function tumourgrow(maxpopsize::Int64,p,mu,selectiongen::Int64,lambda)
    # Tumour grows to maxpopsize given offspring probability distribution (p) and mutation rate mu
    #
    # The simulation of tumour growth is specified by a vector which we call pop
    # and a vector called parent ID. alive designates which crypts are still "alive"
    # 1 being alive 0 being dead and parent ID designates the parents of each cell.
    # The ID simply refers to the index of the vector. For example if we grow to a size of
    # 4 crypts we might get a output that looks as follows. pop = [0,0,0,1,1,1,1],
    # and parentid = [0,1,1,2,2,3,3]. This would correspond to an idealised branching
    # process whereby 1 crypts splits into 2 which splits into 4. We can
    # then calculate the allele frequency given that all mutations are passed onto daughter
    # crypts.

    #population begins with 1 cell, parentid is 0, number of mutations is zero and celltype is 0:
    alive=[1]
    mutations=[0]
    parentid=[0]
    celltype=[0]

    population = popstruct(alive,mutations,celltype,parentid)

    #define a variable selectionevent that increases by 1 at each generation, then use this
    #to specofy the generation at which a mutant allele with fitness advantage enters the population.
    selectionevent = 1
    popsize_event=0
    mutantid = 0

    while sum(population.alive)<maxpopsize
        if sum(population.alive) == 0
            break
        end

        #generate next generation given current generation and probabilities
        popnew = newGenerationGrow2( sum(population.alive) , population.celltype[end-sum(population.alive)+1:end] , p )

        #set parent ID of new population
        x = setID( [ zeros(Int,length(population.alive) - sum(population.alive)) ; popnew ] )

        # append new parentID onto previous parent IDs
        population.parentid = [population.parentid;x]

        #add new population onto old one. Remembering that this simply
        #corresponds to alive (1) or dead (0)
        population.alive = [zeros(Int,length(population.alive)) ; ones(Int,sum(popnew)) ]

        #set celltype of new population
        population.celltype = setcelltype(population.parentid,population.celltype)

        #increase selection event by 1
        selectionevent = selectionevent + 1

        #introduce fitter cell at correct generation
        if selectionevent == selectiongen
            population.celltype[end] = 1
            popsize_event = sum(population.alive)
            mutantid = population.parentid[end]
        end

    end

    # We now calculate the number of mutations that happen at each divison.
    # This can be done at the end.
    population.mutations = rand(Poisson(mu),length(population.alive))

    simstats = simstat(sum(population.celltype),sum(population.alive),selectionevent,popsize_event,selectiongen,mutantid)

    return population , simstats
end

function tumourgrow_nonzero(maxpopsize,p,mu,selectiongen,lambda)

    #often a tumour will not grow as 0 is an abosrbing state, function that calls tumourgrow until we get a tumour that has
    #grown to sufficent size
    alive=[1]
    mutations=[0]
    parentid=[0]
    celltype=[0]

    population = popstruct(alive,mutations,celltype,parentid)
    S = simstat(0,0,0,0,0,0)
    while sum(population.alive) < maxpopsize
        population,S = tumourgrow(maxpopsize,p,mu,selectiongen,lambda)
    end

    return population,S

end

function setID(vec::Array{Int,1}) #get id for parent of new mutants
    new_vec=zeros(Int,sum(vec))
    k=1

    for i = 1:length(vec)
        for j = 1:vec[i]
            new_vec[k] = i
            k = k + 1
        end
    end
    return new_vec
end

function setcelltype(parentID,celltype)
    #cells inheret the type of their parent

    x=findin(celltype,[1])
    y=Array(Int64,1)
    new_celltype = zeros(Int64,length(parentID))

    for i in x
        y=findin(parentID,i)
        new_celltype[y] = 1
    end

    return new_celltype
end

function allele_frequency(population)

    #reconstruct frequency distribution from knowledge of Parent IDs. Traces back through the tree to get allele_frequency of
    #all mutations. Uses pop struct

    maxID=maximum(population.parentid)
    numalive=sum(population.alive)
    counter=length(population.parentid)

    while maxID > 0

        index=Array(Int64,0)

        if population.parentid[counter-1]==maxID
            push!(index,counter-1)
        end

        if population.parentid[counter]==maxID
            push!(index,counter)
        end

        sumID = 0

        for i = 1 : length(index)
            sumID = sumID + population.alive[index[i]]
            counter=counter-1
        end

        population.alive[maxID] = population.alive[maxID] + sumID
        maxID = maxID - 1

    end

    return population,numalive
end

function allele_frequency(pop,parent_ID)
    #reconstruct frequency distribution from knowledge of Parent IDs. Traces back through the tree to get allele_frequency of
    #all mutations. Uses popalive vector and parent_ID vector

    popcopy = copy(pop)
    maxID=maximum(parent_ID)
    numalive=sum(pop)

    while maxID > 0

        index = findin(parent_ID,maxID)
        sumID = 0

        for i = 1 : length(index)
            sumID = sumID + popcopy[index[i]]
        end

        popcopy[maxID] = popcopy[maxID] + sumID
        maxID = maxID - 1

    end

    return popcopy,numalive
end

function analysedata(VAF,max_range,min_range)

    #fit model and data that can be used to plot histogram

    v=min_range
    u=max_range

    steps=u:-0.001:v
    cumsum=Array(Int64,0)
    v=Array(Float64,0)

    for i in steps
        push!(cumsum,length(filter((x)->x>=i,VAF)))
        push!(v,i)
    end
    cumsum=cumsum-cumsum[1]

    DFcumsum = DataFrame(cumsum=map(Float64,cumsum),v=v)
    DFcumsum[:invVAF] = 1./DFcumsum[:v] - 1./u

    x,y=hist(VAF,0.0:0.01:1)

    #fit constrained fit using GLM fit function
    lmfit=fit(LinearModel, cumsum ~ invVAF + 0 , DFcumsum)
    DFcumsum[:prediction] = predict(lmfit)
    DFcumsum[:residuals] = residuals(lmfit)

    #calculate R^2 value
    rsq = 1-(sum(residuals(lmfit).^2)/sum((DFcumsum[:cumsum]-0).^2))

    return DFcumsum,lmfit,rsq
end

function pop_expand(population)
    # adjust the population vector so that
    # we include the mutations, so transforms array of cells into array of alleles.

    pop_new=Array(Int64,sum(population.mutations))
    k=1

    for i=1:length(population.alive)

        for j=1:population.mutations[i]

            pop_new[k]=population.alive[i]
            k=k+1

        end

    end

    population.alive=pop_new

    return population
end

function pop_expand(pop,muts)
    # adjust the population vector so that
    # we include the mutations

    pop_new=Array(Float64,sum(muts))
    k=1

    for i=1:length(pop)

        for j=1:muts[i]

            pop_new[k]=pop[i]
            k=k+1

        end

    end

    return pop_new
end

function idealbranching(mu)
    #calculate popstruct given that all cells divide into 2 daughter cells at each generation

    alive=zeros(Int64,2047)
    alive[end-1023:end]=1

    parentid = zeros(Int64,2047)
    parentid[2:2047] = int(floor([1:0.5:(2047)/2]))

    mutations = rand(Poisson(mu),length(alive))

    celltype = zeros(Int64,2047)

    population = popstruct(alive,mutations,celltype,parentid)
end

function idealbranching(numgen,mu,clonalpeak)

    #an an ideal branching process we do not need much of the above functions. We can easily write down
    #what the parent_id vector would look like for a given number of generations given that every cell divides into 2.

    x=zeros(numgen)
    power=0:1:numgen-1

    for i = 1:numgen
        x[i]=2^power[i]
    end

    #number of alive is given by 2^numgen
    alive=zeros(Int64,round(Int64,sum(x)))
    #alive[end-x[end]+1:end]=1
    alive[round(Int64,length(alive)-x[length(x)]+1:length(alive))]=1

    parentid = zeros(Int64,round(Int64,sum(x)))
    parentid[2:round(Int64,sum(x))] = round(Int64,floor(1:0.5:(round(Int64,sum(x)))/2))

    #draw number of mutations acquired at each division from poisson distribution
    #modify the number of mutations a founder cell has acquired
    #Note mu is the number of mutation per cell division, for a diploid cell divide this by 2
    #to get the number of mutations per genome per cell division.
    mutations = rand(Poisson(mu),length(alive))
    mutations[1]=clonalpeak*mu

    celltype = zeros(Int64,round(Int64,sum(x)))

    population = popstruct(alive,mutations,celltype,parentid)

end

function idealsampling(numgen,mu,clonalpeak,normalcontamination,minrange,maxrange,det_limit,ploidy,read_depth)
    #sample from an "ideal branching population".

    #compute population struct, calculate allele_frequencies and convert array of cells
    # to array of alleles.
    pop=idealbranching(numgen,mu,clonalpeak)
    pop,numalive=allele_frequency(pop)
    pop=pop_expand(pop)

    #set ploidy and adjust population accordingly
    pop.alive=round(Int64,pop.alive/ploidy)
    p=pop.alive/sum(pop.alive);

    #calculate the percentage of alleles we should sample given a target "read depth"
    samp_percent = read_depth/numalive

    #find indices of alleles above detection limit in original unsampled population.
    above_detlimit=findin(pop.alive/numalive,filter((x)->x>det_limit,pop.alive/numalive))

    #sample alleles using multinomial sampling
    samp_alleles=rand(Multinomial(round(Int64,sum(pop.alive)*samp_percent),p))

    #Calculate distribution of values to represent read depth using Binomial sampling.
    depth = rand(Binomial(numalive,samp_percent),length(samp_alleles[above_detlimit]))

    #Add % of normal contamination
    depth = round(depth./(1-normalcontamination))

    #calculate VAF
    VAF=samp_alleles[above_detlimit]./depth

    #data for histogram
    x,y=hist(VAF,0.0:0.01:1)
    DFhist = DataFrame(VAF=x[1:end-1],freq=y)

    DFcumsum,lmfit,rsq = analysedata(VAF,maxrange,minrange)

    return DFhist,DFcumsum,coef(lmfit)*log(2),rsq
end

function fullmodel_sampling(popsize,mu,clonalpeak,normalcontamination,minrange,maxrange,probmat,lambda,det_limit,read_depth)

    pop,S=tumourgrow_nonzero(popsize,probmat,mu,500,lambda)
    pop,numalive=allele_frequency(pop)
    pop.mutations[1]=clonalpeak
    pop=pop_expand(pop)

    #set ploidy and adjust population accordingly
    ploidy=2
    pop.alive=round(Int64,pop.alive/ploidy)
    p=pop.alive/sum(pop.alive);
    samp_percent = read_depth/numalive

    above_detlimit=findin(pop.alive/numalive,filter((x)->x>det_limit,pop.alive/numalive))

    samp_alleles=rand(Multinomial(round(Int64,sum(pop.alive)*samp_percent),p))

    numcells = rand(Binomial(numalive,samp_percent),length(samp_alleles[above_detlimit]))
    numcells = round(numcells./(1-normalcontamination))

    VAF=samp_alleles[above_detlimit]./numcells

    #data for histogram
    x,y=hist(VAF,0.0:0.01:1)
    DFhist = DataFrame(VAF=x[1:end-1],freq=y)

    DFcumsum,lmfit,rsq = analysedata(VAF,maxrange,minrange)

    return DFhist,DFcumsum,coef(lmfit)*lambda,rsq
end

function fullmodel_sampling(popsize,mu,clonalpeak,normalcontamination,minrange,maxrange,probmat,lambda,det_limit,read_depth,selectiongen)
    #this is for the case where we have an offspring probability function that is not p=[0,0,1]

    pop,S=tumourgrow_nonzero(popsize,probmat,mu,selectiongen,lambda)
    pop,numalive=allele_frequency(pop)
    pop.mutations[1]=clonalpeak
    pop=pop_expand(pop)

    #set ploidy and adjust population accordingly
    ploidy=2
    pop.alive=round(Int64,pop.alive/ploidy)
    p=pop.alive/sum(pop.alive);
    samp_percent = read_depth/numalive
    above_detlimit=findin(pop.alive/numalive,filter((x)->x>det_limit,pop.alive/numalive))

    samp_alleles=rand(Multinomial(round(Int64,sum(pop.alive)*samp_percent),p))

    numcells = rand(Binomial(numalive,samp_percent),length(samp_alleles[above_detlimit]))
    numcells = round(numcells./(1-normalcontamination))


    VAF=samp_alleles[above_detlimit]./numcells

    #data for histogram
    x,y=hist(VAF,0.0:0.01:1)
    DFhist = DataFrame(VAF=x[1:end-1],freq=y)

    DFcumsum,lmfit,rsq = analysedata(VAF,maxrange,minrange)

    return DFhist,DFcumsum,coef(lmfit)*lambda,rsq
end


function multsample(numgen,num_samples,mut_rate,clon_peak,norm_cont,min_range,max_range,det_limit,ploidy,read_depth)

    #Run simulation num_samples times to get a distribution of paramater fits

    #initialize arrays where mutation rate and rsq values from fits will be stored
    mu=zeros(Float64,num_samples)
    rsq=zeros(Float64,num_samples)

    for i = 1:num_samples
        DF,DF1,x,r = idealsampling(11,mut_rate,clon_peak,norm_cont,min_range,max_range,det_limit,ploidy,read_depth);
        mu[i] = x[1]
        rsq[i] = r
    end

    statvec = [mean(mu),mean(rsq)]
    return statvec,mu,rsq
end

function tumourgrow_sel(maxpopsize::Int64,p,mu,selectiongen::Int64)
    #variation of tumourgrow function where the mutation rate mu can be changed

    # Tumour grows to maxpopsize given probability convolution matrix and probabilities. And mutation rate mu
    #
    # The simulation of tumour growth is specified by a vector which we call pop
    # and a vector called parent ID. alive designates which crypts are still "alive"
    # 1 being alive 0 being dead and parent ID designates the parents of each crypt.
    # The ID simply refers to the index of the vector. For example if we grow to a size of
    # 4 crypts we might get a output that looks as follows. pop = [0,0,0,1,1,1,1],
    # and parentid = [0,1,1,2,2,3,3]. This would correspond to an idealised branching
    # process whereby 1 crypts splits into 2 which splits into 4. Using FD_reconstruct we can
    # then calculate the allele frequency given that all mutations are passed onto daughter
    # crypts.

    selectionevent=0
    popsize_event=0


    #population begins with 1 crypt, parentID is 0, number of mutations is zero and celltype is 0:
    alive=[1]
    mutations=[0]
    parentid=[0]
    celltype=[0]

    population = popstruct(alive,mutations,celltype,parentid)

    #define a variable selectionevent that increases by 1 at each generation, then use this
    #to specofy the generation at which a mutant allele with fitness advantage enters the population.
    selectionevent = 1
    mutantid = 0

    while sum(population.alive)<maxpopsize
        if sum(population.alive) == 0
            break
        end

        #generate next generation given current generation and probabilities
        popnew = newGenerationGrow2( sum(population.alive) , population.celltype[end-sum(population.alive)+1:end] , p)

        #set parent ID of new population
        x = setID( append!(zeros(Int,length(population.alive) - sum(population.alive)) , popnew ) )

        # append new parentID onto previous parent IDs
        append!(population.parentid,x)

        #add new population onto old one. Remembering that this simply
        #corresponds to alive (1) or dead (0)
        population.alive = append!(zeros(Int,length(population.alive)) , ones(Int,sum(popnew)) )

        #set celltype of new population
        population.celltype = setcelltype(population.parentid,population.celltype)

        #increase selection event by 1
        selectionevent = selectionevent + 1

        #introduce fitter cell at correct generation
        if selectionevent == selectiongen
            population.celltype[end] = 1
            popsize_event = sum(population.alive)
            mutantid = population.parentid[end]
        end

    end
    # We now calculate the number of mutations that happen at each divison.
    # This can be done at the end.

    #use poisson to sample mutations for each cell
    population.mutations = rand(Poisson(mu),length(population.alive))

    simstats = simstat(sum(population.celltype),sum(population.alive),selectionevent,popsize_event,selectiongen,mutantid)

    return population , simstats
end


function idealsampling_selection(popsize,p,mu,selectiongen,clonalpeak,normalcontamination,minrange,maxrange,det_limit,ploidy,read_depth)
    #sample from an "ideal branching population".

    #compute population struct, calculate allele_frequencies and convert array of cells
    # to array of alleles.
    pop,S=tumourgrow_sel(popsize,p,mu,selectiongen)

    pop.mutations[1]=rand(Poisson(clonalpeak*mu))
    allele_frequency(pop)

    numalive=maximum(pop.alive)
    pop=pop_expand(pop)

    #set ploidy and adjust population accordingly
    pop.alive=round(Int64,pop.alive/ploidy)
    prob_vec=pop.alive/sum(pop.alive);

    #calculate the percentage of alleles we should sample given a target "read depth"
    samp_percent = read_depth/numalive

    #find indices of alleles above detection limit in original unsampled population.
    above_detlimit=findin(pop.alive/numalive,filter((x)->x>det_limit,pop.alive/numalive))

    #sample alleles using multinomial sampling
    samp_alleles=rand(Multinomial(round(Int64,(sum(pop.alive))*samp_percent),prob_vec))


    #Calculate distribution of values to represent read depth using Binomial sampling.
    #depth = rand(Binomial(numalive,samp_percent),length(samp_alleles[above_detlimit]))
    numtrials=numalive * length(pop.alive)
    depth=rand(Multinomial(round(Int64,samp_percent*numtrials),length(pop.alive)))

    #Add % of normal contamination
    depth = round(depth./(1-normalcontamination))

    VAF=samp_alleles[above_detlimit]./depth[above_detlimit]

    #data for histogram
    x,y=hist(VAF,0.0:0.01:1)
    DFhist = DataFrame(VAF=x[1:end-1],freq=y)

    DFcumsum,lmfit,rsq = analysedata(VAF,maxrange,minrange)


    return DFhist,DFcumsum,coef(lmfit),rsq,S
end

function tumourgrow_mutation(maxpopsize::Int64,p,mu_old,mu_new,selectiongen::Int64)
    #variation of tumourgrow function where the mutation rate mu can be changed

    # Tumour grows to maxpopsize given probability convolution matrix and probabilities. And mutation rate mu
    #
    # The simulation of tumour growth is specified by a vector which we call pop
    # and a vector called parent ID. alive designates which crypts are still "alive"
    # 1 being alive 0 being dead and parent ID designates the parents of each crypt.
    # The ID simply refers to the index of the vector. For example if we grow to a size of
    # 4 crypts we might get a output that looks as follows. pop = [0,0,0,1,1,1,1],
    # and parentid = [0,1,1,2,2,3,3]. This would correspond to an idealised branching
    # process whereby 1 crypts splits into 2 which splits into 4. Using FD_reconstruct we can
    # then calculate the allele frequency given that all mutations are passed onto daughter
    # crypts.

    selectionevent=0
    popsize_event=0


    #population begins with 1 crypt, parentID is 0, number of mutations is zero and celltype is 0:
    alive=[1]
    mutations=[0]
    parentid=[0]
    celltype=[0]

    population = popstruct(alive,mutations,celltype,parentid)

    #define a variable selectionevent that increases by 1 at each generation, then use this
    #to specofy the generation at which a mutant allele with fitness advantage enters the population.
    selectionevent = 1
    mutantid = 0

    while sum(population.alive)<maxpopsize
        if sum(population.alive) == 0
            break
        end

        #generate next generation given current generation and probabilities
        popnew = newGenerationGrow2( sum(population.alive) , population.celltype[end-sum(population.alive)+1:end] , p)

        #set parent ID of new population
        x = setID( append!(zeros(Int,length(population.alive) - sum(population.alive)) , popnew ) )

        # append new parentID onto previous parent IDs
        append!(population.parentid,x)

        #add new population onto old one. Remembering that this simply
        #corresponds to alive (1) or dead (0)
        population.alive = append!(zeros(Int,length(population.alive)) , ones(Int,sum(popnew)) )

        #set celltype of new population
        population.celltype = setcelltype(population.parentid,population.celltype)

        #increase selection event by 1
        selectionevent = selectionevent + 1

        #introduce fitter cell at correct generation
        if selectionevent == selectiongen
            population.celltype[end] = 1
            popsize_event = sum(population.alive)
            mutantid = population.parentid[end]
        end

    end
    # We now calculate the number of mutations that happen at each divison.
    # This can be done at the end.

    simstats = simstat(sum(population.celltype),sum(population.alive),selectionevent,popsize_event,selectiongen,mutantid)

    #use poisson to sample mutations for each cell
    population.mutations = zeros(Int64,length(population.alive))
    population.mutations[1:simstats.mutantid] = rand(Poisson(mu_old),simstats.mutantid)
    population.mutations[(simstats.mutantid+1):end] = rand(Poisson(mu_new),length(population.mutations)-simstats.mutantid)

    return population , simstats
end


function idealsampling_mutation(popsize,p,mu_old,mu_new,selectiongen,clonalpeak,normalcontamination,minrange,maxrange,det_limit,ploidy,read_depth)

    #sample from an "ideal branching population".

    #compute population struct, calculate allele_frequencies and convert array of cells
    # to array of alleles.
    pop,S=tumourgrow_mutation(popsize,p,mu_old,mu_new,selectiongen)

    pop.mutations[1]=rand(Poisson(clonalpeak*mu_old))
    allele_frequency(pop)

    numalive=maximum(pop.alive)
    pop=pop_expand(pop)

    #set ploidy and adjust population accordingly
    pop.alive=round(Int64,pop.alive/ploidy)
    prob_vec=pop.alive/sum(pop.alive);

    #calculate the percentage of alleles we should sample given a target "read depth"
    samp_percent = read_depth/numalive

    #find indices of alleles above detection limit in original unsampled population.
    above_detlimit=findin(pop.alive/numalive,filter((x)->x>det_limit,pop.alive/numalive))

    #sample alleles using multinomial sampling
    samp_alleles=rand(Multinomial(round(Int64,(sum(pop.alive))*samp_percent),prob_vec))


    #Calculate distribution of values to represent read depth using Binomial sampling.
    #depth = rand(Binomial(numalive,samp_percent),length(samp_alleles[above_detlimit]))
    numtrials=numalive * length(pop.alive)
    depth=rand(Multinomial(round(Int64,samp_percent*numtrials),length(pop.alive)))

    #Add % of normal contamination
    depth = round(depth./(1-normalcontamination))

    VAF=samp_alleles[above_detlimit]./depth[above_detlimit]

    #data for histogram
    x,y=hist(VAF,0.0:0.01:1)
    DFhist = DataFrame(VAF=x[1:end-1],freq=y)

    DFcumsum,lmfit,rsq = analysedata(VAF,maxrange,minrange)


    return DFhist,DFcumsum,coef(lmfit),rsq,S
end

function mixing(popsize,pct1,pct2,μ)

    popsize=1000
    numalive=popsize

    clonalpeak=2
    normalcontamination=0
    minrange=0.12
    maxrange=0.24
    det_limit=0.1
    ploidy=2
    read_depth=200

    clonesize = rand(3:7,3)
    clonesize[1] = 3
    clonesize[2] = 4
    clonesize[3] = 2


    founderclone = popsize*ones(Int64,clonesize[1])
    cloneA = round(Int64,popsize*pct1)*ones(Int64,clonesize[2])
    cloneB = round(Int64,popsize*pct2)*ones(Int64,clonesize[3])

    popmixed = popstruct([founderclone;cloneA;cloneB],rand(Poisson(μ),sum(clonesize)),[0],[0])
    unshift!(popmixed.alive,numalive)
    unshift!(popmixed.mutations,rand(Poisson(μ)))
    popmixed=pop_expand(popmixed)

    #set ploidy and adjust population accordingly
    popmixed.alive=round(Int64,popmixed.alive/ploidy)
    prob_vec=popmixed.alive/sum(popmixed.alive);

    #calculate the percentage of alleles we should sample given a target "read depth"
    samp_percent = read_depth/numalive

    #find indices of alleles above detection limit in original unsampled population.
    above_detlimit=findin(popmixed.alive/numalive,filter((x)->x>det_limit,popmixed.alive/numalive))

    #sample alleles using multinomial sampling
    samp_alleles=rand(Multinomial(round(Int64,(sum(popmixed.alive)*samp_percent)),prob_vec))


    #Calculate distribution of values to represent read depth using Binomial sampling.
    numtrials=numalive * length(popmixed.alive)
    depth=rand(Multinomial(round(Int64,samp_percent*numtrials),length(popmixed.alive)))

    #Add % of normal contamination
    depth = round(depth./(1-normalcontamination))

    VAF=samp_alleles[above_detlimit]./depth[above_detlimit]

    #data for histogram
    x,y=hist(VAF,0.0:0.01:1)
    DFhist = DataFrame(VAF=x[1:end-1],freq=y)

    #apply model fit
    DFcumsum,lmfit,rsq = analysedata(VAF,maxrange,minrange)


    return DFhist,DFcumsum,coef(lmfit),rsq,VAF

end

func(x)="1/$(round(1/(x+(1/0.24)),3))"

################################################################################
################################################################################

end
