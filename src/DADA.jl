
include("DADA\\errormodels.jl")
ρₚₒᵢₛ = StatsFuns.poispdf
get_pval(nⱼλᵢⱼ,aⱼ) = poisccdf(nⱼλᵢⱼ,aⱼ-1)

function get_norm(nⱼ,TAIL_APPROX_CUTOFF =1e-7)
    norm = (1.0 - exp(-nⱼ))
    if norm < TAIL_APPROX_CUTOFF
      norm = nⱼ - 0.5*nⱼ*nⱼ 
    end
    return norm
end

get_normed_pval(nⱼ,λᵢⱼ,aⱼ) = get_pval(nⱼ*λᵢⱼ,aⱼ)/get_norm(nⱼ*λᵢⱼ) #get_norm(nⱼ*λᵢⱼ)?



function compare!(seqs,quals,counts,ind,centres,partition,pvals,n,reads,compI,compV, err,mers,hdists,E_minmax; band_size = 16,kdist_cut = 0.42,gapless = false)
    affinegap = AffineGapScoreModel(match=5,
                                mismatch=-4,
                                 gap_open=-0,
                                 gap_extend=-8)
    base_vals = Base.ImmutableDict(DNA_A=>1,DNA_C=>2,DNA_G=>3,DNA_T=>4)
    l = keys(base_vals)
    mer_dist = 0.0
    λᵢⱼ = 0.0
    nⱼ = sum(counts[partition .== ind])
    N = sum(counts)
    #Threads.@threads for i in 1:n
     for i in 1:n
            hdist = 0.0
                        if !in(i,centres)          
                            mer_dist = kmer_distance(mers[:,ind],mers[:,i],length(seqs[ind]),length(seqs[i]))
                            if mer_dist < kdist_cut
                                if gapless & (kord_distance(seqs[ind], seqs[i]) == mer_dist)
                                    aln = Tuple([(seqs[ind][j],seqs[i][j]) for j in 1:length(seqs[i])])              
                                else
                                    if band_size >0
                                        pa = pairalign(OverlapAlignment(), seqs[ind], seqs[i],affinegap,true,band_size,band_size)  
                                    else   
                                        pa = pairalign(OverlapAlignment(), seqs[ind], seqs[i],affinegap)
                                    end            
                                aln = (collect(pa.aln))
                                end
                                qs = quals[i]
                                q = 0
                                λᵢⱼ = 1.0
                                for j in 1:length(aln)

                                    if in(aln[j][2],l) 
                                        q+=1
                                            j_ind = base_vals[aln[j][2]]
                                            i_ind = in(aln[j][1],l) ? base_vals[aln[j][1]] : base_vals[aln[j][2]]
                                            λᵢⱼ*=  err[(i_ind-1)*4+j_ind,qs[q]+1]
                                        hdist += 1
                                    end
                                end
                            else
                                λᵢⱼ = 0.0
                                hdist = Inf
                            end
                            if λᵢⱼ * N >= E_minmax[i]
                                push!(compI[i],length(centres))
                                push!(compV[i],λᵢⱼ)
                            end
                            if λᵢⱼ * counts[ind] > E_minmax[i]
                                E_minmax[i] = λᵢⱼ * counts[ind]
                            end
                            #λ[i] = λᵢⱼ
                            hdists[i] = hdist
                        end
                    end 
end

boolfunction(pvals,n,ωₐ) = any(pvals .* n .< ωₐ)

function minp(pvals,counts,centres,partition,n = length(pvals))
    minval = 1.0
    min_ind = 1
        @inbounds for i in 1:n
        if !in(i,centres)
            if pvals[i] < minval
                minval = pvals[i]
                min_ind = i
            elseif (pvals[i] == minval) & (counts[i]> counts[min_ind])
                minval = pvals[i]
                min_ind = i
            end
        end
    end
    return minval, min_ind
end

function shuffle(centres,partition,compI,compV,reads,pvals,counts,n,shuff,max_shuffle)
    change = false
    @inbounds for i in (1:n)
        if !in(i,centres) 
            λ,k = findmax(compV[i] .* reads[compI[i]])
            j = compI[i][k]
            prev_j = findfirst(centres .== [partition[i]])
            if j !== prev_j
                change = true
                reads[j] += counts[i]
                reads[prev_j] -= counts[i]
                partition[i] = centres[j]
            end      
            if shuff == max_shuffle
                pvals[i] = get_normed_pval(reads[j],compV[i][k],counts[i])
            end
        end
    end
    return change
end

function clustering(seqs, counts, quals,mers,err,ωₐ;band_size =16,log_p = false,log_λ=false, kdist_cut = 0.42,
    gapless = false, max_shuffle = 10)

    n = length(counts)
    pvals = zeros(n)
    centres = []
    E_minmax = zeros(n)
    hdists = zeros(n)
    compI = [Vector{Int}(undef,0) for i in 1:n]
    compV = [Vector{Float64}(undef,0) for i in 1:n]
    ind = findmax(counts)[2]
    centres = [ind]
    pvals[ind] = 1.0
    partition = fill(ind,n)
    compare!(seqs,quals,counts,ind,centres,partition,pvals,n,[],compI,compV, err,mers,hdists,E_minmax,band_size = band_size,kdist_cut = Inf,gapless = gapless)
    base_vals = Base.ImmutableDict(DNA_A=>1,DNA_C=>2,DNA_G=>3,DNA_T=>4)
    l = keys(base_vals)
    reads=[sum(counts)]

    for i in 1:n
        if !in(i,centres) 
            partition[i] = centres[1]
            pvals[i] = get_normed_pval(reads[1],compV[i][1],counts[i])
        end
    end
    coun = 0
    while true
        coun +=1
        println(coun)
        val, ind = minp(pvals,counts,centres,partition)
        if  boolfunction(val,n,ωₐ) & !in(ind, centres)
            partition[ind] = ind
            pvals[ind] = 1.0
            partition[centres] .= centres
            push!(centres,ind)
            hdists = hcat(hdists,zeros(n))               
            compare!(seqs,quals,counts,ind,centres,partition,pvals,n,reads,compI,compV, 
            err,mers,view(hdists,:,length(centres)),E_minmax,
            band_size = band_size,kdist_cut = kdist_cut,gapless = gapless)
        else
            break
        end
        reads = [sum(counts[partition .== i]) for i in centres]'
        change = true
        shuff = 1
        while shuff < max_shuffle 
            shuff +=1
            if change ==false
               shuff = max_shuffle
            end
            change = false
            change = shuffle(centres,partition,compI,compV,reads,pvals,counts,n,shuff,max_shuffle)         
        end
    end
    reads = [sum(counts[partition .== i]) for i in centres]'
    trans_mat = get_trans(centres,seqs,quals,counts,partition,n)
    cluster_map = Vector{Int}(undef,n)

    center_raws = []
    for i in 1:length(centres)
        boolmask = partition .== centres[i]
        raws = seqs[boolmask]
        abunds = counts[boolmask]
        inds = collect(1:n)[boolmask]
        push!(center_raws,raws[findmax(abunds)[2]])
        centres[i] = inds[findmax(abunds)[2]]
        partition[boolmask] .= centres[i] 
    end
    partition_map = Dict([centres[i] => i for i in 1:length(centres)]...)
    for i in 1:n
        cluster_map[i] = partition_map[partition[i]]
    end

    return (df =DataFrame([center_raws,vec(reads)],[:sequence,:abundance]),trans =  trans_mat,pvals = pvals,cluster_map =cluster_map)
end

function get_trans(centres,seqs,quals,counts,partition,n)
    affinegap = AffineGapScoreModel(match=5,
                                mismatch=-4,
                                 gap_open=-0,
                                 gap_extend=-8)
    base_vals = Base.ImmutableDict(DNA_A=>1,DNA_C=>2,DNA_G=>3,DNA_T=>4)
    l = keys(base_vals)
    trans_mat = zeros(Int,16,41)
    for i in 1:n
        if !in(i,centres)
            pa = pairalign(OverlapAlignment(), seqs[partition[i]], seqs[i],affinegap)
            qs = quals[i]#Int.(round.(derep.quality[i]))
            q = 0
            aln = (collect(pa.aln))
            for j in 1:length(aln)# need to make this start from start of seq[i] in alignment e.g. firt !== DNA_GAP
                    if in(aln[j][1],l) 
                        if in(aln[j][2],l)
                            q+=1        
                            i_ind = base_vals[aln[j][2]]
                            j_ind =  base_vals[aln[j][2]] 
                            trans_mat[(i_ind-1)*4+j_ind, qs[q]+1] +=counts[i]
                        end
                        
                    end
            end
        end
    end 
    return trans_mat
end
function learnErrors(dereps ::Vector{DataFrame}, max_tries = 10;nbases = 1e8, band_size = 16,ωₐ = 1e-40,kdist_cut = 0.42)

    st = time()
    base_count = 0
    dereps_for_errs = []
    mers = []
    for (j,derep) in enumerate(dereps)
        derep.sequence =LongDNASeq.(derep.sequence)
        push!(dereps_for_errs, derep)
        push!(mers,kmer_count(derep))
        base_count += sum(derep.count .* length.(derep.sequence))
        println("$(base_count) bases in $(j) samples")
        if base_count > nbases
            break
        end
    end
    #counts = derep.count
    err = ones(16,41)
    trans_mat = zeros(16,41)
    errs = [err]
    trans_mats =Vector{Array}(undef,length(dereps_for_errs))
    for i in 2:12
        #Threads.@threads for j in 1:length(dereps_for_errs)
            for j in 1:length(dereps_for_errs)
            derep = dereps_for_errs[j]
            trans_mats[j]= clustering(derep.sequence,derep.count,derep.quality,mers[j],err,ωₐ,
             band_size = band_size,log_p = true, kdist_cut = kdist_cut).trans
        end
        trans_mat = sum(trans_mats)
        err = loess_errors(trans_mat)
        push!(errs,err)
        println("Time elapsed: $(time() - st) seconds")
       println(UnicodePlots.heatmap(errs[i] .- errs[i-1], title = "Change from previous error estimate"))
       if any([errs[i][.!isnan.(errs[i])] ≈ errs[j][.!isnan.(errs[j])] for j in 1:i-1])
            println("Successfully converged")
            break
       end
    end
    return errs,trans_mat, trans_mats
end
function learnErrors(dereps ::Vector{String};nbases = 1e8, band_size = 16,ωₐ = 1e-40,kdist_cut = 0.42)
    dereps = dereplicate.(dereps)
    return learnErrors(dereps,nbases = nbases , band_size = band_size,ωₐ = ωₐ,kdist_cut = kdist_cut)
end

function learnErrors(derep ::DataFrame;nbases = 1e8, band_size = 16,ωₐ = 1e-40,kdist_cut = 0.42)
    return learnErrors([derep],nbases = nbases , band_size = band_size,ωₐ = ωₐ,kdist_cut = kdist_cut)
end


function dada(derep,err; band_size = 16,ωₐ = 1e-40,kdist_cut = 0.42,log_λ = false, log_p = true, gapless = false)
    println(derep.sample)
    counts = derep.count
    derep.df.sequence =LongDNASeq.(derep.sequence)
    seqs = LongDNASeq.(derep.sequence)
    qs = derep.quality
    partition = Vector{Int}(undef,length(counts))
    mers = kmer_count(derep.df)
    seq_abund = []
    pvals = []
    partition=[]
    λ= []
    hdists = []
    reads = []
    trans_mat = []
        clust = clustering(seqs,counts,qs,mers,err,ωₐ, 
        band_size = band_size,log_λ = log_λ, log_p = log_p, kdist_cut = kdist_cut, gapless = gapless)
  
    return clust
end

function mergePairs(dadaF,derepF,dadaR,derepR)

    rf = dadaF.cluster_map[derepF.read_map]
    rr = dadaR.cluster_map[derepR.read_map]
    pairsArr = hcat(rf,rr)
    unqs = unique(pairsArr, 1)

    Fseqs = dadaF.df.sequence[unqs[:,1]]
    Rseqs = dadaR.df.sequence[unqs[:,2]]
    affinegap = AffineGapScoreModel(match=5,
    mismatch=-4,
     gap_open=-64,
     gap_extend=-64)
     alignments = pairalign.(OverlapAlignment(), Fseqs,Rseqs,affinegap)
    counts = countmap(eachrow(unqs))
    abundance = [counts[unq] for unq in eachrow(unqs)]
    df = Dataframe(hcat(pairsArr,aln.(alignments),abundance),[:rev,:fwd,:aln,:abundance])
    return df
end


function pairalign(::OverlapAlignment, a::S1, b::S2, score::AffineGapScoreModel{T}, banded::Bool,lower_offset,upper_offset;
    score_only::Bool=false) where {S1,S2,T}
    m = length(a)
    n = length(b)
    if banded
        if m > n
            L = m - n + lower_offset
            U = upper_offset
        else
            L = lower_offset
            U = n - m + upper_offset
        end
        bnw = BioAlignments.BandedNeedlemanWunsch{T}(m, n, L, U)
        score = BioAlignments.run!(bnw, a, b, score.submat,
        T(0), T(0), score.gap_open, score.gap_extend, T(0), T(0),
        T(0), T(0), score.gap_open, score.gap_extend, T(0), T(0),
        )
        if score_only
            return BioAlignments.PairwiseAlignmentResult{S1,S2}(score, true)
        else
            a′ = BioAlignments.traceback(bnw, a, b, (m, n))
            return BioAlignments.PairwiseAlignmentResult(score, true, a′, b)
        end
    end
end





