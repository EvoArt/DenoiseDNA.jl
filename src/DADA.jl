ρₚₒᵢₛ = poispdf
#get_pval(nⱼλᵢⱼ,aⱼ) = (1/(1-ρₚₒᵢₛ(nⱼλᵢⱼ,0)))*(1-poiscdf(nⱼλᵢⱼ,aⱼ))
get_pval(nⱼλᵢⱼ,aⱼ) = poisccdf(nⱼλᵢⱼ,aⱼ-1)
log_pval(nⱼλᵢⱼ,aⱼ)  = poislogccdf(nⱼλᵢⱼ,aⱼ-1)

function get_norm(nⱼ,TAIL_APPROX_CUTOFF =1e-7)
    norm = (1.0 - exp(-nⱼ))
    if norm < TAIL_APPROX_CUTOFF
      norm = nⱼ - 0.5*nⱼ*nⱼ 
    end
    return norm
end
function log_norm(nⱼ,TAIL_APPROX_CUTOFF =log(1e-7))
    norm = log(1.0 - exp(-nⱼ))
    if norm < TAIL_APPROX_CUTOFF
      norm = log(nⱼ - 0.5*nⱼ*nⱼ) 
    end
    return norm
end
get_normed_pval(nⱼ,λᵢⱼ,aⱼ) = get_pval(nⱼ*λᵢⱼ,aⱼ)/get_norm(nⱼ*λᵢⱼ) #get_norm(nⱼ*λᵢⱼ)?
log_normed_pval(nⱼ,λᵢⱼ,aⱼ) = log_pval(nⱼ*λᵢⱼ,aⱼ) -log_norm(nⱼ*λᵢⱼ) #get_norm(nⱼ*λᵢⱼ)?


function loess_errors(mat)
    est = Array{Float64}(undef,0,size(mat)[2])
    qs = float.(1:41)
      for j in 1:4
          group = 4*(j-1) 
          tot = vec(sum(mat[group+1:group+4,:], dims = 1))
          for i in 1:4
          if i != j
              errs  = mat[group+i,:]
              rlogp = log10.((errs .+1) ./tot)  # 1 psuedocount for each err, but if tot=0 will give NA
              rlogp .= inf2nan.(rlogp)
              model = loess(qs,rlogp)
              pred = Loess.predict(model, qs)
              R"df <- data.frame(q=$(qs), errs=$(errs), tot=$(tot), rlogp=$(rlogp))"
              R"mod.lo <- loess(rlogp ~ q, df, weights=tot)" ###!
          #        mod.lo <- loess(rlogp ~ q, df)
              pred = rcopy(R"predict(mod.lo, $(qs))")
              #println(ismissing.(pred)')
              pred[ismissing.(pred)] .= NaN
              maxrli = findlast(x-> !isnan(x),pred)
              minrli = findfirst(x-> !isnan(x),pred)
              pred[maxrli:end] .= pred[maxrli]
              pred[1:minrli] .= pred[minrli]
  
              est = vcat(est, 10 .^pred')
              end
          end
      end
      clamp!(est,1e-7,1.0) # editing out this line give closer output to R
      err = vcat(1 .-sum(est[1:3,:],dims = 1), est[1:3,:],
                 est[4,:]', 1 .-sum(est[4:6,:],dims = 1), est[5:6,:],
                 est[7:8,:], 1 .-sum(est[7:9,:],dims = 1), est[9,:]',
                 est[10:12,:], 1 .-sum(est[10:12,:], dims = 1))
      return err
  end
  
  
function log_compare!(seqs,quals,counts,ind,centres,partition,pvals,n,reads,λ, err,mers,hdists,E_minmax;band_size =16, kdist_cut = 0.42,gapless = false)
    affinegap = AffineGapScoreModel(match=5,
                                mismatch=-4,
                                 gap_open=-0,
                                 gap_extend=-8)
    base_vals = Base.ImmutableDict(DNA_A=>1,DNA_C=>2,DNA_G=>3,DNA_T=>4)
    l = keys(base_vals)
    mer_dist = 0.0
    nⱼ = sum(counts[partition .== ind])#counts[ind]#expected_reads_centre(seqs[ind],counts[ind],quals[ind],err)
    Threads.@threads for i in 1:n
            hdist = 0.0
        if !in(i,centres)          
            mer_dist = kmer_distance(mers[:,ind],mers[:,i],length(seqs[ind]),length(seqs[i]))
            if mer_dist < kdist_cut
                if gapless & (kord_distance(seqs[ind], seqs[i]) == mer_dist)
                    pa = [(seqs[ind][j],seqs[i][j]) for j in 1:length(seqs[i])]               
                elseif band_size >0
                    pa = pairalign(OverlapAlignment(), seqs[ind], seqs[i],affinegap,true,band_size,band_size)  
                else   
                    pa = pairalign(OverlapAlignment(), seqs[ind], seqs[i],affinegap)
                end      
                qs = quals[i]#Int.(round.(derep.quality[i]))
                q = 0
                aln = (collect(pa.aln))
                λᵢⱼ = 0.0

                for j in 1:length(aln)
   
                    if in(aln[j][2],l) 
                        q+=1
                        j_ind = base_vals[aln[j][2]]
                            i_ind = in(aln[j][1],l) ? base_vals[aln[j][1]] : base_vals[aln[j][2]]
                            λᵢⱼ+=  log(err[(i_ind-1)*4+j_ind,qs[q]+1])                    
                    end
                end
                
            else
                λᵢⱼ = -Inf
                hdist = Inf
            end
            if exp(log(counts[ind]) + λᵢⱼ) > E_minmax[i]
                E_minmax[i] = exp(log(counts[ind]) + λᵢⱼ)
                λ[i] = exp(λᵢⱼ)
            else λ[i] =-Inf
            end
            hdists[i] = hdist
        end
    end 
end

function compare!(seqs,quals,counts,ind,centres,partition,pvals,n,reads,λ, err,mers,hdists,E_minmax; band_size = 16,kdist_cut = 0.42,gapless = false)
    affinegap = AffineGapScoreModel(match=5,
                                mismatch=-4,
                                 gap_open=-0,
                                 gap_extend=-8)
    base_vals = Base.ImmutableDict(DNA_A=>1,DNA_C=>2,DNA_G=>3,DNA_T=>4)
    l = keys(base_vals)
    mer_dist = 0.0
    nⱼ = sum(counts[partition .== ind])#counts[ind]#expected_reads_centre(seqs[ind],counts[ind],quals[ind],err)
    Threads.@threads for i in 1:n
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
                qs = quals[i]#Int.(round.(derep.quality[i]))
                q = 0
                λᵢⱼ = 1.0
                    for j in 1:length(aln)#end_aln
                        #if start_aln <= j <= end_aln
                            #hdist += aln[j][1] !== aln[j][2]
                        #end
                        if in(aln[j][2],l) 
                            q+=1
                            #if in(aln[j][1],l) 
                                j_ind = base_vals[aln[j][2]]
                                i_ind = in(aln[j][1],l) ? base_vals[aln[j][1]] : base_vals[aln[j][2]]
                                λᵢⱼ*=  err[(i_ind-1)*4+j_ind,qs[q]+1]
                            hdist += 1

                            #end
                        end
                    end
                #end               
            else
                λᵢⱼ = 0.0
                hdist = Inf
            end
            if n * λᵢⱼ > E_minmax[i]
                E_minmax[i] = counts[ind] * λᵢⱼ
                λ[i] = λᵢⱼ
            else λ[i] = 0.0
            end
        
            #λ[i] = λᵢⱼ
            hdists[i] = hdist
        end
    end 
end

log_boolfunction(pvals,n,ωₐ) = any(pvals .+ log(n) .< log(ωₐ))

boolfunction(pvals,n,ωₐ) = any(pvals .* n .< ωₐ)

function minp(pvals,counts,centres,n = length(pvals))
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



function clustering(seqs, counts, quals,mers,err,ωₐ;band_size =16,log_p = false,log_λ=false, kdist_cut = 0.42,gapless = false)
    if log_p
        boolfunc = log_boolfunction
    else
        boolfunc = boolfunction
    end
    if log_λ
        compfunc! = log_compare!
    else
        compfunc! = compare!
    end

    n = length(counts)
    N = sum(counts)
    pvals = zeros(n)
    centres = []
    clusters = []
    expected_reads = zeros(n)
    E_minmax = zeros(n)
    hdists = zeros(n)
    λ = zeros(n)
    # first centre is most abundant sequence 
    ind = findmax(counts)[2]
    centres = [ind]
    pvals[ind] = log_p ? 0.0 : 1.0
    partition = fill(ind,n)
    compfunc!(seqs,quals,counts,ind,centres,partition,pvals,n,[],λ, err,mers,hdists,E_minmax,band_size = band_size,kdist_cut = Inf,gapless = gapless)
    base_vals = Base.ImmutableDict(DNA_A=>1,DNA_C=>2,DNA_G=>3,DNA_T=>4)
    l = keys(base_vals)
    reads=[sum(counts)]
    expected_reads .= λ .* reads
    E_minmax .= λ .* counts[ind]
    if log_p
        for i in 1:n
            if !in(i,centres) # need to actually do all the alignments as wel!!!!
                partition[i] = centres[1]
                pvals[i] = log_normed_pval(reads[1],λ[i],counts[i])
                #pvals[i] = isnan(pvals[i]) ? -Inf : pvals[i]         
                pvals[i] = counts[i] ==1 ? 0.0 : pvals[i]
            end
        end
    else
        for i in 1:n
            if !in(i,centres) 
                partition[i] = centres[1]
                pvals[i] = get_normed_pval(reads[1],λ[i],counts[i])
                #pvals[i] = isnan(pvals[i]) ? 0.0 : pvals[i]         
                pvals[i] = counts[i] ==1 ? 1.0 : pvals[i]
            end
        end
    end

    coun = 0
    while true#(boolfunc(pvals,n,ωₐ)) #& (length(centres) < 356)
        coun +=1
        #ind = findmin(pvals)[2]
        val, ind = minp(pvals,counts,centres)
        if  boolfunc(val,n,ωₐ) & !in(ind, centres)
        #if  boolfunc(val,N,ωₐ) & !in(ind, centres)
            partition[ind] = ind
            pvals[ind] = log_p ? 0.0 : 1.0
            partition[centres] .= centres
            push!(centres,ind)
            expected_reads = hcat(expected_reads,zeros(n))

            λ = hcat(λ,zeros(n))
            hdists = hcat(hdists,zeros(n))               
            compfunc!(seqs,quals,counts,ind,centres,partition,pvals,n,reads,view(λ,:,length(centres)), 
            err,mers,view(hdists,:,length(centres)),E_minmax,
            band_size = band_size,kdist_cut = kdist_cut,gapless = gapless)
        else
            break
        end
        reads = [sum(counts[partition .== i]) for i in centres]'
       # println([length(centres) sum(pvals .+ log(n) .< log(ωₐ)) sum(reads .==1)])
        change = true
        shuff = 1
        while shuff <10   
            shuff +=1
            if change ==false
               shuff = 10
            end
            change = false

            expected_reads .=  λ .* reads
            for bi in 1:length(centres)
                @inbounds for i in (1:n)[partition .== centres[bi]]
                    if !in(i,centres) 
                        j = findmax(expected_reads[i,:])[2]# slow!
                        prev_j = findfirst(centres .== [partition[i]])
                        if j !== prev_j
                            change = true
                            reads[j] += counts[i]
                            reads[prev_j] -= counts[i]
                            partition[i] = centres[j]
                            #expected_reads[:,j] .=  λ[:,j] .* reads[j]
                            #expected_reads[:,prev_j] .=  λ[:,prev_j] .* reads[prev_j]

                        end
                        if log_p & (shuff ==10)
                            pvals[i] = log_normed_pval(reads[j],λ[i,j],counts[i])
                            pvals[i] = counts[i] ==1 ? 0.0 : pvals[i]
                        elseif shuff ==10
                            pvals[i] = get_normed_pval(reads[j],λ[i,j],counts[i])
                            pvals[i] = counts[i] ==1 ? 1.0 : pvals[i]
                        end         
                    end
                end  
            end        
        end
        if coun /30 == coun ÷ 30
            println([length(centres) sum(pvals .+ log(n) .< log(ωₐ)) sum(reads .==1)])
            println(reads)
        end
       # println([sum(pvals .== -0.0) sum(0.0 .< pvals .< 1.0) sum(pvals .== 1.0)] )
    end
    reads = [sum(counts[partition .== i]) for i in centres]'
    trans_mat = get_trans(centres,seqs,quals,counts,partition,n)
    cluster_map = Vector{Int}(undef,n)
    partition_map = Dict([centres[i] => i for i in 1:length(centres)]...)
    for i in 1:n
        cluster_map[i] = partition_map[partition[i]]
    end
    inds = 1:n
    cs = []
    for i in 1:length(centres)
        boolmask = partition .== centres[i]
        raws = seqs[boolmask]
        abunds = counts[boolmask]
        push!(cs,raws[findmax(abunds)[2]])
    end

    #return (df =DataFrame([seqs[centres],vec(reads)],[:sequence,:abundance]),trans =  trans_mat,pvals = pvals,cluster_map =cluster_map)
    return (df =DataFrame([cs,vec(reads)],[:sequence,:abundance]),trans =  trans_mat,pvals = pvals,cluster_map =cluster_map)
end

function get_trans(centres,seqs,quals,counts,partition,n)
    affinegap = AffineGapScoreModel(match=5,
                                mismatch=-4,
                                 gap_open=-0,
                                 gap_extend=-8)
    base_vals = Base.ImmutableDict(DNA_A=>1,DNA_C=>2,DNA_G=>3,DNA_T=>4)
    l = keys(base_vals)
    trans_mat = zeros(Int,16,41)
    Threads.@threads for i in 1:n
        if !in(i,centres)
            pa = pairalign(OverlapAlignment(), seqs[partition[i]], seqs[i],affinegap)
            qs = quals[i]#Int.(round.(derep.quality[i]))
            q = 0
            aln = (collect(pa.aln))
            for j in 1:length(aln)# need to make this start from start of seq[i] in alignment e.g. firt !== DNA_GAP
                #if aln[j][1] !== aln[j][2]
                    if in(aln[j][2],l) 
                        q+=1
                        
                            j_ind = base_vals[aln[j][2]]
                            i_ind = in(aln[j][1],l) ? base_vals[aln[j][1]] : base_vals[aln[j][2]]
                            trans_mat[(i_ind-1)*4+j_ind, qs[q]+1] +=counts[i]
                        
                    end
                #end
            end
        end
    end 
    return trans_mat
end

function learnErrors(dereps ::Vector{String};nbases = 1e8, band_size = 16,ωₐ = 1e-40,kdist_cut = 0.42)
    dereps = dereplicate.(dereps)
    return learnErrors(dereps,nbases = nbases , band_size = band_size,ωₐ = ωₐ,kdist_cut = kdist_cut)
end

function learnErrors(derep ::DataFrame;nbases = 1e8, band_size = 16,ωₐ = 1e-40,kdist_cut = 0.42)
    return learnErrors([derep],nbases = nbases , band_size = band_size,ωₐ = ωₐ,kdist_cut = kdist_cut)
end

function learnErrors(dereps ::Vector{DataFrame};nbases = 1e8, band_size = 16,ωₐ = 1e-40,kdist_cut = 0.42)

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
    counts = derep.count
    err = ones(16,41)
    trans_mat = zeros(16,41)
    errs = [err]
    partition = Vector{Int}(undef,length(counts))
    partition=[]
    trans_mats =Vector{Array}(undef,length(dereps_for_errs))
    for i in 2:12
        Threads.@threads for j in 1:length(dereps_for_errs)
            derep = dereps_for_errs[j]
            _,trans_mats[j],_,_,_,_,_ = clustering(derep.sequence,derep.count,derep.quality,mers[j],err,ωₐ, band_size = band_size,log_p = true, kdist_cut = kdist_cut)
        end
        trans_mat = sum(trans_mats)
        err = loess_errors(trans_mat)
        push!(errs,err)
       # println(errs[i] .- errs[i-1])
        println(mean(abs.(errs[i] .- errs[i-1])))
        println(any([errs[i][.!isnan.(errs[i])] ≈ errs[j][.!isnan.(errs[j])] for j in 1:i-1]))
       println(UnicodePlots.heatmap(errs[i] .- errs[i-1]))
    end
    return errs,trans_mat, trans_mats
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
