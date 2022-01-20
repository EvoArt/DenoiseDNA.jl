
base_vals = Base.ImmutableDict(DNA_A=>0,DNA_C=>1,DNA_G=>2,DNA_T=>3)

function get_mer_idx(mer,k=5)
    idx = 0
    @turbo for i in 1:k
        idx = 4*idx + base_vals[mer[i]]
    end
    return idx +1
end

function count_mers_4_dist(seq,k = 5)
    counts = zeros(Int,4^k)
    for mer in each(DNAMer{k}, seq)
        current_mer = mer.fw
        current_mer_idx = get_mer_idx(current_mer)
        counts[current_mer_idx] += 1
    end
    return counts
end
function count_mers_4_dist!(seq,mer_array,i,k = 5)
    for mer in each(DNAMer{k}, seq)
        current_mer = mer.fw
        current_mer_idx = get_mer_idx(current_mer)
        mer_array[current_mer_idx,i] += 1
    end
end
function kord_distance(n₁,n₂,k = 5) 
    l1,l2 = length(n₁), length(n₂)
    if l1 == l2
        dist = 0.0
        step = k-1
        for i in 1:l1-step
            dist += n₁[i:i+step] == n₂[i:i+step]
        end
        return 1-dist/(l1-step)
    else
        return -1.0
    end
end

kmer_distance(n₁,n₂,L1,L2,k=5) = 1-sum(min.(n₁, n₂))/(min(L1,L2) - k + 1)
function kmer_distance(derep,k = 5)
    seqs = derep.sequence
    n = length(seqs)
    mer_array = zeros(Int,4^k,n)
    for i in 1:n
        count_mers_4_dist!(seqs[i],mer_array,i)
    end
    dist_mat = Array{Float64}(undef,n,n)
    dist_mat = zeros(n,n)
    for i in 1:n
        for j in 1:i
            dist_mat[i,j] = kmer_distance(mer_array[:,i],mer_array[:,j],length(seqs[i]),length(seqs[j]),k)
        end
    end
    return dist_mat
end
function kmer_count(derep,k = 5)
    seqs = derep.sequence
    n = length(seqs)
    mer_array = zeros(Int,4^k,n)
    for i in 1:n
        count_mers_4_dist!(seqs[i],mer_array,i)
    end
    
    return mer_array
end
