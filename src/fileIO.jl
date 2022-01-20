# stats function

# get mean quality scores
function Q_mean(qs,RC = size(qs))
    mean_qs = zeros(Float64,RC[1])
   # mean_qs = Vector{UInt8}(undef,RC[1])
    @turbo for col in 1:RC[2]
        #mean_qs[col] = 0
        for row in 1:RC[1]
            mean_qs[row] += qs[row,col]
        end
    end
    return round.(mean_qs ./ RC[2])
end

# calculate expected errors
function EE(q,n = length(q))    
    ee = 0.0
    @turbo for i in 1:n
        ee+=10.0^(-Int(q[i])/10.0)
    end
    return ee
end

# readers 
function get_reader(file)
    if endswith(file,"gz")
        if contains(file,"fastq")
            reader = FASTQ.Reader(GzipDecompressorStream(open(file)))
        else
            reader = FASTA.Reader(GzipDecompressorStream(open(file)))
        end
    else
        if contains(file,"fastq")
            reader = open(FASTQ.Reader, file)
        else
            reader = open(FASTA.Reader, file)
        end
    end
    return reader
end

function get_files(dir,identifier; just_names = false)
    file_list = []
    sample_names =[]
    for (root, dirs, files) in walkdir(dir)
        for (i,path) in enumerate(joinpath.(root, files))
            if contains(path,identifier)
            push!(file_list,path)
            push!(sample_names,files[i])
            end
        end
    end
    return just_names ? sample_names : file_list
end

# dereplicate

function dereplicate(file; sample_name = "sample")
   # println(sample_name)
    try
    seqs = Vector{String}(undef,0) #LongDNASeq
    
    qs = Vector{Vector{UInt8}}(undef,0)
    reader = get_reader(file)
    for record in reader
        push!(seqs,sequence(record))
        push!(qs,quality(record))
    end
    close(reader)
    sp = sortperm(seqs)
    srt = seqs[sp]#sort(seqs) # map from rank to sequence/quality/id
    unqs = unique(srt)
    qsrt=qs[sp]
    counts = countmap(srt) # map from rank to abundance (much faster than table)
    mnq = Vector{Vector{Int}}(undef,length(unqs))
    abunds= Vector{Int}(undef,length(unqs))
    sample = fill(sample_name,length(unqs))
    row = 1
    for (i,seq) in enumerate(unqs)
        abund = counts[seq]
        #mnq[i] = vec(mean(hcat(qsrt[row:row+abund-1]),dims =2 ))[1]
        mnq[i] = Q_mean(hcat(qsrt[row:row+abund-1]...))
        abunds[i] = abund
        row += abund
    end
    df =DataFrame([sample,unqs,abunds,mnq],[:sample,:sequence,:count,:quality])
    sort!(df,:count,rev = true)
    unq_map = Dict([df.sequence[i] => i for i in 1:length(unqs)])
    read_map = Vector{Int}(undef,length(seqs))
    for i in 1:length(seqs)
        read_map[i] = unq_map[seqs[i]]
    end

    return (sample = sample_name,df =df, read_map = read_map, sequence = df.sequence,count = df.count,quality = df.quality)
    catch e
        println("Error encountered with file $(sample_name)")
        println(e)
    end
end


function dereplicate(dir, identifier)
    file_list = []
    sample_names =[]
    for (root, dirs, files) in walkdir(dir)
        for (i,path) in enumerate(joinpath.(root, files))
            if contains(path,identifier)
            push!(file_list,path)
            push!(sample_names,files[i])
            end
        end
    end
    n = length(file_list)
    dfs = []#Vector{DataFrame}(undef,0)
    for i in 1:n #  not sure why multi-threading slows this down 

        try
            push!(dfs,dereplicate(file_list[i], sample_name = sample_names[i]))
        catch e
            println("Error encountered with file $(file_list[i])")
            println(e)
        end
    end
return vcat(dfs...)
end


#trim and filter

function trim(file,out_path ="my-out.fastq";trunc_quality=0, trunc_length=Inf,max_length=Inf,min_length=0,trim_left=0,
    trim_right=0,max_N =0,min_quality = 0, max_expected_errors = Inf)

    trim_left+=1
    reader = get_reader(file)
    println(out_path)
    writer = open(FASTQ.Writer, out_path)
    for record in reader    
        seq = sequence(record)
        q = quality(record)
        if length(seq) <=max_length                 #max length
            seq = seq[trim_left:(end-trim_right)]     # trim left/right
            q = q[trim_left:(end-trim_right)]
            if trunc_quality > 0                    # trunc quality
                tq = findfirst(x-> x <trunc_quality,q)
                if tq !== nothing
                    seq = seq[1:tq]
                    q = q[1:tq]
                end
            end
            if trunc_length <= length(seq)          # trunc length
                seq = seq[1:trunc_length]
                q = q[1:trunc_length]
                if min_length <= length(seq)        #min length
                    if count(isambiguous, seq) <= max_N # max n
                        if minimum(q) >= min_quality # min quality
                            if EE(q) <= max_expected_errors # max expected errors 
                                write(writer,FASTQ.Record("id", seq, q))
                            end
                        end
                    end
                end
            end
        end
    end
    close(writer)
    close(reader)
end

function trim(dir,identifier,out_dir;trunc_quality=0, trunc_length=Inf,max_length=Inf,min_length=0,trim_left=0,
    trim_right=0,max_N =0,min_quality = 0, max_expected_errors = Inf)
    file_list = []
    sample_names =[]
    for (root, dirs, files) in walkdir(dir)
        for (i,path) in enumerate(joinpath.(root, files))
            if contains(path,identifier)
            push!(file_list,path)
            push!(sample_names,files[i])
            end
        end
    end
    n = length(file_list)
    for i in 1:n
        
        try
        file = file_list[i]
        out_path = joinpath(splitpath(out_dir)..., sample_names[i][1:end-3])
        trim(file,out_path,trunc_quality=trunc_quality, trunc_length=trunc_length,
        max_length=max_length,min_length=min_length,trim_left=trim_left,
        trim_right=trim_right,max_N =max_N,min_quality = min_quality, max_expected_errors = max_expected_errors)
        catch e
            println("Error encountered with file $(file_list[i])")
            println(e)
        end
    end
end
