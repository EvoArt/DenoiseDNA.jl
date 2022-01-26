nanmean(x) = mean(filter(!isnan, x))
nanquantile(x,p) = quantile(filter(!isnan, x),p)

function quality_plot(file::String)
    qs = Vector{Vector{Float64}}(undef,0)
    reader = get_reader(file)
    for record in reader
        push!(qs,Vector{Float64}(quality(record)))
    end
    close(reader)
    q_array = float.(rstack(qs; fill=NaN))
    Iq_array = Int.(rstack(qs))
    max_len = maximum(length.(qs))
    q25,q50,q75 = Vector{Float64}(undef,max_len),Vector{Float64}(undef,max_len),Vector{Float64}(undef,max_len)
    Threads.@threads for i in 1:max_len
        #q25[i],q50[i],q75[i] = nanquantile(q_array[i,:],[0.25,0.5,0.75])
        q50[i] = nanmean(q_array[i,:])
    end
    rows,cols = size(q_array)
    hm_mat = zeros(max_len,41)
     @turbo for row in 1:rows
        for col in 1:cols
                hm_mat[row,Iq_array[row,col]+1] +=1
        end
    end
    fig = Figure()
    ax = Axis(fig[1,1])
    CairoMakie.heatmap!(ax,hm_mat[:,2:end], colormap = (:viridis,0.7))
    CairoMakie.lines!(ax,q50, color = :red)
    return fig
end

function quality_plot(file::String,ax)
    qs = Vector{Vector{Float64}}(undef,0)
    reader = get_reader(file)
    for record in reader
        push!(qs,Vector{Float64}(quality(record)))
    end
    close(reader)
    q_array = float.(rstack(qs; fill=NaN))
    Iq_array = Int.(rstack(qs))
    max_len = maximum(length.(qs))
    q25,q50,q75 = Vector{Float64}(undef,max_len),Vector{Float64}(undef,max_len),Vector{Float64}(undef,max_len)
    Threads.@threads for i in 1:max_len
        #q25[i],q50[i],q75[i] = nanquantile(q_array[i,:],[0.25,0.5,0.75])
        q50[i] = nanmean(q_array[i,:])
    end
    rows,cols = size(q_array)
    hm_mat = zeros(max_len,41)
     @turbo for row in 1:rows
        for col in 1:cols
                hm_mat[row,Iq_array[row,col]+1] +=1
        end
    end
    
    CairoMakie.heatmap!(ax,hm_mat[:,2:end], colormap = (:viridis,0.7))
    CairoMakie.lines!(ax,q50, color = :red)
end

function quality_plot(dir::String,identifier::String,out_path = "out.pdf";cols = 2, title_len = 28, title_from =false)
    
    file_list = []
    titles =[]
    for (root, dirs, files) in walkdir(dir)
        for (i,path) in enumerate(joinpath.(root, files))
            if contains(path,identifier)
            push!(file_list,path)
            push!(titles,files[i])
            end
        end
    end
    #n = length(file_list) 
    tmps = []
    n = cols*3
    k = length(file_list) ÷ n
    for j in 1:k
        fig = Figure()
        set = (j-1)*n#+1
        xs = repeat(1:cols,(n+cols)÷cols)
        ys = vcat([fill(i,cols) for i in 1:((n+cols)÷cols)]...)
        for i in 1:max(n,length(file_list)-(k*n))
            try
            if title_from == false
                start = max(1,length(titles[i+set])-title_len)
                ax = Axis(fig[ys[i],xs[i]],title =  titles[i+set][start:end])
            else
                start = title_from
                ax = Axis(fig[ys[i],xs[i]],title =  titles[i+set][start:start+title_len])
            end
            quality_plot(file_list[i+set],ax)
            catch e
                println("Error encountered with file $(file_list[i+set])")
                println(e)
            end
        end
        save("tmpqualityplt$(j).pdf",fig)
        push!(tmps,"tmpqualityplt$(j).pdf")
    end
    read(`$(Poppler_jll.pdfunite()) $(tmps) $(out_path)`, String)# join temp pdfs together
    rm.(tmps); # remove temp pds after uniting them
    println("plots saved at $(out_path)")

end
