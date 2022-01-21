inf2nan(x) = x == Inf ? NaN : x

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
      clamp!(est,1e-7,0.25) # editing out this line give closer output to R
      err = vcat(1 .-sum(est[1:3,:],dims = 1),
                 est[1:3,:],
                 est[4,:]', 
                 1 .-sum(est[4:6,:],dims = 1), 
                 est[5:6,:],
                 est[7:8,:], 
                 1 .-sum(est[7:9,:],dims = 1), 
                 est[9,:]',
                 est[10:12,:], 
                 1 .-sum(est[10:12,:], dims = 1))
      return err
  end

