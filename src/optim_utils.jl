using ForwardDiff, SparseArrays

function newton_optim(f, df, d2f, u0; verbose=false, tols=(1e-12,1e-8), niter=1000, decent=false)
    tol1,tol2 = tols
    iter = 0
    u = copy(u0)
    while true # Newton iterations
        iter += 1
        f0,df0 = f(u), df(u)
        if !decent
            d2f0 = d2f(u)
            du_newton = d2f0 \ df0
        end
        du_decent = df0

        if decent
            step_types = [ du_decent ]
        else
            step_types = [ du_newton, du_decent]
        end
        
        delta,f1 = 0.0,0.0
        success = false
        for du = step_types
            delta = 1.0
            f1 = f0
            while true # Line-search
                u1 = u - delta*du
                f1 = f(u1)
                
                if f1 < f0
                    u = u1
                    success = true
                    break
                end
                delta /= 2.0
                
                if delta < tol1
                    # Too small stepsize
                    break
                end
            end
            if success break end
        end
        if !success
            error("Too small stepsize in every search direction, exiting")
        end

        
        if verbose; println("  Iter $iter, d = $delta, fobj = $f1"); end
        if abs(f1 - f0) / abs(f0) < tol2
            if verbose; println("Small changes, terminating"); end
            break
        end
        if iter >= niter
            println("Too many iterations, terminating")
            break
        end
    end
    
    u,iter
end

function optimize_mesh!(m::HighOrderMesh{D}, fobj; newton_args...) where {D}
    np = size(m.x,1)
    nn,nel = size(m.el)
    
    bndnodes = boundary_nodes(m)
    intnodes = setdiff(1:np, bndnodes)
    dofmap = reshape(1:D*np, np, D)
    dpdu = sparse(dofmap[intnodes,:][:], 1:D*length(intnodes), 1, D*np, D*length(intnodes))

    p0 = copy(m.x)
    
    d2fobj = quad -> ForwardDiff.hessian(fobj, quad)
    
    function Fobj(p,el)
        F = 0
        for iel = 1:nel
            F += fobj(p[el[:,iel],:], iel)
        end
        F
    end
    
    function dFobj(p,el)
        dF = zeros(size(p))
        for iel = 1:nel
            dF[el[:,iel],:] .+= ForwardDiff.gradient(q->fobj(q,iel), p[el[:,iel],:])
        end
        dF
    end

    function d2Fobj(p,el)
        dofmap = reshape(1:D*np, np, D)
        ix = reshape(dofmap[el,:], (nn*D,1,nel))
        ix = repeat(ix, 1, nn*D, 1)
        iy = permutedims(ix, (2,1,3))
        
        data = zeros(nn*D,nn*D,nel)
        for iel = 1:nel
            data[:,:,iel] = ForwardDiff.hessian(q->fobj(q,iel), p[el[:,iel],:])
        end
        
        sparse(ix[:], iy[:], data[:])
    end

    u2p(u) = p0 + reshape(dpdu * u, :, D)
    f_optim(u) = Fobj(u2p(u), m.el)
    df_optim(u) = dpdu' * dFobj(u2p(u), m.el)[:]
    d2f_optim(u) = dpdu' * d2Fobj(u2p(u), m.el) * dpdu

    u = zeros(size(dpdu,2))
    u,niter = newton_optim(f_optim, df_optim, d2f_optim, u; verbose=true, newton_args...)
    println("  --> Finished optimizing: Iterations = $niter")
    m.x[:] = u2p(u)
    nothing
end
