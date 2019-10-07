@everywhere push!(LOAD_PATH,"./mdp")
@everywhere push!(LOAD_PATH, ".")
@everywhere push!(LOAD_PATH, pwd())
@everywhere using LocalFunctionApproximation
@everywhere using GridInterpolations
@everywhere using POMDPs
@everywhere using HDF5
@everywhere using POMDPModelTools
#@everywhere using HCAS
@everywhere include("../UAM_Speed/mdp/HCAS.jl")
@everywhere using SharedArrays

function get_Q(tau, init_qmat; factor=1.0)
	mdp = HCAS_MDP()
	mdp.currentTau = tau
    mdp.factor = factor

	ns = n_interpolating_points(interp)
	na = length(actions(mdp))

	state_chunks = split_states(ns, nprocs()-1)

    qmat_begin  = SharedArray{Float64}((ns, na), init = S -> S[localindices(S)] = init_qmat[localindices(S)])
    qmat_end  = SharedArray{Float64}((ns, na), init = S -> S[localindices(S)] = init_qmat[localindices(S)])

    interp_points = get_all_interpolating_points(interp)
    S = statetype(typeof(mdp))
    interp_states = Vector{S}(undef, ns)
    for (i,pt) in enumerate(interp_points)
        interp_states[i] = POMDPs.convert_s(S, pt, mdp)
    end
    shared_states = SharedArray{S}(ns, init = S -> S[localindices(S)] = interp_states[localindices(S)])
    pool = CachingPool(workers())

    results = pmap(x -> solve_chunk(mdp, shared_states, qmat_begin, qmat_end, x), pool, state_chunks)
    return qmat_end
end

@everywhere function solve_chunk(mdp::M, 
                    states::SharedArray{S, 1}, 
                    qmat_begin::SharedArray{Float64, 2},
                    qmat_end::SharedArray{Float64, 2},
                    state_indices::UnitRange{Int64}
                    ) where {M <: Union{MDP, POMDP}, S}

    discount_factor = discount(mdp)
    for istate=state_indices
        s = states[istate]
        sub_aspace = actions(mdp)

        iaction = 0
        for a in sub_aspace
            iaction += 1
            dist = transition(mdp, s, a) # creates distribution over neighbors
            u = 0.0
            for (sp, p) in weighted_iterator(dist)
                p == 0.0 ? continue : nothing # skip if zero prob
                r = reward(mdp, s, a)
                sp_point = POMDPs.convert_s(Vector{Float64}, sp, mdp)
                sps, probs = interpolants(interp.grid, sp_point)
                for (spi, probi) in zip(sps,probs)
                    u += probi * p * (r + discount_factor * maximum(qmat_begin[spi,:]))
                end
            end
            qmat_end[istate, iaction] = u
        end # action
    end # state loop
    return 
end

function split_states(ns::Int64, n_procs::Int64)
    state_chunks = Vector{UnitRange{Int64}}(undef, n_procs)
    stride = div(ns, n_procs)
    for i=0:n_procs-1
        si = i*stride + 1
        ei = (i + 1)*stride
        if i == n_procs-1
            ei = ns
        end
        state_chunks[i+1] = si:ei
    end 
    return state_chunks
end
