# State transition function
function POMDPs.transition(mdp::HCAS_MDP, s::stateType, ra::actType)
    r = s[1]; t = s[2]; p = s[3]; vown = s[4]; vint = s[5]; pra = s[6];
    
    # Computation is faster when using vector of static size
    nextStates = MVector{9, stateType}(undef)
    nextProbs = @MVector(zeros(9))
    next_pra = ra
    ind=1

    # Compute probabilities of next states using sigma point sampling
    ownProbs, ownAccels = mdp.turns[pra]
    intProbs, intAccels = mdp.turns[-1]
    for i = 1:3
        for j = 1:3
            next_r,next_t,next_p,next_vown,nextAccels = dynamics(r,t,p,vown,vint,ownAccels[i],intAccels[j],pra)
            nextStates[ind] = (next_r,next_t,next_p,next_vown,next_vint,next_pra)
            nextProbs[ind]  = ownProbs[i]*intProbs[j]
            ind+=1
        end
    end

    return SparseCat(nextStates,nextProbs)
end

# Dynamic equations
function dynamics(r::Float64,t::Float64,p::Float64,vown::Float64,vint::Float64,ownAccel::Float64, intAccel::Float64, ra::Int)
    vown_new = vown + ownAccel
    vint_new = vint + intAccel

    x₁, y₁ = r*cos(t), r*sin(t)

    x₀_new = vown + 0.5ownAccel
    y₀_new = 0.0
    x₁_new = x₁ + vint*cos(p) + 0.5intAccel*cos(p)
    y₁_new = y₁ + vint*sin(p) + 0.5intAccel*sin(p)

    r_new = √((x₁_new - x₀_new)^2 + (y₁_new - y₀_new)^2))
    t_new = atan(y₁_new - y₀_new, x₁_new - x₀_new)
    
    return r_new, t_new, p, vown_new, vint_new
end