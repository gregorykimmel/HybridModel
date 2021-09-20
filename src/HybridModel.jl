#=

This will create the structure HybridModel and simulate

    HybridModel will require:
        Stochastic version      (using Gillespie)
        Deterministic version   (using DifferentialEquations)
        HybridModel params (e.g. StochasticThreshold, finaltime, important events,etc.)


=#

#=

    include(Gillespie)
    include(Hybridmodel)

=#

module HybridModel

using DifferentialEquations

export Model,Solution,DeterministicModel,StochasticModel,deterministic_step

abstract type Model end
abstract type AbstractReactionNetwork end

struct ReactionMatrix{P} <: AbstractReactionNetwork
    sto_mat::Matrix{Int}        # How does population change?
    compute_rates::Function     # output should be same length as sto_mat
    rate_params::P              # To be fed into compute_rates function
end

struct ReactionNetwork{P} <: AbstractReactionNetwork
    sto_mat::Vector{Vector{Int}}    # How does population change?
    compute_rates::Function         # output should be same length as sto_mat
    rate_params::P                  # To be fed into compute_rates function
end

struct GillespieProblem{T<:Number,R<:AbstractReactionNetwork}
    initial_population::Vector{Int}
    time_domain::Tuple{T,T}
    rn::R
end

struct StochasticModel{S<:Number,T<:Number,R<:AbstractReactionNetwork} <: Model
    stochastic_population::Vector{Int}      # Population that is considered stochastic
    deterministic_population::Vector{S}     # Population that is large enough to be deterministic
    time_domain::Tuple{T,T}                 # Time to simulate
    rn::R                                   # Reaction network that governs population evolution
end

struct DeterministicModel{P,T,U<:Number} <: Model
    ODEmodel::Function              # Deterministic version of the Model
    stochastic_threshold::Number    # When to switch to/from the StochasticModel
    progression_threshold::Number   # A factor that determines progression
    initial_population::Vector{T}
    params::P                       # Parameters to be fed to the ODEmodel
    time_domain::Tuple{U,U}         # Max length of simulation go. (typically (0, tf) where tf is final time)
    in_stochastic_region::Bool      # Modification to ODE if in this region
end

mutable struct Solution{T<:Number,U}
    time_array::Vector{T}
    population_array::U
    patient_outcome::String
    stochastic_compartment::Vector{Int}
end

Solution(time,pop) = Solution([time],[pop],"Unknown",[3])

"find_event(rates) computes a tuple of the population rates for simulate."
function find_event(rates)

    # Get total to compute time to next event
    rate_total = sum(rates)

    # Compute tau (time till next event)
    τ = -log(rand())/rate_total

    # Normalize rates to determine which event occurs
    cum_rate = cumsum(rates./rate_total)

    # Roll the dice
    roll = rand()

    # find out which individual process occurred
    which_event = findfirst(x-> x >= roll,cum_rate)

    return τ, which_event

end

"Get initial time from GillespieProblem."
t₀(p::GillespieProblem) = p.time_domain[begin]

"Get final time from GillespieProblem."
t₁(p::GillespieProblem) = p.time_domain[end]

function simulate_single_event(p::GillespieProblem)

    rates = p.rn.compute_rates(t₀(p),p.initial_population,p.rn.rate_params)

    if sum(rates) == 0.0
        @warn "zero propagation!"
        return Inf,nothing
    end

    τ, which_event = find_event(rates)

    return τ, which_event

end

function deterministic_step(dm::Model)

    dm.in_stochastic_region = true

    # Check to see if the tumor has entered the stochastic regime
    enter_stochastic_region(u,t,integrator) = dm.stochastic_threshold-u[3]
    affect!(integrator) = begin 
        
        terminate!(integrator)
    end
    cb1 = ContinuousCallback(enter_stochastic_region,affect!,nothing,abstol=1e-9,save_positions=(false,false))

    # Check to see when the tumor exceeds some c*original size. We define 
    # this as progression and assume the third compartment is the tumor compartment...
    # Potentially we can use a component array instead so we don't need to track the compartment
    progression(u,t,integrator) = u[3] - dm.progression_threshold*dm.initial_population[3]
    affect2!(integrator) = terminate!(integrator)
    cb2 = ContinuousCallback(progression,affect2!,save_positions=(false,false))
    cbset = CallbackSet(cb1,cb2)

    prob = ODEProblem(dm.ODEmodel,dm.initial_population,dm.time_domain,dm.params,callback=cbset)

    sol = solve(prob,Tsit5(),callback=cbset,abstol=1e-9,reltol=1e-6)

    return sol

end

function stochastic_step(sm::Model)

    Gillespie.simulate_single_event(sm.gp)

end

end