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


include("Gillespie.jl")

using .Gillespie

using DifferentialEquations

export Model,Solution,DeterministicModel,StochasticModel,deterministic_step

abstract type Model end

struct DeterministicModel{P,T,U<:Number} <: Model
    ODEmodel::Function              # Deterministic version of the Model
    stochastic_threshold::Number    # When to switch to/from the StochasticModel
    progression_threshold::Number   # A factor that determines progression
    initial_population::Vector{T}
    params::P                       # Parameters to be fed to the ODEmodel
    time_domain::Tuple{U,U}         # Max length of simulation go. (typically (0, tf) where tf is final time)
    in_stochastic_region::Bool      # Modification to ODE if in this region
end

struct StochasticModel <: Model
    gp::GillespieProblem
end

mutable struct Solution{T<:Number,U}

    time_array::Vector{T}
    population_array::U
    patient_outcome::String
    stochastic_compartment::Vector{Int}

end

Solution(time,pop) = Solution([time],[pop],"Unknown",[3])

function deterministic_step(dm::Model)

    # Check to see if the tumor has entered the stochastic regime
    enter_stochastic_region(u,t,integrator) = dm.stochastic_threshold-u[3]
    affect!(integrator) = terminate!(integrator)
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