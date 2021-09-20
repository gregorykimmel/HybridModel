using DifferentialEquations, Plots

"""
    plot(sol;stride,linecolor,label)
    
# Arguments
- `sol::Solution`:  The structure returned by simulate.
- `stride::Int`:    The number of steps to skip when plotting.
- `linecolor`:      Specifies the color of each population (e.g. [:red :blue]).
- `label`:          Specifies the name of each population (e.g. ["red" "blue"]).

"""
function Plots.plot(sol::Solution;stride::Int=1,linecolor=:black,label=nothing, yaxis=:normal)

    plot(
        sol.time_array[1:stride:end],
        hcat(sol.population_array[1:stride:end]...)',
        linecolor=linecolor,
        label=label,
        yaxis=yaxis
        )

end

abstract type Model end
abstract type AbstractReactionNetwork end

struct Options
    saveat
    abstol::Float64
    reltol::Float64
end

Options(;saveat=[],abstol=1e-6,reltol=1e-3) = Options(saveat,abstol,reltol)

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

struct StochasticModel{T<:Number,R<:AbstractReactionNetwork}
    current_time::T
    current_population::Vector{T}
    rn::R
end

struct DeterministicModel{P,T,U<:Number} <: Model
    ODEmodel::Function              # Deterministic version of the Model
    initial_population::Vector{T}
    params::P                       # Parameters to be fed to the ODEmodel
    stochastic_threshold::Number    # When to switch to/from the StochasticModel
    progression_threshold::Number   # A factor that determines progression
    time_domain::Tuple{U,U}         # Max length of simulation. (typically (0, tf) where tf is final time)
    options::Options
end

DeterministicModel(ode,init_pop,pars,sto_thres,prog_thres,time_domain) = begin
    DeterministicModel( ode,
                        init_pop,
                        pars,
                        sto_thres,
                        prog_thres,
                        time_domain,
                        Options()
                        )
end

mutable struct Solution{T<:Number,U}
    time_array::Vector{T}
    population_array::U
    patient_outcome::String
    in_stochastic_region::Bool      # Modification to ODE if in this region
    stochastic_compartment::Vector{Int}
end

Solution(time::Number,pop) = Solution([time],[pop],"Unknown",false,Vector{Int}())
Solution(time::Vector{<:Number},pop,is_stochastic::Bool) = Solution(time,pop,"Unknown",is_stochastic,[3])
Solution(ode,is_stochastic::Bool) = Solution(ode.t,ode.u,"Unknown",is_stochastic,[3])

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

function simulate_single_event(sm::StochasticModel)

    rates = sm.rn.compute_rates(sm.current_time,sm.current_population,sm.rn.rate_params)

    if sum(rates) == 0.0
        @warn "zero propagation!"
        return Inf,nothing
    end

    τ, which_event = find_event(rates)

    return τ, which_event

end

### FIXME
function round!(x::Vector)
    x[3] = round(x[3])
    return x
end

function deterministic_step(dm::Model)

    in_stochastic_region = false

    saveat_callback = (length(dm.options.saveat)==0) ? false : true

    # Check to see if the tumor has entered the stochastic regime
    enter_stochastic_region(u,t,integrator) = dm.stochastic_threshold-u[3]
    affect!(integrator) = begin
        in_stochastic_region = true
        terminate!(integrator)
    end
    cb1 = ContinuousCallback(enter_stochastic_region,affect!,nothing,save_positions=(false,saveat_callback))

    # Check to see when the tumor exceeds some c*original size. We define 
    # this as progression and assume the third compartment is the tumor compartment...
    # Potentially we can use a component array instead so we don't need to track the compartment
    progression(u,t,integrator) = u[3] - dm.progression_threshold*dm.initial_population[3]
    affect2!(integrator) = terminate!(integrator)
    cb2 = ContinuousCallback(progression,affect2!,save_positions=(false,saveat_callback))
    cbset = CallbackSet(cb1,cb2)

    prob = ODEProblem(dm.ODEmodel,dm.initial_population,dm.time_domain,dm.params,callback=cbset,
                        abstol=dm.options.abstol,reltol=dm.options.reltol,saveat=dm.options.saveat)

    odesol = solve(prob,Tsit5(),callback=cbset,abstol=1e-9,reltol=1e-6)

    return odesol,in_stochastic_region

end

function deterministic_step(dm::Model,sol::Solution)

    saveat_callback = (length(dm.options.saveat)==0) ? false : true

    # Check to see if the tumor has entered the stochastic regime
    enter_stochastic_region(u,t,integrator) = dm.stochastic_threshold-u[3]
    affect!(integrator) = begin
        sol.in_stochastic_region = true
        terminate!(integrator)
    end
    cb1 = ContinuousCallback(enter_stochastic_region,affect!,nothing,save_positions=(false,saveat_callback))

    # Check to see when the tumor exceeds some c*original size. We define 
    # this as progression and assume the third compartment is the tumor compartment...
    # Potentially we can use a component array instead so we don't need to track the compartment
    progression(u,t,integrator) = u[3] - dm.progression_threshold*dm.initial_population[3]
    affect2!(integrator) = terminate!(integrator)
    cb2 = ContinuousCallback(progression,affect2!,save_positions=(false,saveat_callback))
    cbset = CallbackSet(cb1,cb2)

    prob = ODEProblem(dm.ODEmodel,dm.initial_population,dm.time_domain,dm.params,callback=cbset,
                        abstol=dm.options.abstol,reltol=dm.options.reltol,saveat=dm.options.saveat)

    odesol = solve(prob,Tsit5(),callback=cbset,abstol=1e-9,reltol=1e-6)

    push!(sol.population_array,odesol.u[end])

    return sol

end

function simulate(dm::Model,rn::ReactionNetwork)

    # FIXME: We initially assume that the system is deterministic
    odesol,in_stochastic_region = deterministic_step(dm)

    sol = Solution(odesol,in_stochastic_region)

    if in_stochastic_region
        round!(sol.population_array[end])
    end
    current_time = sol.time_array[end]
    current_pop = sol.population_array[end]

    while current_time <= dm.time_domain[end]

        # Build the stochastic model
        sm = StochasticModel(current_time,current_pop,rn)

        # Get time to simulate forward and what occurs
        τ, which_event = simulate_single_event(sm)

        newparams = (in_stochastic_region = sol.in_stochastic_region, modelparams = dm.params.modelparams)

        # Simulate forward the deterministic components to the new time τ
        dm_2 = DeterministicModel(dm.ODEmodel,
                                    current_pop,
                                    newparams,
                                    dm.stochastic_threshold,
                                    dm.progression_threshold,
                                    (current_time,current_time+τ),
                                    dm.options)

        sol = deterministic_step(dm_2,sol)

        # solution with stochastic component and update time array
        sol.population_array[end] += rn.sto_mat[which_event]
        current_pop = sol.population_array[end]
        current_time += τ
        push!(sol.time_array,current_time)

        if floor(current_pop[3]) == 0.0
            sol.patient_outcome = "cure"
            break
        end
        
    end

    return sol

end
