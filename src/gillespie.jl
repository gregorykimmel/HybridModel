module Gillespie

using Parameters, Plots

export Solution, GillespieParameters,simulate

mutable struct Solution{S,T<:Number,U}

    current_population::S
    current_time::T
    population_array::U  #Vector{Array{Int}}
    time_array::Vector{T}

end

Solution(pop,time) = Solution(pop,time,[pop],[time])

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

struct GillespieProblem{S<:Integer,T<:Number,R<:AbstractReactionNetwork}
    max_events::S                   # Maximum number of events that can occur
    initial_population::Vector{Int}
    time_domain::Tuple{T,T}
    rn::R
end

struct MonteCarloProblem
    n::Int
    p::GillespieProblem
end


function Plots.plot(sol::Solution,stride::Int=1)

    plot(
        sol.time_array[1:stride:end],
        hcat(sol.population_array[1:stride:end]...)',
        linecolor=:black,label=nothing
        )

end

function Plots.plot!(sol::Solution,stride::Int=1)

    plot!(
        sol.time_array[1:stride:end],
        hcat(sol.population_array[1:stride:end]...)',
        linecolor=:black,label=nothing
        )

end

function find_event(rates)

    # Get total to compute time to next event
    rate_total = sum(rates)

    # Compute tau (time till next event)
    τ = -log(rand())/rate_total

    # Normalize rates to determine which event occurs
    cum_rate = cumsum(rates./rate_total)

    if true
        # println(cum_rate)
    end

    # Roll the dice
    roll = rand()

    # find out which individual process occurred
    which_event = findfirst(x-> x >= roll,cum_rate)

    return τ, which_event

end

t₀(p::GillespieProblem) = p.time_domain[begin]
t₁(p::GillespieProblem) = p.time_domain[end]

function simulate(p::GillespieProblem)

    #=
    max_events::S                   # Maximum number of events that can occur
        initial_population::Vector{Int}
        time_domain::Tuple{T,T}
        reaction_network::R
    =#

    # Initialize the elements that will fill the solution struct
    current_population = p.initial_population
    current_time = t₀(p)

    # Initialize the solution array
    pop_array, time_array = Vector{Array{Int}}(),Vector{Float64}()

    push!(time_array,current_time)
    push!(pop_array,current_population)

    # Begin loop to run over Gillespie until max_events or final time is reached
    for event = 1 : p.max_events

        # println(current_population)

        # Compute population rates
        rates = p.rn.compute_rates(current_time,current_population,p.rn.rate_params)

        # Zero propagation
        if sum(rates) == 0
            break
        end

        τ, which_event = find_event(rates)

        if τ < 0
            @warn "something wrong"
            @show rates
        end

        # update time and population
        current_time += τ
        current_population += p.rn.sto_mat[which_event]

        push!(time_array,current_time)
        push!(pop_array,current_population)

        if current_time >= t₁(p)
            break
        end

    end

    return Solution(
        current_population,
        current_time,
        pop_array,
        time_array
        )


end

simulate(mcp::MonteCarloProblem) = map(1:mcp.n) do i
    run_gillespie(mcp.p)
end

end