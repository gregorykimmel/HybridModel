using Parameters

function plotfig(sol::solution,stride::Int=1)

    plot(
        sol.time_array[1:stride:end],
        vcat(sol.population_array...)[1:stride:end],
        linecolor=:black,label=nothing
        )

end

function plotfig!(sol::solution,stride::Int=1)

    plot!(
        sol.time_array[1:stride:end],
        vcat(sol.population_array...)[1:stride:end],
        linecolor=:black,label=nothing
        )

end

mutable struct solution{T<:Number}

    current_population::Vector{Int}
    current_time::T
    population_array::Vector{Array{Int}}
    time_array::Vector{T}

end

struct gillespie_parameters{T<:Number}
    max_events::Int                 # Maximum number of events that can occur
    final_time::Number              # Final time of simulation
    store_array::Bool               # Check whether we fill the pop/time Array
    initial_population::Vector{Int}
    initial_time::T
    sto_mat::Vector{Vector{Int}}    # How does population change?
    # compute_rates::Function         # output should be same length as sto_mat
end

# struct gillespie_computations



# end

function compute_rates(current_time,current_population)#,rate_params)

    λ = 0.2
    μ = 0.1
    K = 1000

    return λ*current_population[1]*(1.0 - current_population[1]/K), μ*current_population[1]

end


function find_event(rates)

    # Get total to compute time to next event
    rate_total = sum(rates)

    # Compute tau (time till next event)
    tau = -log(rand())/rate_total

    # Normalize rates to determine which event occurs
    cum_rate = cumsum(rates./rate_total)

    if true
        # println(cum_rate)
    end

    # Roll the dice
    roll = rand()

    # find out which individual process occurred
    which_event = findfirst(x-> x >= roll,cum_rate)

    return tau, which_event

end

function run_gillespie(parameters::gillespie_parameters)

    # unpack parameters into variables
    @unpack (max_events,final_time,store_array,initial_population,
     initial_time,sto_mat) = parameters

    # Initialize the elements that will fill the solution struct
    current_population = initial_population
    current_time = initial_time

    # Initialize the solution array
    pop_array, time_array = Vector{Array{Int}}(),Vector{Float64}()

    push!(time_array,current_time)
    push!(pop_array,current_population)

    # Begin loop to run over Gillespie until max_events or final time is reached
    for event = 1 : max_events

        # Compute population rates
        rates = compute_rates(current_time,current_population)#,rate_params)

        tau, which_event = find_event(rates)

        if which_event == nothing
            break
        end

        # update time and population
        current_time += tau
        current_population += sto_mat[which_event]

        push!(time_array,current_time)
        push!(pop_array,current_population)

        if current_time >= final_time
            break
        end

    end

    return solution(
        current_population,
        current_time,
        pop_array,
        time_array
        )


end



