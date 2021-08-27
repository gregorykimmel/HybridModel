include("Gillespie.jl")

using Plots
using DifferentialEquations

function compute_rates(current_time,current_population,params::Vector{T}) where {T<:Number}

    λA,λB,μA,μB,K0,σA,σB,Σ,ϕ = params;

    K = K0*(1.0 + Σ*sin(ϕ*current_time))

    tot_pop = sum(current_population)

    return (
        # max(λA*current_population[1]*(1.0 - sum(current_population)/K),0.0),
        # max(λB*current_population[2]*(1.0 - sum(current_population)/K),0.0), 
        λA*current_population[1],
        λB*current_population[2],
        μA*current_population[1] + λA*current_population[1]*tot_pop/K,
        μB*current_population[2] + λB*current_population[2]*tot_pop/K,
        σA*current_population[1],
        σB*current_population[2]
        )

end

function ode(du,u,p,t)

    λA,λB,μA,μB,K0,σA,σB,Σ,ϕ = p;
    A,B = u
    T = sum(u)

    K = K0*(1.0 + Σ*sin(ϕ*t))

    du[1] = (λA*(1.0 - T/K) - μA - σA)*A + σB*B
    du[2] = (λB*(1.0 - T/K) - μB - σB)*B + σA*A

end

sto_mat = [[1,0],[0,1],[-1,0],[0,-1],[-1,1],[1,-1]]
# sto_mat = [[1],[-1]]
max_events = 100000
initialPopulation = [50;10]
rate_params = [0.5,0.1,0.05,0.03,800,0.3,0.2,0.99,0.5]
# rate_params = [0.2,0.1,500,0.7]
time_domain = (0,200)
reaction_network = Gillespie.ReactionNetwork(sto_mat,compute_rates,rate_params)

par = Gillespie.GillespieProblem(max_events,initialPopulation,time_domain,reaction_network)

res = Gillespie.simulate(par)

prob = ODEProblem(ode,initialPopulation,time_domain,rate_params)
sol = solve(prob,Tsit5())