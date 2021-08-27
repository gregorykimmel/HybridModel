using Parameters, ComponentArrays

# include("Gillespie.jl")
include("HybridModel.jl")

function ctDNAmodel(du,u,params,t)

	# Grab dependent variables from soln array
	N,C,B,Z,I = u
	@unpack rN,λC,b,a,kN,kC,kB,λB,γB,θ,δZ,δB,α,β,ϕ,τ = params

    # Tumor-killing function
	γC = γB*C/(kB + C)

	T = N + C
	rC = λC + b*(T - kN)^2/(a*T^2 + (T - kN)^2)

	# RHS
	du[1] = -rN*N*log(T/kN)
	du[2] = -rC*C*log(T/kC)
	du[3] = (λB - δB)*B - γC*B
	du[4] = θ*(α*δB + β*γC)*B - δZ*Z		# Normal T cells is in there δz
	du[5] = ϕ*γC*B - τ*I

end

function compute_rates(t,u,params)

    N,C,B,Z,I = u
    @unpack rN,λC,b,a,kN,kC,kB,λB,γB,θ,δZ,δB,α,β,ϕ,τ = params

    # Get total Tcell population
    T = N+C

    # Tumor-killing function
	γC = γB*C/(kB + C)

	T = N + C
	rC = λC + b*(T - kN)^2/(a*T^2 + (T - kN)^2)

    # du[1] = -rN*N*log(T/kN)
	# du[2] = -rC*C*log(T/kC)
	# du[3] = (λB - δB)*B - γC*B
	# du[4] = θ*(α*δB + β*γC)*B - δZ*Z		# Normal T cells is in there δz
	# du[5] = ϕ*γC*B - τ*I

    return (
        λB*B,
        (δB + γC)*B
        )

end

# function main()

    params=(rN = 0.16,		# normal T cell net growth rate
    kN = 500.0,		# normal T cell carrying capacity
    λC = 0.037,	# basal CAR growth rate
    kC = 139.0,		# CAR T cell carrying capacity
    a  = 0.423,		# signal inefficiency factor
    b  = 0.525,		# immune reconstition factor
    λB  = 0.25,		# tumor birth rate (0.03 - 0.28) Wilkins 2000
    γB = 1.15,		# tumor killing rate (0.85-1.15) Kimmel
    kB = 4.05,		# Half-maximal tumor-killing rate Kimmel
    δZ = 0.1,
    δB = 0.1,
    α = 0.1,
    β = 0.1,
    ϕ = 0.1,
    θ = 10.0,
    τ = 0.1)

    sto_mat = [[1],[-1]]

    u0 = [6.0;0.36;1e6;1e2;0]

    # struct DeterministicModel{P,T,U<:Number} <: Model
    #     ODEmodel::Function              # Deterministic version of the Model
    #     stochastic_threshold::Number    # When to switch to/from the StochasticModel
    #     progression_threshold::Number   # A factor that determines progression
    #     initial_population::Vector{T}
    #     params::P                       # Parameters to be fed to the ODEmodel
    #     time_domain::Tuple{U,U}         # Max length of simulation go. (typically (0, tf) where tf is final time)
    # end

    t_final = 365.0
    t_current = 0.0
    u_current = u0

    # Build the deterministic model step
    model = HybridModel.DeterministicModel(ctDNAmodel,100,1.5*u0[3],u_current,params,(t_current,t_final),false)
    sol = HybridModel.deterministic_step(model)

    # The reaction network that describes the transitions of the populations
    reaction_network = Gillespie.ReactionNetwork(sto_mat,compute_rates,model.params)

    # # If the population never goes stochastic, we won't enter into the hybrid model part
    # while t_current < t_final

    #     in_stochastic_region = true

    #     # Build the stochastic model step
    #     model = Gillespie.GillespieProblem(1,round.(Int,u),(t_current,t_final),reaction_network)
        
    #     # Get new time and which event occurred in the stochastic compartment
    #     τ,which_event = Gillespie.simulate_single_event(model)

    #     # Simulate forward the ODE (holding the stochastic compartment constant)
    #     # 
    #     HybridModel.DeterministicModel(ctDNAmodel,100,1.5*u0[3],u_current,params,min(t_current+τ,t_final))
    #     sol = HybridModel.deterministic_step(model)


    # end

#     return sol

# end