using Parameters, ComponentArrays

# include("Gillespie.jl")
# include("HybridModel.jl")

function ctDNAmodel(du,u,p,t)

	# Grab dependent variables from soln array
	N,C,B,Z,I = u

	in_stochastic_region, params = p.in_stochastic_region,p.modelparams

	@unpack rN,λC,b,a,kN,kC,kB,λB,γB,θ,δZ,δB,α,β,ϕ,τ = params

    # Tumor-killing function
	γC = γB*C/(kB + C)

	T = N + C
	rC = λC + b*(T - kN)^2/(a*T^2 + (T - kN)^2)

	# RHS
	du[1] = -rN*N*log(T/kN)
	du[2] = -rC*C*log(T/kC)
	du[3] = in_stochastic_region ? 0.0 : (λB - δB)*B - γC*B
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

function main()

	params=(rN = 0.16,	# normal T cell net growth rate
	kN = 500.0,			# normal T cell carrying capacity
	λC = 0.037,			# basal CAR growth rate
	kC = 119.0, #139.0,			# CAR T cell carrying capacity
	a  = 0.423,			# signal inefficiency factor
	b  = 0.525,			# immune reconstition factor
	λB  = 0.25,			# tumor birth rate (0.03 - 0.28) Wilkins 2000
	γB = 1.15,			# tumor killing rate (0.85-1.15) Kimmel
	kB = 4.05,			# Half-maximal tumor-killing rate Kimmel
	δZ = 0.1,
	δB = 0.1,
	α = 0.1,
	β = 0.1,
	ϕ = 0.1,
	θ = 10.0,
	τ = 0.1)

	sto_mat = [[0,0,1,0,0],[0,0,-1,0,0]]

	t_current = 0.0
	t_final = 200.0
	u0 = [6.0;0.36;1e6;1e2;0]

	opts = Options()

	u_current = u0

	dmparams = (modelparams=params,in_stochastic_region=false)

	# Build the deterministic model step
	deterministic_model = DeterministicModel(ctDNAmodel,u_current,dmparams,100,1.2*u0[3],(t_current,t_final),opts)

	reaction_network = ReactionNetwork(sto_mat,compute_rates,params)

	sol = simulate(deterministic_model,reaction_network)

	return sol
end