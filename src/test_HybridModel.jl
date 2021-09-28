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

	V = 5.0*10^6        # Volume of blood in the body
	VT = 10^9           # Average number of cancer cells per cm^3
	freqInPeriph = 0.01

	params=(rN = 0.16,	# normal T cell net growth rate
	kN = 500.0*V/freqInPeriph,			# normal T cell carrying capacity
	λC = 0.037,			# basal CAR growth rate
	kC = 139.0*V/freqInPeriph, #139.0,			# CAR T cell carrying capacity
	a  = 0.423,			# signal inefficiency factor
	b  = 0.525,			# immune reconstition factor
	λB  = 0.25,			# tumor birth rate (0.03 - 0.28) Wilkins 2000
	δB = 0.13,
	γB = 1.25,			# tumor killing rate (0.85-1.35) Kimmel
	kB = 4.05*V/freqInPeriph,			# Half-maximal tumor-killing rate Kimmel
	δZ = 0.5,
	α = 0.01,
	β = 0.01,
	ϕ = 43.3/VT,
	θ = 1.137/VT,
	τ = 0.267*(1.0 + 0.25*randn()) #0.267
	)

	sto_mat = [[0,0,1,0,0],[0,0,-1,0,0]]

	t_current = 0.0
	t_final = 365.0

	N0,C0,B0,Z0,I0 = 6.0*V/freqInPeriph,0.36*V/freqInPeriph,94.84*VT,1e2*V/freqInPeriph,0

	u0 = [N0;C0;B0;Z0;I0]

	possible_stochastic_compartments = [3]

	opts = Options(abstol=1e-15,reltol=1e-12)
	
	u_current = u0

	dmparams = (modelparams=params,in_stochastic_region=false)

	# Build the deterministic model step
	deterministic_model = DeterministicModel(ctDNAmodel,u_current,dmparams,100,1.2*u0[3],(t_current,t_final),
	possible_stochastic_compartments,opts)

	reaction_network = ReactionNetwork(sto_mat,compute_rates,params)

	sol = simulate(deterministic_model,reaction_network)

	if sol.patient_outcome == "cure"
		if sol.population_array[end][5] > 1
			days_till_not_pet_pos = log(sol.population_array[end][5])/params.τ
			println("No longer PET+ after $(ceil(sol.time_array[end] + days_till_not_pet_pos))")
		else
			println("PET- on day $(ceil(sol.time_array[end]))")
		end
	else
		detectable_idx=findall(z->z[3]/VT > 1,sol.population_array)
		detectable_time = sol.time_array[[idx for idx in diff(detectable_idx) if idx > 1]][1]
		println("Progression at day $(ceil(detectable_time))")
	end

	return sol
end