include("../Random-Pertubation-on-Constrained-Manifolds/RandPertOnConstraindManifolds.jl")
include("../../example_cases/BasinStabilityForPD/pd_basin_stability.jl")
op = find_operationpoint(powergrid, initial_guess(powergrid))

ω_idx = findall(:ω .∈ symbolsof.(powergrid.nodes))

θ_idx = findall(:θ .∈ symbolsof.(powergrid.nodes))

PerturbedState, dz = RandPertWithConstrains(powergrid, op, 3, [-10,10])
check_eigenvalues(powergrid, PerturbedState)

result = sim(powergrid, PerturbedState, (0.0,50.0))
plot_res(result, powergrid, 3)

plot(PowerGridSolution(result.dqsol, powergrid), θ_idx, :θ, legend=false)


##
throws = 100
node = 5
v_vec = Array{Float64}(undef, throws)
θ_vec = Array{Float64}(undef, throws)
ω_vec = Array{Float64}(undef, throws)
for i in 1:throws
        PerturbedState, dz  = RandPertWithConstrains(powergrid, op, node, [-100,100])

        #sim(powergrid, PerturbedState,(0.0,20.0))
        v_vec[i] = abs(dz[3])
        θ_vec[i] = abs(dz[1])
        ω_vec[i] = PerturbedState[node,:ω]
end

scatter(θ_vec, v_vec, marker = (:circle, 8))
xlabel!("theta [rad]")
ylabel!("v [p.u.]")

###

pg = powergrid

u0 = op

timespan = (0.0,20.0)
ode = ODEProblem(rhs(pg), u0.vec, timespan, nothing)

function g(dx, x, e_s, e_d, p, t)
        i = total_current(e_s, e_d) + Y_n * (x[1] + x[2] * 1im)
        u = x[1] + x[2] * 1im
        θ = x[3]
        ω = x[4]
        resid[1] = (1im) * i * exp((-1im) * θ)
        resid[2] = (1im) * u * exp((-1im) * θ)
end

cb = ManifoldProjection(g)

# Something seems to be going wrong in the ForwardDiff package... the same problem as before the duals
sol = solve(ode,Rodas5(),save_everystep=false, callback = cb)
