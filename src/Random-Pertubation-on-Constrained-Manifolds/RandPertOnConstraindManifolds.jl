using PowerDynamics
using ForwardDiff
using Distributions:Uniform
import PowerDynamics:rhs
using DifferentialEquations
using LinearAlgebra
using DiffEqCallbacks
include("../PDpatches.jl")
include("system.jl")


"""
A method to find random pertubations in the case of systems with algebraic
constraints. There is then a set of valid initial conditions which lies on a so
called constrained manifold.
The system must still lie on this manifold after the random pertubation.
This could be applied to Basin Stability Methods in the case of complexer
Generator Methods which contain algebraic constraints.
"""
function RandPertWithConstrains(pg::PowerGrid, operationpoint,node, interval)
    g(φ,i,u) = [1im * i * exp(-1im * φ), 1im * u * exp(-1im * φ)] # Constraints
    #g(θ,i,u) = [1im * u * exp(-1im * θ)] # Constraints

    z = [operationpoint[node,:φ],operationpoint[node,:iabs],operationpoint[node,:v]]

    Jacg_imag = ForwardDiff.jacobian(x -> imag(g(x[1],x[2],x[3])), z)
    Jacg_real = ForwardDiff.jacobian(x -> real(g(x[1],x[2],x[3])), z)

    Jacg = Jacg_real + 1im .* Jacg_imag  # I could not find a package which is able to handle complex valued function

    #Jacg = grad = ForwardDiff.gradient(x -> real(g(x[1],x[2],x[3])), z) .+ 1im .* ForwardDiff.gradient(x -> imag(g(x[1],x[2],x[3])), z)

    kerJacg = nullspace(Jacg) # kernel of Jacg, returns orthogonal vectors

    Pro_kerg = Projvw(z, kerJacg) # Projection of z on ker(∇g(z))

    Frand = rand(Uniform(interval[1],interval[2]), length(z)) # random pertubation "force"

    dz = Pro_kerg .* Frand # Random Pertubation along the Projection of the Constraind Manifold
    PerturbedState = copy(operationpoint) # otherwise operationpoint is overwritten
    PerturbedState[node, :θ] = abs(dz[1]) # thats not correct but its a quick workaround for now
    PerturbedState[node, :u] = dz[3]

    return PerturbedState,dz
end

"""
Calculates and returns the Projection of a vector w on the vector v.
"""
function Projvw(v,w)
    dot(w,v) * v * (1 / norm(v)^2)
end

function rhs(pg::PowerGrid)
    network_dynamics(map(construct_vertex, pg.nodes), map(construct_edge, pg.lines), pg.graph)
end

"""
A method to assures that the conditions given from g stay on the contrained
manifold. A DifferentialEquations Callback is used.
This DOES NOT work because of a type missmatch regarding complex numbers
Inputs:
    pg: Power grid, a graph containg nodes and lines
    endtime: Endtime of the simulation
    u0: current state of the powergrid pg
Outputs:
    solution: a PowerGrid solution containg the timeseries and the Powergrid pg
"""
function SimManifoldCd(pg::PowerGrid, u0::State, endtime)
    ode = ODEProblem(rhs(pg), u0.vec, (0.0, endtime), nothing)

    function g(resid, x, e_s, e_d, p, t)
            i = total_current(e_s, e_d) + Y_n * (x[1] + x[2] * im)
            u = x[1] + x[2] * im
            θ = x[3]
            ω = x[4]
            resid[1] = (1im) * i * exp((-1im) * θ)
            resid[2] = (1im) * u * exp((-1im) * θ)
    end
    cb = ManifoldProjection(g)


    dqsol = solve(ode,Rodas5(autodiff=false),save_everystep=false,callback=cb)
    solution = PowerGridSolution(dqsol, pg)
    return solution
end
