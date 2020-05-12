using PowerDynamics: simulate, Perturbation, Inc, PowerGrid
using Distributions: Uniform
using Statistics
using NetworkDynamics
using LaTeXStrings

include("../new node types/plotting.jl")
include("../new node types/ThirdOrderEq.jl")
include("../../src/Random-Pertubation-on-Constrained-Manifolds/RandPertOnConstraindManifolds.jl")
include("../../src/PDpatches.jl")

"""
Evaluates the result of a simulation after a pertubation.
Checks if conditions for a stable synchronous state where met.
Conditions implemented so far:
1. Is ω returning to zero? (Realized through checking the deviations from 0)
#2. Did one of the machines get too far away from a stable state? -> Survivability
Inputs:
       pg: Power grid, a graph containg nodes and lines
       result: Solution of a power grid after a pertubation
       nodesarray: Array containg the index of all nodes that contain varible
       endtime: Endtime of the simulation
Outputs:
       Returns 1 if the system fullfills the condition for a stable state
       Returns 0 if the condition was not met
"""
function StableState(result::PowerGridSolution,nodesarray,variable,endtime, interval)
    for node in nodesarray
        summation = 0
        #=
        maxdeviation = maximum(result(0:endtime,node,:ω))
        if abs(maxdeviation) > threshold, (5)   # this could be used for Survivability later on...
            return 0
        end
        =#

        for time in endtime-10:0.1:endtime
            summation =+ abs.(result(time,node,:ω))
        end
        if summation > 100
            return 0
        end
    end

    # this is just an arbitrary value that i have choosen so far
    return 1
end

function sim(pg::PowerGrid, x0::State, timespan)
    problem = ODEProblem(rhs(pg),x0.vec,timespan)
    solution = solve(problem, Rodas5(autodiff=false), abstol=1e-8,reltol=1e-6,tspan=Inf)
    PowerGridSolution(solution, pg)
end


"""
Removes all nodes from the from the calculation procedure that do not contain
variable as a symbol.
Returns a list containg the indexes of all the remaning nodes.
Inputs:
    pg: Power grid, a graph containg nodes and lines
Outputs:
    nodesarray: Array containg the index of all nodes with the symbol varibale
"""
function RemoveNodes(pg::PowerGrid, variable::Symbol)
    nodesarray = findall(variable.∈ symbolsof.(pg.nodes))
    return nodesarray
end

"""
Produces a plot of the calulated basin stability against the nodes index.
Inputs:
    nodesarray: Array containg the index of all nodes that contain varible
    basinstability: Array containing the basin stability of each node
    standarddev: Array of the standard deviations of the basin stability
"""
function PlotBasinStability(nodesarray, basinstability, standarddev)
    plot(
        nodesarray,
        basinstability,
        #yerror = standarddev,
        #title = "...",
        xlabel = "Node Index",
        ylabel = "Basin Stability",
        marker = (:circle, 8, "red"),
        line = (:path, 2,"gray"),
        legend = false,
        grid = false,
        show = true
    )
end

"""
Calculates the single node basin stability of the nodes in a power grid using
PowerDynamics.jl.
The calculations were performed on the basis of the three following papers:
1. Menck, P., Heitzig, J., Marwan, N. et al. How basin stability complements
   the linear-stability paradigm. Nature Phys 9, 89–92 (2013)
2. Menck, P., Heitzig, J., Kurths, J. et al. How dead ends undermine
   power grid stability. Nat Commun 5, 3969 (2014).
3. C Mitra, A Choudhary, S Sinha, J Kurths, and R V Donnerk,
   Multiple-node basin stability in complex dynamical networks
   Physical Review E, (2017)

To calculate the basin stability "numtries" random pertubations are
drawn according to the interval.
The behaviour of the power grid after the pertubation is then calculated.
In the function StableState the number of times the system returns
to the stable synchronous state is counted.
If numtries -> ∞: Then BasinStability ≈ stablecounter / numtries.
Inputs:
    pg: Power grid, a graph containg nodes and lines
    endtime: Endtime of the simulation
    variable: Internal variable of the node types to be perturbed
    numtries: Number of random pertubations within interval to be evaluated
    interval: Lower and Upper Bound for the pertubation
Outputs:
    basinstability: Array conating the basinstability of each node
    standarddev: Array of the standard deviations of the basin stability
"""
function BasinStability(pg::PowerGrid, endtime::Int, variable::Symbol,
                        numtries::Int, interval::Vector{Float64})
    @assert interval[2] >= interval[1] "Upper bound should be bigger than lower bound!"

    operationpoint = find_operationpoint(pg, initial_guess(pg))
    nodesarray = RemoveNodes(pg,variable)
    basinstability = []
    standarddev =  []

    for node in nodesarray
        stablecounter = Array{Int64}(undef,numtries)
        # do some parralelization\ multithreading here to speed up calculation?
        for tries = 1:numtries
            println("NODE:",node," TRY:",tries)
            PerturbedState, dz = RandPertWithConstrains(powergrid, operationpoint,node, interval)

            λ, stable = check_eigenvalues(powergrid, PerturbedState)

            if stable != true
                println("Positive eigenvalue detected. The system might be unstable.")
                try
                    result = sim(pg, PerturbedState, (0.0, endtime))
                    display(plot_res(result, pg, node))
                    fn = "node$(node)_try$(tries)"
                    #png(fn)
                    stablecounter[tries] = StableState(result, nodesarray, variable, endtime, interval)
                    println("StableState Counter:  ",stablecounter[tries])
                catch err
                    if isa(err, GridSolutionError)
                        stablecounter[tries] = 0 # counted as an unstable solution
                        println("A GridSolutionError occured.")
                    end
                end
            else
                result = sim(pg, PerturbedState, (0.0, endtime))
                stablecounter[tries] = StableState(result,nodesarray,variable,endtime,interval)
            end
        end
        append!(basinstability, mean(stablecounter))
        append!(standarddev, std(stablecounter)) # is this even a useful measure here?
    end

    display(PlotBasinStability(nodesarray, basinstability, standarddev))
    #png("BasinStability_Plot")
    return basinstability, standarddev
end
