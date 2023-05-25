###############################################################
# Load dependencies                                           #
###############################################################
using Distributed
addprocs(2, exeflags="--project=$(Base.active_project())")
@everywhere using JuliaSimModelOptimizer, OptimizationOptimJL, OptimizationBBO

using CSV
using DataFrames
using Dates
using JuliaSimModelOptimizer: get_trials, get_model, get_search_space,
                              get_initial_assignments, get_petab_problem,
                              export_petab, get_states
using OrdinaryDiffEq: TRBDF2, remake, solve
using Plots


###############################################################
# Import the whole problem as a PEtab model                   #
###############################################################
using PumasQSP
# using PumasQSP: import_petab
model_dir = joinpath(@__DIR__, "..", "Erdem_PLOSComputBiol2021")
petabyaml = joinpath(model_dir, "petab.yaml")  # Use petablint to validate the PEtab model.
invprob = PumasQSP.import_petab(petabyaml)

# Optional: updating the ODE solver for the trials
trials = []
sys = get_model(invprob)
for trial in get_trials(invprob)
    t = remake(trial, sys, alg = TRBDF2())  # Check out Chris JuliaCon2022 talk.
    push!(trials, t)
end
ss = get_search_space(invprob)
invprob = InverseProblem(trials, sys, ss)


###############################################################
# Running the optimization                                    #
###############################################################
algos = [SingleShooting(maxiters=10000, maxtime=100,
                        optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()),
         SingleShooting(maxiters=10000, maxtime=100,
                        optimizer=Optim.LBFGS()),
         MultipleShooting(maxiters=10000, maxtime=100, trajectories=2,
                        optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()),
         MultipleShooting(maxiters=10000, maxtime=100, trajectories=2,
                        optimizer=Optim.LBFGS()),
         MultipleShooting(maxiters=10000, maxtime=100, trajectories=4,
                        optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()),
         MultipleShooting(maxiters=10000, maxtime=100, trajectories=4,
                        optimizer=Optim.LBFGS()),
         MultipleShooting(maxiters=10000, maxtime=100, trajectories=8,
                        optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited()),
         MultipleShooting(maxiters=10000, maxtime=100, trajectories=8,
                        optimizer=Optim.LBFGS())]

using ParallelDataTransfer
sendto(workers(), invprob=invprob)
@everywhere function parres(alg)
    res = calibrate(invprob, alg)
    @info "Finished $alg"
    return res
end

@info "Starting parallel computation"
results = pmap(parres, algos[3:4])
@info "Finished parallel computation"

for (i, vp) in enumerate(results)
    dn = joinpath(dirname(@__FILE__), "plots", "$(now())_$(i)")
    mkpath(dn)
    df = DataFrame([vp.u], :auto)
    CSV.write(joinpath(dn, "res.csv"), df)

    for trial in get_trials(invprob)
        p = plot(trial, invprob, vp, show_data=true, title = nameof(trial), legend=:outertopright)
        savefig(joinpath(dn, nameof(trial) * ".png"))
    end
    p1 = convergenceplot(vp)
    savefig(joinpath(dn, "conv.png"))
end
@info "Finished writing files"


# df = CSV.read(joinpath("Erdem_PLOSComputBiol2021", "parameters.tsv"), DataFrame)
# lb = df[!, :lowerBound]
# ub = df[!, :upperBound]
# randnom = lb .+ (ub .- lb) .* rand(length(lb))
# nominal = df[!, :nominalValue]
# dfr = copy(df)
# dfr[:, :nominalValue] = randnom
# dfr
# CSV.write(joinpath("Erdem_PLOSComputBiol2021", "parameters_rand.tsv"), dfr, delim='\t')
