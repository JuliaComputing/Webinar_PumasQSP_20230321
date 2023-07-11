module WebinarPumasQSP20230321

###############################################################
# Load dependencies                                           #
###############################################################
# using Distributed
# addprocs(10, exeflags="--project=$(Base.active_project())")
# @everywhere using JuliaSimModelOptimizer

# using CSV
# using DataFrames
# using JuliaSimModelOptimizer: get_trials, get_model, get_search_space,
#                               get_initial_assignments, get_petab_problem,
#                               export_petab, report, get_states
# using OrdinaryDiffEq: TRBDF2, remake, solve
# using Plots


###############################################################
# Import the whole problem as a PEtab model                   #
###############################################################
model_dir = joinpath(@__DIR__, "..", "Erdem_PLOSComputBiol2021")
petabyaml = joinpath(model_dir, "petab.yaml")  # Use petablint to validate the PEtab model.

using PumasQSP
# create_petab_template(; force=true)
invprob_default = import_petab(petabyaml)

# Set ODE solver
using OrdinaryDiffEq
invprob = remake_trials(invprob_default, alg=TRBDF2(),
                        abstol=1e-12, reltol=1e-9, dtmin=1e-15)  # Increasing tolerance will noticable chance results.

# Check what we have imported
using Plots
for trial in get_trials(invprob)
    p = plot(trial, invprob, show_data=true,  # @Sebastian: initial conditions that are formulas of formulas cause troubles here.
             legend=:topleft, title=nameof(trial))
    savefig(joinpath(@__DIR__, "results", "default_$(nameof(trial)).png"))
    display(p)
end


###############################################################
# Running the optimization                                    #
###############################################################
# vp = vpop(invprob, StochGlobalOpt(maxiters = 10000),
#           population_size=50,
#           parallel_type = EnsembleDistributed())
# CSV.write(joinpath(model_dir, "results", "vpop.csv"), vp)
vp = CSV.read(joinpath(model_dir, "results", "vpop.csv"), DataFrame)
vp = import_vpop(vp, invprob)

for trial in invprob.trials
    p = plot(vp, trial, title=nameof(trial), show_data=true, legend=:topleft)
    savefig(joinpath(@__DIR__, "results", "vp_$(nameof(trial)).png"))
    display(p)
end

end  # module