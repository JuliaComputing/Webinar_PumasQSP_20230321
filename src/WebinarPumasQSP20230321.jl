module WebinarPumasQSP20230321

###############################################################
# Import the whole problem as a PEtab model                   #
###############################################################
model_dir = joinpath(@__DIR__, "..", "Erdem_PLOSComputBiol2021")
petabyaml = joinpath(model_dir, "petab.yaml")  # Use petablint to validate the PEtab model.

using PumasQSP, ModelingToolkit
invprob_default = import_petab(petabyaml)

# Set ODE solver
using OrdinaryDiffEq
invprob = remake_trials(invprob_default, alg=TRBDF2())

# Check what we have imported
using Plots
for trial in get_trials(invprob)
    p = plot(trial, invprob, show_data=true,  # @Sebastian: initial conditions that are formulas of formulas cause troubles here.
             legend=:topleft, title=nameof(trial))
    savefig(joinpath(model_dir, "results", "default_$(nameof(trial)).png"))
    display(p)
end


###############################################################
# Running the optimization/Importing the results              #
###############################################################
# vp = vpop(invprob, StochGlobalOpt(maxiters = 10000),
#           population_size=50,
#           parallel_type = EnsembleDistributed())
# CSV.write(joinpath(model_dir, "results", "vpop.csv"), vp)
using CSV, DataFrames
vp = CSV.read(joinpath(model_dir, "results", "vpop.csv"), DataFrame)
vp = import_vpop(vp, invprob)

for trial in invprob.trials
    p = plot(vp, trial, title=nameof(trial), show_data=true, legend=:topleft)
    savefig(joinpath(model_dir, "results", "vp_$(nameof(trial)).png"))
    display(p)
end

end  # module