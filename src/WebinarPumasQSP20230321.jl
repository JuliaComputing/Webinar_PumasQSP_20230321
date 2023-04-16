module WebinarPumasQSP20230321

using JuliaSimModelOptimizer
using SBMLToolkit

using CSV
using DataFrames
using Plots

model_dir = joinpath(@__DIR__, "..", "Erdem_PLOSComputBiol2021")
sbmlfile = joinpath(model_dir, "model_sbml.xml")

## Import BioNetGen model
# add `writeSBML()` to the bottom of your BNGL file and run it

## Import SBML
sys = readSBML(sbmlfile, ODESystemImporter())

using ModelingToolkit
sys = structural_simplify(sys)

using OrdinaryDiffEq
prob = ODEProblem(sys, [], (0, 1700))
sol = solve(prob)

## Import Data
trial1 = Trial(nothing, ssys)

## Import the whole problem as a PEtab model
petabyaml = joinpath(model_dir, "petab.yaml")
invprob = import_petab(petabyaml)

# using OptimizationPolyalgorithms
# vp = vpop(invprob, StochGlobalOpt(solver=PolyOpt(), maxiters = 20), population_size=20)
# CSV.write(joinpath(model_dir, "results", "vpop.csv"), vp)
vp = CSV.read(joinpath(model_dir, "results", "vpop.csv"), DataFrame)
vp = import_vpop(vp, invprob)

states = ["obs_pRecTot_free", "obs_pAkt308_free",
          "obs_pRPS6K_free", "obs_pERK_free"]
for trial in JuliaSimModelOptimizer.get_trials(invprob)
    # p = plot(vp, trial, show_data=true, title = nameof(trial), legend=:outertopright)
    p = plot(vp, trial; title = nameof(trial),
             states=states, legend=:outertopright)
    display(p)
    savefig(joinpath(model_dir, "results", nameof(trial) * ".png"))
end

end # module WebinarPumasQSP20230321
