module WebinarPumasQSP20230321

using Distributed
addprocs(10, exeflags="--project=$(Base.active_project())")
@everywhere using JuliaSimModelOptimizer

using CSV
using DataFrames
using JuliaSimModelOptimizer: get_trials, get_model, get_search_space
using OrdinaryDiffEq: TRBDF2, remake, QNDF
using Plots

## Import the whole problem as a PEtab model
# using PumasQSP
using JuliaSimModelOptimizer
model_dir = joinpath(@__DIR__, "..", "Erdem_PLOSComputBiol2021")
petabyaml = joinpath(model_dir, "petab.yaml")
invprob = import_petab(petabyaml)

# Optional: updating the ODE solver
using Setfield
trials = []
sys = get_model(invprob)
for trial in get_trials(invprob)
    # t = remake(trial, sys, alg = QNDF(), abstol=1e-8, reltol=1e-8)  # , abstol=1e-12, reltol=1e-12)
    t = @set trial.config.alg = TRBDF2()
    push!(trials, t)
end
ss = get_search_space(invprob)
invprob = InverseProblem(trials, sys, ss)

# Running the optimization
vp = vpop(invprob, StochGlobalOpt(maxiters = 10000, parallel_type = JuliaSimModelOptimizer.EnsembleDistributed()),
          population_size=50)
CSV.write(joinpath(model_dir, "results", "vpop.csv"), vp)
vp = CSV.read(joinpath(model_dir, "results", "vpop.csv"), DataFrame)
vp = import_vpop(vp, invprob)

using ModelingToolkit
@variables t obs_pRecTot_free(t) obs_pAkt308_free(t) obs_pRPS6K_free(t) obs_pERK_free(t)
states = [obs_pRecTot_free, obs_pAkt308_free, obs_pRPS6K_free, obs_pERK_free]
p = plot(layout=(2, 2), dpi=600, size = (900, 500))
for trial in get_trials(invprob)
    if nameof(trial) == "c1"
        plot!(vp, trial , subplot=1, show_data=true, subplotstates=states,
              title = nameof(trial), legend=false, xlabel="", ylabel="Signal (AU)",
              bottom_margin=5*Plots.mm, left_margin=5*Plots.mm)
    elseif nameof(trial) == "c2"
        plot!(vp, trial , subplot=2, show_data=true, subplotstates=states,
              title = nameof(trial), legend=false, xlabel="", ylabel="",
              bottom_margin=5*Plots.mm, left_margin=5*Plots.mm)
    elseif nameof(trial) == "c3"
        plot!(vp, trial , subplot=3, show_data=true, subplotstates=states,
              title = nameof(trial), legend=false, xlabel="Time (s)", ylabel="Signal (AU)",
              left_margin=5*Plots.mm)
    else
        plot!(vp, trial , subplot=4, show_data=true, subplotstates=states,
              title = nameof(trial), legend=:topright, xlabel="Time (s)", ylabel="",
              left_margin=5*Plots.mm, legendfontsize=7)
    end
    display(p)
    savefig(joinpath(model_dir, "results", "alltrials" * ".png"))
end

## Other ways of importing models (not optimization problems)
# Importing SBML
using SBMLToolkit
sbmlfile = joinpath(model_dir, "model_sbml.xml")
rsys = readSBML(sbmlfile, ReactionSystemImporter())
sys = convert(ODESystem, rsys)

sys = readSBML(sbmlfile, ODESystemImporter())

sys = structural_simplify(sys)  # PumasQSP `petab_import` does that for you.

# Importing CellML
using CellMLToolkit
cellmlfile = joinpath(@__DIR__, "..", "CellML", "model_cellml")
sys = CellModel(cellmlfile)

# Import BioNetGen model
# add `writeSBML()` to the bottom of your BNGL file and run it
using SBMLToolkit
sbmlfile = joinpath(model_dir, "model_sbml.xml")
sys = readSBML(sbmlfile, ODESystemImporter())

# Create a ReactionSystem de novo
using Catalyst
@variables t S(t)=0.99 I(t)=0.01 R(t)=0.0
@parameters α=1e-4 β=1e-2
rsys = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end
sys = convert(ODESystem, rsys)

# Create an ODESystem de novo
using ModelingToolkit
@variables t S(t)=0.99 I(t)=0.01 R(t)=0.0
@parameters α=0.1 β=0.01
D = Differential(t)
eqs = [D(S) ~ -α*S*I,
       D(I) ~  α*S*I - β*I,
       D(R) ~  β*I]
sys = ODESystem(eqs, name=:SIR)

## Create Trials and InverseProblems from the imported or de novo created systems
data1 = CSV.read(joinpath(@__DIR__, "..", "SIR", "sirdata.csv"), DataFrame)
trial1 = Trial(data1, sys)
invprob = InverseProblem([trial1], sys, [α => (0.0, 0.2),  I=> (0.0, 0.02)])

end # module WebinarPumasQSP20230321
