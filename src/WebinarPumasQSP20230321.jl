module WebinarPumasQSP20230321

###############################################################
# Load dependencies                                           #
###############################################################
using Distributed
addprocs(10, exeflags="--project=$(Base.active_project())")
@everywhere using JuliaSimModelOptimizer

using CSV
using DataFrames
using JuliaSimModelOptimizer: get_trials, get_model, get_search_space,
                              get_initial_assignments, get_petab_problem,
                              export_petab, report, get_states
using OrdinaryDiffEq: TRBDF2, remake, solve
using Plots


###############################################################
# Import the whole problem as a PEtab model                   #
###############################################################
using PumasQSP
model_dir = joinpath(@__DIR__, "..", "Erdem_PLOSComputBiol2021")
petabyaml = joinpath(model_dir, "petab.yaml")  # Use petablint to validate the PEtab model.
invprob = import_petab(petabyaml)

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
# vp = vpop(invprob, StochGlobalOpt(maxiters = 10000),
#           population_size=50,
#           parallel_type = EnsembleDistributed())
# CSV.write(joinpath(model_dir, "results", "vpop.csv"), vp)
vp = CSV.read(joinpath(model_dir, "results", "vpop.csv"), DataFrame)
vp = import_vpop(vp, invprob)


###############################################################
# Exporting the results                                       #
###############################################################
export_petab(petabyaml, vp)


###############################################################
# Visualizing the results in a PDF report                     #
###############################################################
report([vp], output=joinpath(model_dir, "reports"), title="Erdem_PLOSComputBiol2021", author="Jay Doe")  # Next release.


###############################################################
# Other ways of importing models (not optimization problems)  #
###############################################################
# Importing SBML
using SBMLToolkit
sbmlfile = joinpath(model_dir, "model_sbml.xml")
rsys = readSBML(sbmlfile, ReactionSystemImporter())

using ModelingToolkit
equations(rsys)

sys = convert(ODESystem, rsys)
equations(sys)

sys = readSBML(sbmlfile, ODESystemImporter())
ssys = structural_simplify(sys)  # PumasQSP `petab_import` does that for you.

observed(sys)  # Access the sys.observed
observed(ssys)


# Importing CellML
using CellMLToolkit
cellmlfile = joinpath(@__DIR__, "..", "CellML", "model_cellml.xml")
mdl = CellModel(cellmlfile)
sys = getsys(mdl)  # Extract the ODESystem from the CellModel
equations(sys)

ssys = structural_simplify(sys)
prob = ODEProblem(ssys, [], (0.0, 100.0))
sol = solve(prob)
plot(sol, idxs=(1,3))

# Import a MATLAB SimBiology(R) model
# Export the model as an SBML file and import it as above.
# Check if exported SBML is valid using an SBML compliant simulator of your choice.

# Import BioNetGen model
# add `writeSBML()` to the bottom of your BNGL file and run it to generate an SBML file.
using SBMLToolkit
sbmlfile = joinpath(model_dir, "model_sbml.xml")
sys = readSBML(sbmlfile, ODESystemImporter())

# Create a ReactionSystem de novo
using Catalyst
rsys = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end
sys = convert(ODESystem, rsys)

@variables t
@species S(t) I(t) R(t)
@parameters α β

ssys = structural_simplify(sys)
prob = ODEProblem(ssys, [S=>0.99 I=>0.01 R=>0.0], (0.0, 250.0), [α=>0.1, β=>0.01])
sol = solve(prob)
plot(sol)

# Create an ODESystem de novo
using ModelingToolkit
@variables t S(t)=0.99 I(t)=0.01 R(t)=0.0
@parameters α=0.1 β=0.01
D = Differential(t)
eqs = [D(S) ~ -α*S*I,
       D(I) ~  α*S*I - β*I,
       D(R) ~  β*I]
sys = ODESystem(eqs, name=:SIR)

ssys = structural_simplify(sys)
prob = ODEProblem(ssys, [], (0.0, 250.0))
sol = solve(prob)
plot(sol)


###############################################################
# Create Trials and InverseProblems from the                  #
# imported or de novo created systems                         #
###############################################################
data1 = CSV.read(joinpath(@__DIR__, "..", "SIR", "sirdata.csv"), DataFrame)
trial1 = Trial(data1, sys)
# trail2 = Trial(data2, sys, u0 => [I => 0.02], params = [β => 0.02])
invprob = InverseProblem([trial1], sys, [α => (0.0, 0.2),  I=> (0.0, 0.02)])


###############################################################
# Override a PEtab problem with your ODESystem                #
###############################################################
# Create an ODESystem (ignore the next two lines)
mdl = get_petab_problem(petabyaml).models[1]
initial_assignments = get_initial_assignments(mdl)
sys = ODESystem(mdl, defaults=initial_assignments)

# Override the PEtab problem with your ODESystem
invprob2 = import_petab(petabyaml, models=[sys])

# Optional: updating the ODE solver for the trials
trials = []
sys = get_model(invprob2)
for trial in get_trials(invprob2)
    t = remake(trial, sys, alg = TRBDF2())  # Check out Chris JuliaCon2022 talk.
    push!(trials, t)
end
ss = get_search_space(invprob2)
invprob2 = InverseProblem(trials, sys, ss)

vp1 = vpop(invprob2, StochGlobalOpt(maxiters = 1),
           population_size=1)  # Just to confirm that we can solve it
for trial in get_trials(invprob2)
    p = plot(vp, trial, show_data=true)  # Still the same plots
    display(p)
end


###############################################################
# Add dosing to your PEtab problem                            #
###############################################################
# Import the PEtab problem
invprob3 = import_petab(petabyaml)

# Create doses
sys = get_model(invprob3)
using UnPack
@unpack S1, S2 = sys
doses = [Bolus(; state = S1, amount = 100000, t = [300, 600]),
         Infusion(; state = S2, rate = 100, t = 1000, duration = 400)]

# Remake the trials with doses.
trials = []
for trial in get_trials(invprob3)
    t = remake(trial, sys, doses = doses, alg = TRBDF2())  # Next release.
    push!(trials, t)
end
ss = get_search_space(invprob3)
invprob3 = InverseProblem(trials, sys, ss)

end # module WebinarPumasQSP20230321
