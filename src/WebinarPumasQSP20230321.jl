module WebinarPumasQSP20230321

using PumasQSP
using SBMLToolkit

model_dir = joinpath(@__DIR__, "..", "Erdem_PlOSComputBiol2021")
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


end # module WebinarPumasQSP20230321
