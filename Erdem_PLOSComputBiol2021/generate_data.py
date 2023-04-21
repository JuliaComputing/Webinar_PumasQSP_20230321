from amici.petab_simulate import PetabSimulator
import petab
import os

### Options #################################################
file_dir = os.path.dirname(__file__)
out_file = os.path.join(file_dir, 'measurements')
yaml_config = os.path.join(file_dir, 'petab.yaml')
#############################################################

# Simulation
petab_problem = petab.Problem.from_yaml(yaml_config)
simulator = PetabSimulator(petab_problem)
# simulation_df = simulator.simulate()
# simulation_df.to_csv(out_file + '_noisefree.tsv', sep='\t', index=False)
simulation_df = simulator.simulate(noise=True)
simulation_df.to_csv(out_file + '.tsv', sep='\t', index=False)
