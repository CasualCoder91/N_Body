import os # for relative paths

simulation_id = 1
timestep = 0
pixelfactor = 0.004 #depends on instrument used with scopesim -> print(cmd["!INST.pixel_scale"])

exposure_time = 3600 #[s]

#paths
output_base_path = os.path.join(os.path.abspath(__file__ + r"\..\..\..\.."), r"Output")
output_path = os.path.join(output_base_path, "Simulation" + str(simulation_id))
database_path = os.path.join(output_base_path,r"Database\Default.db")
fits_path = os.path.join(output_base_path, "Simulation" + str(simulation_id)+r"\scopesim_t"+str(timestep)+".fits")

#images
save_img = False
n_pixel = 14976 #Whole picture: 14976 | 1 fov: 4096 | test: 1000 | 1kM stars: 8.192

#analysis
#eps_magnitude = 0.00001