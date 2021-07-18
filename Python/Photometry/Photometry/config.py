import os # for relative paths

simulation_id = 1
timestep = 1
pixelfactor = 0.004 #depends on instrument used with scopesim -> print(cmd["!INST.pixel_scale"])

exposure_time = 60 #[s]

#paths
output_base_path = os.path.join(os.path.abspath(__file__ + r"\..\..\..\.."), r"Output")
output_path = os.path.join(output_base_path, "Simulation" + str(simulation_id))
database_path = os.path.join(output_base_path,r"Database\Default.db")
fits_path = os.path.join(output_base_path, "Simulation" + str(simulation_id)+r"\scopesim_t"+str(timestep)+".fits")

#images
save_img = False
n_pixel = 4096 #Whole picutre: 14976 | 1 fov: 4096

#analysis
eps_magnitude = 0.05