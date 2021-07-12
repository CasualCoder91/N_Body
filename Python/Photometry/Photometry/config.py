import os # for relative paths

simulation_id = 1
timestep = 0
pixelfactor = 0.004 #depends on instrument used with scopesim -> print(cmd["!INST.pixel_scale"])

#paths
output_base_path = os.path.join(os.path.abspath(__file__ + r"\..\..\..\.."), r"Output")
output_path = os.path.join(output_base_path, "Simulation" + str(simulation_id))
database_path = os.path.join(output_base_path,r"Database\Default.db")
fits_path = os.path.join(output_base_path, "Simulation" + str(simulation_id)+r"\scopesim.fits")

#images
save_img = True
n_pixel = 14976