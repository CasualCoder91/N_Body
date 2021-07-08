import os # for relative paths

simulationID = 1

#paths
outputBasePath = os.path.join(os.path.abspath(__file__ + r"\..\..\..\.."), r"Output")
outputPath = os.path.join(outputBasePath, "Simulation" + str(simulationID))
databasePath = os.path.join(outputBasePath,r"Database\Default.db")
fitsPath = os.path.join(outputBasePath, "Simulation" + str(simulationID)+r"\scopesim.fits")
