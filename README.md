# GeoDT
# Developed by Luke P. Frash

This Geothermal Design Tool (GeoDT) is a fast multi-well flow and heat transfer model intended to aid high-level decision making 
for enhanced geothermal systems - geothermal energy development. This tool: 
(1)  generates a 3D geometry that includes wells and fractures
(2)  assigns dimensionally-scaled properties these wells and fractures
(3)  creates a mesh of 1D pipes and nodes to represent hydraulic connectivity in the 3D well and fracture network
(4)  solves this 1D network for fluid flow based on user assigned boundary conditions
(5)  predicts natural fracture and hydraulic fracture stimulation by fluid injection
(6)  solves this 1D network for time-dependent heat production
(7)  estimates transient net electrical power production from the network
(8)  outputs a csv file that summarizes the input and output parameters
(9)  outputs vtk files for visualizing the system geometry
(10) provides statistical data visualization example scripts and plots
* This code is in active development. We appreciate comments and questions that will help to improve this project.

For main/lfrash: 
+ watch: https://youtu.be/NgLch9r4XJQ
+ run using example files (e.g., "example_GLADE_2023.py")

For main/spio: operating instructions forthcoming 
+ watch: tbd
+ run using GeoDT_standard.py then select an settings file (e.g., settings "DEFAULT.gts")
+ copy, edit, or create new settings files to adjust the model to your needs (comma delimited text)
+ outputs in "setup_payoff.csv", vtk files module forthcoming
