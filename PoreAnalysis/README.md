# PoreAnalysis
Scripts to analyse the pore formation in the membranes. 

AnalysisLib.py
The basic library of the analysis. Functions include:
a. findcentre: find the centre of mass of a particular residue
b. findNeighour: find 4/8 neighbours of a particular cell in a 2D matrix
c. using First Depth Search find the longest connnected zero-valued cell

BilayerLib.py
The advanced library to calculate particular structural properties of the bilayer. require AnalysisLib to execute

SingleAPL.py - a script
Calculate the area per lipid (APL) of lower/upper leaflet at away and close to the protein region. 
Input: a given pdb
Method:
Get the 2D head group distribution of one leaflet - present as a 2D heatmap. Colour is the number of lipids at this local area
Look for the longest connected zero cells - should be where PC is replaced by the protein
The border cells of the area are used to calculate the APL close to the protein
The rest areas are used to calculate the APL away from the protein
APL.png: A graph illustration


SinglePore.py - a script
Calculate the diameter of the pore/ion channel at a certain frame
Input: a given pdb
Method:
Draw a water histogram along z-axis. Look for the narrowest region of the water bins. 
This water region is used for further analysis.
Apply PCA on this water cylinder. As the cylinder is round and flat. only the PC3 is representing the height of the cylinder.
The angle of the pore is calculated
nWAT * volWAT = pi* r^2 * height of the cylinder, where volWAT is from a pure WAT simulation at given temperature
diameter = r * 2

RealTimeAnalysis.py - a script
Used to get structural properties for several frames. 
Requirement:
AnalysisLib and BilayerLib
Input: image centred trajectory file, normal files required for cpptraj
Output:
results.dat - a summary of the structural properties
