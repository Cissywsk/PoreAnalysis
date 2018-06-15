# PoreAnalysis
Scripts to analyse the pore formation in the membranes. 
<<<<<<< HEAD

AnalysisLib.py
 <br />The basic library of the analysis. Functions include:
 <br />a. findcentre: find the centre of mass of a particular residue
 <br />b. findNeighour: find 4/8 neighbours of a particular cell in a 2D matrix
 <br />c. using First Depth Search find the longest connnected zero-valued cell

BilayerLib.py
 <br />The advanced library to calculate particular structural properties of the bilayer. require AnalysisLib to execute

SingleAPL.py - a script
 <br />Calculate the area per lipid (APL) of lower/upper leaflet at away and close to the protein region. 
 <br />Input: a given pdb
 <br />Method:
 <br />Get the 2D head group distribution of one leaflet - present as a 2D heatmap. Colour is the number of lipids at this local area
 <br />Look for the longest connected zero cells - should be where PC is replaced by the protein
 <br />The border cells of the area are used to calculate the APL close to the protein
 <br />The rest areas are used to calculate the APL away from the protein
 <br />APL.png: A graph illustration 
 <br />

SinglePore.py - a script
 <br />Calculate the diameter of the pore/ion channel at a certain frame
 <br />Input: a given pdb
 <br />Method:
 <br />Draw a water histogram along z-axis. Look for the narrowest region of the water bins. 
 <br />This water region is used for further analysis.
 <br />Apply PCA on this water cylinder. As the cylinder is round and flat. only the PC3 is representing the height of the cylinder.
 <br />The angle of the pore is calculated
 <br />nWAT * volWAT = pi* r^2 * height of the cylinder, where volWAT is from a pure WAT simulation at given temperature
 <br />diameter = r * 2 

RealTimeAnalysis.py - a script
 <br />Used to get structural properties for several frames. 
 <br />Requirement:
 <br />AnalysisLib and BilayerLib
 <br />Input: image centred trajectory file, normal files required for cpptraj
 <br />Output:
 <br />results.dat - a summary of the structural properties
=======
Please go to folder PoreAnalysis for more detailed and tidy information
>>>>>>> e29f56efb21c82a40ef6f10aeeca169c8a8acc61
