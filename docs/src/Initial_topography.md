# Generate Initial topography

Many carbonate platforms grow on a pre-exisiting abandoned carbonate platform. Therefore, it is reasonable to start the simulation with a platform-shaped initial topography.

The workflow to generate the topography is:
1) run the CarboKitten.jl to get the results.
2) export the sedimentation for each cell.
3) calculate the elvelation for each cell by adding subsidence. 

In the first step, the code is listed below. One thing should be noticed is that the grids are 100 by 70, with scale of 170 m. These values should be same as your runs lalter.

For the second step, the disintegration, production, deposition and sedimentation are exported respectively. The starting bathymetry for this run is set to be a slight slope (slope: 1/300)

For the third step, the elevation is calculated through by substracting the sedimentation with subsidence. 

The resultant initial topography that's ready for your run is stored in csv format. 

The next step is to import the csv file, through `initial_topography(path)`, where you may need to secify the location of the csv file. 

