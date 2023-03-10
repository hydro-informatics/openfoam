This chapter is developed by {{ scolari }}.

### Running the case

After having correctly created and filled in all the necessary dictionaries, the case can finally be run. 

* The first step consists in initializing a region containing water as defined in the *setFieldsDict* file.

```
user@user123:~/OpenFOAM-9/channel/Simulation$ setFields
```

* In the case of parallel runs, run the decomposePar command, as done during the meshing process, to decompose the geometry into individual geometries for each MPI process.

```
user@user123:~/OpenFOAM-9/channel/Simulation$ decomposePar

```

* Start running the simulation by typing the following in the terminal window:

    * In the case of parallel runs (substitute "x" with the number of cores):

```
user@user123:~/OpenFOAM-9/channel/Simulation$ mpirun -np x interFoam -parallel
```

    * Alternatively:

```
user@user123:~/OpenFOAM-9/channel/Simulation$ interFoam
```

Below, an example of the output screen of the *interFoam* solvers is shown. The main aspects to be mentioned are the following:

* The first line shows the mean and maximum flow {term}`CFL` condition.
* Line 2 shows instead the interface {term}`CFL` condition which is more restrictive than the previous and should be kept below 1 when solving multiphase flows.
* "MULES: Correcting alpha.water" depends on the value set to *nAlphaCorr* set in the *fvSolution* file.
* lines 7 and 9 refer to the *nAlphaSubCycles* which was set to 1, meaning only one loop.
* lines 11, 13 and 15 refer to the three pressure correctors and no non-orthogonal corrections are present.
* As shown in line 16 a tighter tolerance is applied only to this iteration (p-rghFinal)
  
```
1    Courant Number mean: 0.00233217 max: 0.961058
2    Interface Courant Number mean: 0.000313967 max: 0.234243
3    deltaT = 0.000461857
4    Time = 148.004
5 
6    smoothSolver:  Solving for alpha.water, Initial residual = 2.35103e-05, Final residual = 1.00306e-08, No Iterations 1
7    Phase-1 volume fraction = 0.144966  Min(alpha.water) = -1.2641e-05  Max(alpha.water) = 1
8    MULES: Correcting alpha.water
9    MULES: Correcting alpha.water
10   Phase-1 volume fraction = 0.144966  Min(alpha.water) = -2.71177e-05  Max(alpha.water) = 1
11   DICPCG:  Solving for p-rgh, Initial residual = 0.000249542, Final residual = 1.17508e-05, No Iterations 6
12   time step continuity errors : sum local = 8.05179e-08, global = 6.39537e-10, cumulative = -2.18471e-08
13   DICPCG:  Solving for p-rgh, Initial residual = 2.05956e-05, Final residual = 1.01515e-06, No Iterations 58
14   time step continuity errors : sum local = 6.95406e-09, global = -1.17766e-09, cumulative = -2.30247e-08
15   DICPCG:  Solving for p-rgh, Initial residual = 2.95498e-06, Final residual = 9.64318e-08, No Iterations 75
16   time step continuity errors : sum local = 6.60043e-10, global = 5.94012e-11, cumulative = -2.29653e-08
17   smoothSolver:  Solving for epsilon, Initial residual = 0.0002476, Final residual = 5.94123e-06, No Iterations 1
18   bounding epsilon, min: -3.63319e-06 max: 21370 average: 8.9958
19   smoothSolver:  Solving for k, Initial residual = 0.000103198, Final residual = 1.17097e-06, No Iterations 1
20   ExecutionTime = 63.45 s  ClockTime = 64 s
```

### Post-processing of the simulation results

In the case in which the simulations were run in parallel, before post-processing the data, the first step consists in reconstructing all solution steps of the analyzed case. This can either be done for all time steps or only for a specific one. The commands that need to be typed in the terminal window are shown below:

* To reconstruct all solution steps:

```
user@user123:~/OpenFOAM-9/channel/Simulation$ reconstructPar
```
  
* To reconstruct a specific time step (substitute "x" with the time step):

```
user@user123:~/OpenFOAM-9/channel/Simulation$ reconstructPar -time x
```

Once the case has been reconstructed, as for the meshing process, the following command should be used to visualize the case in ParaView:

```
user@user123:~/OpenFOAM-9/channel/Simulation$ paraFoam
```

The *channel.OpenFOAM* should now be present in the Pipeline Browser and to visualize it in the layout the *Apply* button can be used. Additionally, in the *Fields* section, the various fields that can be visualized are shown and can be selected/deselected according to the focus of the analysis.

```{figure} ../img/openfoam/interFoam/Paraview/channelOpenFOAM.png
:alt: openfoam 
:name: of-channelOpenFOAM

Visualization of the case results in ParaView.
```

In order to visualize the air and water phases, *alpha.water* should then be selected in the drop-down menu as shown in the image.

```{figure} ../img/openfoam/interFoam/Paraview/view-alpha-water.png
:alt: openfoam 
:name: of-view-alphawater

Enabling the setting for viewing the air and water phases in ParaView.
```

To change the shown time step, the arrows that can be seen in the area highlighted in red can be used.


```{figure} ../img/openfoam/interFoam/Paraview/final-time-step.png
:alt: openfoam timestep time step
:name: of-final-time-step

Options for changing the time step to be visualized.
```

Next, to visualize only the water phase, the *Clip* filter is used. This can either be found in the *Filters* section in the menu or alternatively the shortcut can be used. The *Clip Type* should be set to *Scalar*, selecting *alpha.water* as scalar and setting the value to 0.5, which represents the interface between air and water. To view the air phase the *Invert* option should be selected whereas for the water phase it should be deselected.

```{figure} ../img/openfoam/interFoam/Paraview/clip-water.png
:alt: openfoam clip water interFoam
:name: of-clip-water


Clip filter used for viewing the water phase in ParaView.
```

Finally, to also add the walls and patches to the view, the *Extract Block* filter can be implemented (click on the *channel.OpenFOAM* file before applying it).

```{figure} ../img/openfoam/interFoam/Paraview/extract-block.png
:alt: openfoam 
:name: of-extract-block

List of filters available in ParaView, highlighting ExtractBlock.
```

The patches of interest can the be either selected or deselected and the *Coloring* can be set to Solid Color.

```{figure} ../img/openfoam/interFoam/Paraview/choose-patches.png
:alt: openfoam 
:name: of-choose-patches

Available options for selecting the patches and changing the color.
```

The resulting view of the water phase and block extraction is shown below:

```{figure} ../img/openfoam/interFoam/Paraview/alpha-water.png
:alt: openfoam 
:name: of-alpha-water

Simulation results highlighting the water phase.
```

Different parameters can also be viewed, such as the flow velocity, and this can be done in the *Coloring* section by selecting *U*. The *preset* can be modified to better view the results by selecting the corresponding icon (highlighted in green).

```{figure} ../img/openfoam/interFoam/Paraview/flow-velocity.png
:alt: openfoam 
:name: of-flow-velocity

Simulation results highlighting the flow velocity.
```






