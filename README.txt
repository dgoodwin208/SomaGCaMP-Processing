This is a set of scripts we used in the SomaGCaMP paper, still in the process of being properly commented for publication.

It uses NormCorre for motion correction, and CaImAn for putative neuron segmentation. Both are written and maintained by Eftychios Pnevmatikakis and Andrea Giovanucci at the Flatiron Institute.

NoRMCorre can be found here:
https://github.com/flatironinstitute/NoRMCorre

and CaImAn is available here:
https://github.com/flatironinstitute/CaImAn-MATLAB

I also attached a simple load3DTif_uint16 function for loading 3D tif data, which is also a wrapper around CaImAn functions for file I/O.  It assumes that the data is uint16 format.

If you want to do motion correction, the normcorre_demo.m should walk you through how to do it.

The processCNMF.m uses the CaImAn package to extract neurons and get resulting spatial footprints and time traces. It requires tuning a few parameters, but it should produce results from the first run. The CaImAn plot_components_GUI() feature is a great tool to see the spatial and temporal component of the extracted ROI.  

Once all the variables have been computed and are still in memory, you can use the line: 
plot_components_GUI(M1,A,C,b,f,Cn,options);

The postProcessCNMF uses the result to make the correlation heatmap and visualize the location of the neurons.

Finally, this is "research grade" code, which means it has been used by N=1 humans (Dan Goodwin), and so is surely going to have issues. If you hit any problems, don't hesitate to email me at dgoodwin@mit.edu!
 
