This repository contains the Matlab files to implement the G-ETM and G-ETMV methods
described in :

Robust point-process Granger causality analysis in presence of exogenous
temporal modulations and trial-by-trial variability in spike trains.

by Casile A., Faghih R. T. & Brown E. N.


For computational reasons, we use Matlab's parallel toolbox. If that is
not available the "parfor" commnands in the routines fitGLM_G_ETM, fitGLM_G_ETMV,
runGranger_G_ETM and runGranger_G_ETMV must be substituted with a "for".

------------------ Directory Matlab/G-ETM ------------------

The file runGranger_Examples_G_ETM.m runs two examples. Which of the
two examples is run depends on how the user sets the variable "dataSet"
in the code.

The first example runs our G-ETM Granger-causality method on a simple
network consisting of two units with one functional connection from
unit 1 to unit 2

The second example runs G-ETM on a data set consisting of 12 neurons recorded
from the monkey pre-motor cortex (Figs. 7 and S1 in the paper). This is a
computation-intensive example and it will take quite some time to complete.


------------------ Directory Matlab/G-ETMV ------------------

The file runGranger_Examples_G_ETMV generates the results plotted in Fig. 8
of the paper.


