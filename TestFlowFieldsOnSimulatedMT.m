%% Simulate MT responses to velocity noise stimulus
% How well can we resolve the true spatiotemporal kernels for typical MT
% cells using a velocity noise stimulus? Motion reverse correlation
% presents a difficult problem: Motion IS correlation in space and time
% and correlations are bad for reverse correlation

% Several papers have invented different reverse correlation stimuli for
% measuring Spike-Triggered Average kernels in area MT

% Papers on the subject:

% 1. Borghuis, B. G., Perge, J. A., Vajda, I., van Wezel, R. J., van de Grind, W. A., & Lankheet, M. J. (2003). 
% The motion reverse correlation (MRC) method:: A linear systems approach in the motion domain. 
% Journal of neuroscience methods, 123(2), 153-166.

% 2. Perge, J. A., Borghuis, B. G., Bours, R. J., Lankheet, M. J., & van Wezel, R. J. (2005). 
% Temporal dynamics of direction tuning in motion-sensitive macaque area MT. 
% Journal of Neurophysiology, 93(4), 2104-2116.

% 3. Richert, M., Albright, T. D., & Krekelberg, B. (2013). 
% The complex structure of receptive fields in the middle temporal area. 
% Frontiers in systems neuroscience, 7.

% 4. Cui, Y., Liu, L. D., Khawaja, F. A., Pack, C. C., & Butts, D. A. (2013). 
% Diverse Suppressive Influences in Area MT and Selectivity to Complex Motion Features. 
% The Journal of Neuroscience, 33(42), 16715-16728.

% Papers 1 and 2 use a grid of points that steps in diferent directions on
% each frame (white in velocity), uniform over space... no spatial RF.
% Paper 3 uses a grid of drifting dots and reverse correlates over patches
% of motion... spatial resolution limited to grid size and temporal
% resolution is fixed to rate of update for the motion patches.

%% Set up some global parameters
clear all; close all; clc
frameRate    = 60; % this is what we've got
nTotalFrames = 10e3; % pick a big number... this is how much data you have
ny = 10;
nx = 10;
spatres = 1;

% stimulus parameters
temporalSmoothness = 2; % time constant of eponential (frames)
spatialSmoothness  = 3; % sigma of 2d gaussian (0 = white noise)


% simulated MT parameters
directionPreference        = 45;



%% Test Noisy Velocity Field
testNoisyFlowField

%% Play movie

%% Simulate MT's response
simulateMT

%% Do regression to recover the RF from spikes

% regression parameters

% ridge parameter
rho = 2e5; % bigger means more shrinkage (closer to STA than regression)
% number of principal components to use 
nPrincipalComponents = 1e3; 
nkt = 25;   % number of temporal frames in temporal RF

RecoverReceptiveField





