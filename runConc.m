clear all

[cfg.load.filename, cfg.load.dir] = uigetfile('*.lfd', 'Select lfd file'); %AG edit
% cfg.load.dir = 'C:\Documents and Settings\Andrew Goldfine\My Documents\EEGResearch\LBEX\loc2010'; %AG edit
% cfg.load.filename = 'hdm1.lfd';
% cfg.load.filename = 'N_R1.lfd';
cfg.save.dir=pwd; %AG edit
% cfg.save.dir = cfg.load.dir ;
[pathstr, name, ext] = fileparts( cfg.load.filename ); 
cfg.save.filename = [name '.cev'];
%cfg.save.filename = 'hdm1.con';
cfg.show = 0;

% Load the leadfield (lf) & its configuration file (lfcfg)
% 
load([ cfg.load.dir filesep cfg.load.filename ], '-mat')

kernel = lf' ;
voxelCentroids = lfcfg.srcSpace.grid( lfcfg.srcSpace.inside, : );

% Set the discrete-discrete concentration parameters 
%
conc.roiVolume = 0.25; % as percentage of overall head volume
% conc.cutOff = 0.75 / 100; 
% conc.threshold = 0.95;
conc.dir = cfg.save.dir;
conc.filename = cfg.save.filename;
conc.show.eigenvalues = 0;  % Show conc eig plots
conc.show.msg = 0;

conc = ddConcentration( kernel, conc, voxelCentroids );
%params.conc = conc;

cfg.date = datestr(clock);

