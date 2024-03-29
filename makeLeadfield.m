
clear all;

% Leadfield Library 
%
library = '/home/chronux/cornell/data/lbex_data/loc2010/leadfields'; %AG edit
addpath( genpath( library ) );


% Load from & Save to ...
%
[lfcfg.load.filename, lfcfg.load.dir] = uigetfile('*.prx', 'Select prx file'); %AG edit
% lfcfg.load.dir = 'C:\Documents and Settings\Andrew Goldfine\My
% Documents\EEGResearch\LBEX\loc2010'; %AG edit but think it's wrong so choose above
% lfcfg.load.filename = 'hdm1.prx';

% lfcfg.load.filename = 'N_R1.prx'; %AG - this file exists in loc2010
% folder, but not likely what I want
lfcfg.save.dir=pwd;
% lfcfg.save.dir = lfcfg.load.dir ;
[pathstr, name, ext] = fileparts( lfcfg.load.filename ); 
lfcfg.save.filename = [name '.lfd'];
%lfcfg.save.filename = 'hdm1.lfd';

% Display intermediate results? 
%
lfcfg.show = 1; % Display intermediate plots & messages
hdm.displayMsg = 1;
srcSpace.plot = 1; srcSpace.display = 1;

% Set all the parameters
%

% Re-reference the leadfield as per EEG data referencing
% 'average' or 'channel' for a specific channel
%
lfcfg.refScheme = 'average'; % or 'channel' for a specific channel

% Source/Dipole space parameters
%
srcSpace.units = 'cm';
srcSpace.voxelSize = 2.5; % This to some degree controls resolution

scalpRadius = 8.0; % units = cm
% Typical scalp conductivity = 0.33 (ohm-meter)^(-1) = 0.0033 (ohm-cm)^(-1)
scalpConductivity = 0.0033; % units = 1/(ohm-cm)

% Can use preset models or user provided models
% For both, *.model.radii & *.model.conductivity, the 
% options are: 'rush-driscoll' or 'cuffin-cohen' or 'stok' or 'user' 
% However, I'd suggest that we stay with 'user' option since conductivity
% normalization issue can get confusing. 
%
hdm.model.radii = 'user';   
hdm.model.conductivity = 'user';
% If 'user', specify radii & conductivity ratios relative to scalp
% 1=scalp, 2=bone, 3=CSF fluid, 4=cortex
% 
r(1) = 1.0; r(2) = 0.9467; r(3) = 0.8667; r(4) = 0.84;
c(1) = 1.0; c(2) = 0.0125; c(3) = 3.0; c(4) = 1.0;


% Load the data file to get the sensor locations - nothing else reqd from
% data file
% 
load([ lfcfg.load.dir filesep lfcfg.load.filename ], '-mat')
sensor.Location = ts.channel.location;
system = ts.channel.system;


% Change sensor locations to centimeters since standard 10-20 system file
% we used gave sensor locations based on a sphere of radius = 0.1 meters
%
if strcmpi(system.units, 'meters'), sensor.Location = 100 * sensor.Location; end


% Check whether the sensor & voxel Cartesian systems are the same
%
% What coordinate system are the sensor locations reported in? 
% -- Sphere of radius 0.1 meters
% -- They are Cartesian, (x,y,z)
% -- The +z axes passes thru the top of the head, channel Cz=(0,0,0.1)
% -- The +x axes passes thru extreme right channel T10=(0.0865,0,-0.05)
% -- The +y axes passes thru nose-front channel FPz=(0,0.1,-0.0021)
% Therefore in summary, +x towards right ear, +y towards nose, +z towards
% head top. 

% --> This won't matter since voxelizeSphere.m will utilize sensor coord
% system, ie set head coordinates based on sensor system


% Head model parameters -- see help for eegDipoleInSphere3d3.m
% 4-shell spherical model: r=radii, c=conductivity. 
% Conductivity is assumed isotropic (no difference between tangential &
% radial conductivities). 
% 
% Note: A 3-shell model can be simulated by merging 2 layers, ie 
% specify 2 radii and conductivities the same. 
%
% scalpRadius, scalpConductivity ARE dimensioned units. 
% The output should be interpreted carefully. If length unit in meters,
% conductivity in 1/(ohm-meter), the unit dipole can be interepreted as
% ampere-meter and potential (leadfield) as volts(V) (or alternatively, the
% unit-dipole as mA-meter and potential (leadfield) as mV, here m=milli).
% If the units of length are in centimeter (cm), for the example above, the
% unit-dipole will be either A-cm or mA-cm, potential (leadfield) remains
% the same. 
%

% scalpRadius set above or here
% scalpRadius = 8.0; % units = cm

% scalpConductivity set above or here
% Typical scalp conductivity = 0.33 (ohm-meter)^(-1) = 0.0033 (ohm-cm)^(-1)
% scalpConductivity = 0.0033; % units = 1/(ohm-cm)

% Can use preset models or user provided models
% For both, *.model.radii & *.model.conductivity, the 
% options are: 'rush-driscoll' or 'cuffin-cohen' or 'stok' or 'user' 
% However, I'd suggest that we stay with 'user' option since conductivity
% normalization issue can get confusing. 
%
% hdm.model.radii = 'user';   
% hdm.model.conductivity = 'user';

% Set above or here
% If 'user', specify radii & conductivity ratios relative to scalp
% 1=scalp, 2=bone, 3=CSF fluid, 4=cortex
% 
% r(1) = 1.0; r(2) = 0.9467; r(3) = 0.8667; r(4) = 0.84;
% c(1) = 1.0; c(2) = 0.0125; c(3) = 3.0; c(4) = 1.0;


hdm.radius = r * scalpRadius; 
hdm.conductivity = c * scalpConductivity;
params.eegHeadModel = hdm; 
% hdm.displayMsg = 1;
lfcfg.hdm = hdm;


% Source/Dipole space
%
srcSpace.gridLimits(:,1) = [-hdm.radius(1),hdm.radius(1)];
srcSpace.gridLimits(:,2) = [-hdm.radius(1),hdm.radius(1)];
srcSpace.gridLimits(:,3) = [-hdm.radius(1),hdm.radius(1)];
srcSpace.radius = hdm.radius(1);
srcSpace.origin = [0,0,0];
% srcSpace.plot = 0; srcSpace.display = 1;

% Set above or here
% srcSpace.units = 'cm';
% srcSpace.voxelSize = 1.5; % This to some degree controls resolution

% srcSpace.grid( :, 3 ) = Cartesian coordinates of the voxel centroids of entire cube.
% srcSpace.inside = Array containing indices of voxels within sphere(s) radius (radii).
% srcSpace.outside = Array containing indices of voxels outside sphere(s) radius (radii).
% srcSpace.nVoxels = Total number of voxels within or on the sphere. 
%
srcSpace = voxelizeSphere( srcSpace );
% Transfer the voxelization data too, since centroids will be required
% during the concentration problem evaluation.
%
lfcfg.srcSpace = srcSpace; 


% Evaluate the leadfield
%
tic, [lf, tmp] = eegDipoleInSphere3d3( sensor.Location, srcSpace.grid(srcSpace.inside, :), params ); toc


% Re-reference the leadfield as per EEG data referencing
%
% Set above or here
% lfcfg.refScheme = 'average'; % or 'channel' for a specific channel

if strcmpi( lfcfg.refScheme, 'average')
    reflf = mean(lf,2); % average across channels 
elseif strcmpi( lfcfg.refScheme, 'channel') 
    reference = ts.channel.reference;
    if isempty(reference.location), error('Reference channel location is [].'); end
    [reflf, tmp] = eegDipoleInSphere3d3( reference.location, srcSpace.grid(srcSpace.inside, :), params );
else
    error( [mfilename ' need to specify a referencing scheme'] );
end
lf = lf - repmat(reflf, 1, size(lf,2));


% Save to file
%
lfcfg.date = datestr(clock);
save([ lfcfg.save.dir filesep lfcfg.save.filename ], 'lf', 'lfcfg');


% Visualize
%
if lfcfg.show
    figure; imagesc( lf(1:3:end,:)); colorbar('vert'); title( ['eegDipoleInSphere3d:: X-Leadfield Matrix'] );
    figure; imagesc( lf(2:3:end,:)); colorbar('vert'); title( ['eegDipoleInSphere3d:: Y-Leadfield Matrix'] );
    figure; imagesc( lf(3:3:end,:)); colorbar('vert'); title( ['eegDipoleInSphere3d:: Z-Leadfield Matrix'] );
%    imagesc( lf ); colorbar('vert'); title( ['eegDipoleInSphere3d:: Leadfield Matrix'] );
end

if lfcfg.show, sensorDistribution(sensor.Location, srcSpace.grid(srcSpace.inside, :)); end
