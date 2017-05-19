clear all;

% Framework
%
%

% You can modify/remove the following "addpath" statements 
library = 'X:\Aquarids\code\chronux_code_2_00\chronux\spectral_analysis';
genpath library; 

if 0
library = 'X:\Aquarids\code\ftrip\fieldtrip-20100320\fieldtrip-20100320';
genpath library; 
fieldtripdefs;
end

% Load the crx preprocessed data into matlab session
%
cfg.load.dir = 'X:\Aquarids\code\loc2010';
cfg.load.filename = 'hdm1.lrx';
load([ cfg.load.dir filesep cfg.load.filename ], '-mat')


% Select recording segments
%

% Load preprocessing parameters
%

% Load localization parameters
%

% Load subject- & experimental-specific localization information 
%
% (1) Channel locations
% (2) Conductivity profile 
% (3) Conductivity delineating surfaces information 
%


% Load post-processing parameters
%

% Begin leadfield computations
% 

% Save core/unaltered leadfield 
%

% Begin data preprocessing computations 
%
% Temporal preprocessing: detrend, band-pass, artifact removal, averaging etc
%

% Spatial preprocessing: surface Laplacian 
% 

% Save preprocessed data


% Begin localization-specific preprocessing
%
% (1) Time or frequency domain 
% (2) Covariance matrices 

% Save localization-specific preprocessed data 

