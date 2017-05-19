% Script for extracting and Non-analytic preprocessing of crx data 
%
% (1) Load the data into structures suitable for localization 
% (2) Load the channel relabelling xml file & relabel
% (3) Remove non-neurological channels, eg EKG, EMG, EOG, etc 
% (4) Remove any bad channels, eg with drift or noticeable artifacts etc 
% (5) Add on the reference channel (and get it's position too)
% (6a) Match the current data labels with the 10-20 standard labels, get the
%     electrode location 
% (6b) Check whether coordinate system for location specification is the
%      same (consistent) with voxelizeSphere.m, or else leadfield evaluation
%      will be incorrect. 
% (6c) Group channels based on established brain regions, get 
%      corresponding indices
% (7) Save the processed data to a mat file
%

clear all; 

debug = 0;
cfg.show = 1; % Show the time-series plots

% (1) Load the data
%
% cfg.data.dir = 'X:\Aquarids\code\loc2010';
% cfg.data.filename = 'hdm1.crx';
[cfg.data.filename, cfg.data.dir] = uigetfile('*.crx', 'Select Crx file'); %AG edit
% cfg.data.dir = 'C:\Users\j2m172\Desktop\21May2010\loc2010\NvsV_CRX';
% cfg.data.filename = 'N_R1.crx';
cfg.data.debug = debug;

cfg.save.do = 1;
cfg.save.dir=pwd; %AG edit
% cfg.save.dir = 'X:\Aquarids\code\loc2010';
[pathstr, name, ext] = fileparts( cfg.data.filename ); 
cfg.save.filename = [name '.lrx'];


% (2)Load the channel relabelling xml file & relabel
%
cfg.relabel.do = 1;
[cfg.relabel.filename, cfg.relabel.dir] = uigetfile('*.xml', 'Select xml relabel file for crx'); %AG edit
% cfg.relabel.dir = 'C:\Users\j2m172\Desktop\21May2010\loc2010\AllRelabelForHM';
% cfg.relabel.filename='Bardin~ Jon_bc6edf96-a434-4e81-92bb-88affad40960.xml';
cfg.relabel.debug = debug;


% (3) Remove non-neurological channels, eg EKG, EMG, EOG, etc 
% (4) Remove any bad channels, eg with drift or noticeable artifacts etc 
%
%AG - I don't need to do this since exported data only has channels I want
% cfg.rmChannel.nonneuro = { 'EMGREF', 'EMG', 'EMGp', 'EMGd' 'EKG', 'ECG', 'EOG' }; 
% cfg.rmChannel.other = { 'REF' };
% %cfg.rmChannel.other = { };
% cfg.rmChannel.debug = debug;


% (5) Add on the reference channel (and get it's position too)
% Set to empty cell, {}, or to {'AVG'} for average referencing
%
 cfg.addReference.channel = {'AVG'}; %AG edit
% cfg.addReference.channel = {'FCz'};


% (6a) Match the current data labels with the 10-20 standard labels, get the
%     electrode location
%
cfg.stdlabel.dir='/home/chronux/cornell/data/lbex_data/loc2010';
% cfg.stdlabel.dir = 'C:\Users\j2m172\Desktop\21May2010\loc2010';
% Sphere2.dat: Units are meters, coordinates in Cartesian 
cfg.stdlabel.filename = 'electrodePositions_Sphere2.dat' ;
cfg.stdlabel.debug = debug;


% (7) Save the processed data to a mat file
%
% Process & Save to file
%
ts=crxProcess( cfg ); 


