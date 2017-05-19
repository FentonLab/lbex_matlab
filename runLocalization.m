clear all;

% Begin localization of data. This script needs data from runConc,
% runfTransform
%  

% Load concentration EVP from ... 
% This is NOT a mat-file and file/dir needs to be passed on localizePower
% directly. 
%
[cevCfg.filename, cevCfg.dir] = uigetfile('*.cev', 'Select cev file'); %AG edit
% cevCfg.dir = 'C:\Documents and Settings\Andrew Goldfine\My
% Documents\EEGResearch\LBEX\loc2010'; %AG edit though likely wrong
% % cevCfg.filename = 'hdm1.cev';
% cevCfg.filename = 'N_R1.cev';

% Load Fourier transformed data from ... 
% This IS a mat-file and file/dir needs to be passed on localizePower
% directly. 
%
%AG - this seems to want the data that I create so will change below
[ftCfg.load.filename, ftCfg.load.dir] = uigetfile('*.frx', 'Select frx file'); %AG edit
% ftCfg.load.dir = 'C:\Users\j2m172\Desktop\21May2010\loc2010';
% % ftCfg.load.filename = 'hdm1.frx';
% ftCfg.load.filename = 'N_R1.frx';
lzCfg.save.dir=pwd; %AG edit
% lzCfg.save.dir = ftCfg.load.dir ;
[pathstr, name, ext] = fileparts( ftCfg.load.filename ); 
lzCfg.save.filename = [name '.lzp']; % lzp = localize power, lzt = localize time-series
lzCfg.show = 1; % Show plots & messages

% Set the concentration EVP parameters
%

% Eigenspectrum cutoff for concentration EVP
% This regularizes the EVP, but is NOT the concentration eigenvalue
% retention criteria. Specify from [0,1]
%
cevCfg.cutOff = 0.90; 

% Criteria for number of concentrated spatial tapers retained
% This is applied for each ROI considered. 
%
cevCfg.threshold = 0.75; % consider voxel if lambda(1)>threshold

% Load the FFT processed data into matlab session
%
load([ ftCfg.load.dir filesep ftCfg.load.filename ], 'ts', '-mat')

% Extract basic parameters
%f = ts.TF.freqRange{i}; indf = find( f > fCutOff(1) & f < fCutOff(2)); f = f(indf);

% Localize power
%
nTrials = length(ts.TF.data);
nTrials = 1;
for i = 1:1
% for i = 1:nTrials
    [ lz.power{i} ] = localizePower( ts.TF.data{i}, cevCfg );
end

% Save to file
%
lz.date = datestr(clock);
lz.cfg = lzCfg;
save([ lzCfg.save.dir filesep lzCfg.save.filename ], 'lz');


if lzCfg.show
    
    fCutOff = [ 0 100 ];
    speedOfDisplay = 0.2; % seconds before next plot appears
    baseSubtract = 0; % Use average spectrum as baseline
    
    % Fix the time axis limits since trial lengths may vary
    %
    tmax = 0; for i=1:nTrials, tmax = max( tmax, max(ts.TF.timeRange{i}) ); end
    
    % Max/min power
    %
    pmax = realmin; for i=1:nTrials, pmax = max( pmax, max(lz.power{i}(:)) ); end
    pmin = realmax; for i=1:nTrials, pmin = min( pmin, min(lz.power{i}(:)) ); end

    % Use this value to set a fixed colormap
    % Allows for a consistent map for visual comparison
    %
    clims =[log10(pmin), log10(pmax)];
    
    figure; pwFig = axes('Parent', gcf); set(gca,'YDir','reverse');
    
    for i = 1:nTrials
        t = ts.TF.timeRange{i};
        f = ts.TF.freqRange{i}; indf = find( f > fCutOff(1) & f < fCutOff(2)); f = f(indf);
            
        for c = 1:size(lz.power{i},3)
            q = squeeze( lz.power{i}(indf,:, c) ); % tapers x frequency x windows/time
            imagesc( [0,tmax], fCutOff, log10( q ),  'Parent', pwFig );
            set(gca, 'CLim', clims); % Fix the color scale across trials & channels
            axis xy; colorbar;
            title( ['Trial: ' num2str(i) ', Voxel #' num2str(c) ] );
            pause(speedOfDisplay);
        end
    end
end

