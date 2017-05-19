% Script for analytic preprocessing of time-series data
% Operations: detrend, band-pass, line-noise removal, and re-referencing 

%% Load from & Save to ...
%
[cfg.load.filename, cfg.load.dir] = uigetfile('*.lrx', 'Select lrx file'); %AG edit
% cfg.load.dir = 'C:\Users\j2m172\Desktop\21May2010\loc2010';
% % cfg.load.filename = 'hdm1.lrx';
% cfg.load.filename = 'N_R1.lrx';
cfg.save.dir=pwd; %AG edit
% cfg.save.dir = cfg.load.dir ;
[pathstr, name, ext, versn] = fileparts( cfg.load.filename ); 
cfg.save.filename = [name '.prx'];
cfg.show = 1; % Show the time-series plots

%% Load the crx preprocessed data into matlab session
%
load([ cfg.load.dir filesep cfg.load.filename ], '-mat')

% Extract basic parameters
frequency = ts.frequency; 
trials = ts.trials;

%% Set operations to be performed
%
ops.detrend = 1;
ops.filter = 0; % Don't use
ops.lineRemoval = 0; % Don't use
ops.referencing = 1; 
params.operations = ops; 

%% Set referencing parameters
%
r.scheme = 'average'; % 'none', 'average' or 'channel' for a specific channel
params.referencing = r;

%% Set data parameters
%
d.samplingFrequency = frequency;    % In Hertz
params.data = d; 

%% Set detrend parameters
% 
dt.type = 'linear'; % constant, linear
params.detrend = dt; 

%% Set Filter parameters
%

f.band = [ 2.5, 75 ];     % In Hertz, bandPass(2) > bandPass(1)
f.type = 'Rc-BW';
f.dpsspro.parms = []; % struct, not specified yet but plan for it...
params.filter = f; 

%% Set line noise removal paramters
% 

lr.type = 'bwNotch'; % 'bwNotch' or 'thomson'
lr.band = [59, 61]; % if notch remove all, if thomson search within band
lr.bwnotch.parms.order = 3;
lr.thomson.parms = []; % struct, not specified yet but plan for it...
params.lineRemoval = lr;

%% Process  
%
for i = 1:ts.trials
    [ ts.data{i}, params ] = preProcessTS( ts.data{i}, params );
end

%% Re-reference data 
%

if ops.referencing
    reference = ts.channel.reference;
    if strcmpi( r.scheme, 'average')
        for i=1:ts.trials
            ts.data{i} = ts.data{i} - repmat( mean(ts.data{i},2), 1, size(ts.data{i},2) );
        end
    elseif strcmpi( r.scheme, 'none')

    else
        error( [mfilename ' Re-referencing not implemented, currently only average reference possible.'] );
    end
end


%% Save to file
%
ts.date = datestr(clock);
ts.preprocess = params;
save([ cfg.save.dir filesep cfg.save.filename ], 'ts');


%% Display basic time-series & sensor distribution
%
if cfg.show
    for i=1:ts.trials
        figure;
        plot( [1:size(ts.data{i},1)]/ts.frequency, ts.data{i}(:,:) );
        title( ['Snippet ' num2str(i) ' of ' num2str(ts.trials)] );
    end
end


