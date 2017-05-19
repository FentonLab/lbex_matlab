clear all;

% Create a time-frequency matrix from the data by running FFT on the data.
% This is NOT the same as spectrogram, since this does not evaluate the
% power - it's plain Fourier transform of the data to be used eventually
% for power localization. 
%
% To be called after a runcrx2lrx or runPreProcess run
%

% Load from & Save to ...
%
[cfg.load.filename, cfg.load.dir] = uigetfile('*.prx', 'Select prx file'); %AG edit
% cfg.load.dir = 'C:\Users\j2m172\Desktop\21May2010\loc2010';
% % cfg.load.filename = 'hdm1.prx';
% % cfg.load.filename = 'N_R1.lrx';
% cfg.load.filename = 'N_R1.prx';
cfg.save.dir=pwd; %AG edit
% cfg.save.dir = cfg.load.dir ;
[pathstr, name, ext, versn] = fileparts( cfg.load.filename ); 
cfg.save.filename = [name '.frx'];
cfg.show = 1; % Show the power spectrogram plots

% Load the crx preprocessed data into matlab session
%
load([ cfg.load.dir filesep cfg.load.filename ], '-mat')

% Extract basic parameters
Fs = ts.frequency; 
nTrials = ts.trials;

% Spectral Parameters
% 
% Specify tapers in the 3-form, [W T p] where W is the bandwidth (Hz), 
% T is the duration of the data, in secs, and p is an integer such that  
% 2TW-p tapers are used. Eg, if T=100 ms(millisecs)= 0.1s, W=5 Hz, 2TW=1
%
% To allow for further analysis in frequency domain ensure that the
% frequency grid is the same for all trial sets

W = 5; % Taper bandwidth in frequency domain. In Hertz, Hz, if Fs in Hz
% moving win: units = secs, consistent with T in tapers = [W T p]
movingwin = [400 25] / 1000; 
T = movingwin(1); % to be consistent - don't change!! 
p = 1;

tapers = [W T p];
% Don't remove this test
%
K  = floor(2*T*W - p); % as tested in getparams.m 
if ~K
    error( [mfilename, 'Number of tapers retained, K = 0']);
else
    disp( [mfilename, ' 2TW = ' num2str(2*T*W)]);
    disp( [mfilename, ' TW -p = ' num2str(2*T*W - p)]);
    disp( [mfilename, ' Number of tapers retained, K = ' num2str(K)]);
end

% Remaining parameters for spectral computations
% Note: 
% At this point fpass parameter inputs are best left unspecified so that
% they default to entire range, [0,Fs/2]. This since at this point we don't
% know which frequency band is of interest. However, it appears that some
% datasets are recorded at 1000 Hz, so best if bandlimited here...
%
params.pad = 0;
params.tapers = tapers;
params.Fs = Fs;
params.fpass = [0 100];

for i = 1:nTrials
    [Q,t,f] = mtTFgramc( ts.data{i}, movingwin, params );
    ts.TF.data{i} = Q;  % channels x tapers x frequency x windows/time
    % keep them trial dependent, hence unequal lengths
    ts.TF.freqRange{i} = f;
    ts.TF.timeRange{i} = t;    % same as # of window sliding steps
end

% Save to file
%
ts.TF.params = params;
ts.TF.cfg = cfg; % so that on reload it doesn't erase any existing cfg
ts.TF.date = datestr(clock);
save([ cfg.save.dir filesep cfg.save.filename ], 'ts');


