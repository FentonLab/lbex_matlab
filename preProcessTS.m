function [ data, params ] = preProcessTS( data, params )
%   Routine to preprocess a time-series. In particular, it allows for the
%   following preprocessing to be conducted: de-trend the time-series,
%   filter the time-series, and finally, remove the line-noise from the
%   time-series. Each of these pre-processing steps is optional.
%
%  Input: 
%   data(t,c): t is the time or sample index; c the channel/trial index
%
%   params has the following fields: 
%   params.{operations, detrend, filter, lineRemoval}
%   params.data
%
%   Desription of the params fields: 
%
%  params.operations:
%   params.operations further has the fields: detrend, filter, lineRemoval
%
%       params.operations.{detrend,filter,lineRemoval} = 0,1
%
%   depending on whether the data should detrended, filtered & line-noise
%   removed. 0 = no, 1=yes.
%
%
%  params.detrend:
%   params.detrend holds the detrend controlling parameters.
%
%       params.detrend.type = 'constant' or 'linear'
%
%   see matlab 'help detrend' for what these parameters do.
%
%
%  params.filter:
%   params.filter holds the filter choice & controlling parameters.
%
%       params.filter.type = 'rc-bw' or 'dpsspro'
%       params.filter.band() = [ lower_frequency, upper_frequency ]
%
%   Above, 'rc-bw' is an single-stage RC high-pass filter combined with a
%   2nd-order Butterworth low-pass filter; 'dpss-pro' is a dpss-based
%   projection filter. band(1)= lower_frequency, is the high-pass cutoff 
%   frequency. band(2) = upper_frequency is the low-pass cutoff frequency. 
%   Both these cutoff frequencies should be un-normalized and in Hertz(Hz). 
%   parms() is an array holding filter specific parameters (to be
%   implemented yet.) In addition to the above parameters, to allow for 
%   dpss projection filter parameters (when 'dpsspro' selected), the 
%   following fields must be specified: 
%
%       params.filter.dpsspro.parms =  struct, not specified yet but plan for it...
%
%
%  params.lineRemoval: 
%   params.lineRemoval holds the line-removal choice & parameters.
%
%       params.lineRemoval.type = 'notch' or 'thomson'
%       params.lineRemoval.band() = [ lower_frequency, upper_frequency ]
%
%   Above, 'notch' is the conventional notch-filter type approach to remove
%   the line-noise around 50Hz (Europe) or 60Hz (US). 'thomson' is a
%   sophisticated technique for line-noise removal (see Walden). band() is
%   the un-normalized frequency band (in Hz) within which the above 
%   algorithms are applied. For notch operations, it will remove all
%   frequency components with the specified band; for thomson operation, it
%   will use the band() to detect the precise line-component and then
%   extract the contribution due to line-noise. band() is un-normalized &
%   must be specified in Hertz(Hz). In addition to the above parameters, 
%   to allow for specific parameters (eg, when 'thomson' selected), the 
%   following fields must be specified: 
%
%       params.lineRemoval.thomson.parms =  struct, not specified yet but plan for it...
%       params.lineRemoval.bwnotch.parms.order = 6;
%
%
%  params.data:
%   params.data will hold acquisition specific metadata. Currently, only
%   one field is required: 
%
%       params.data.samplingFrequency = sampling rate/frequency in Hz. 
%
%
%  Outputs: 
%
%   data(t,c) = same as input, but containing the pre-processed data. 
%
%   params: Same as input, but may be supplemented by default values. 
%   If user doesn't provide details, some defaults values are set. Hence
%   modified params returned in call.

%%

rName = 'preProcessTS';
try
    % Check inputs
    %
    params = checkInputs( data, params );
    if ~isempty(lastwarn), [wmsg, wmsgId] = lastwarn; disp( [rName '::' wmsgId] ); disp( wmsg );end
    
    ops = params.operations;
    
    % De-Trend
    %
    if (ops.detrend == 1)
        dt = params.detrend;
        switch lower( dt.type )
            case 'constant'
                data = detrend( data, 'constant' );
            case 'linear'
                data = detrend( data );
        end
    end
    
    % Filter
    %
    if (ops.filter == 1)
        d = params.data;
        f = params.filter;
        switch lower( f.type )
            case {'rc-bw'}
                % Hi-Pass the data
                data = rcHighPass( data, f.band(1), d.samplingFrequency );
                % Low-Pass the data
                data = bwLowPass( data, 2, f.band(2), d.samplingFrequency );
            case 'dpsspro'
                disp('dpssPro not implemented...')
                %data = dpssPro( data, params );
        end
    end
    
    % Line-Removal
    %
    if (ops.lineRemoval == 1)
        d = params.data;
        lr = params.lineRemoval;
        switch lower( lr.type )
            case 'bwnotch'  % Butterworth bandstop 
                data = bwBandStop( data, lr.bwnotch.parms.order, ...
                                   lr.band, d.samplingFrequency );
            case 'thomson'
                disp('thomson line-removal not implemented...')
        end
    end

catch 
    err = lasterror; disp( [rName '::' err.identifier] ); disp( err.message );
    return;
end

function [ params ] = checkInputs( data, params )
%   Check inputs confirm to specifications in calling routines. Also set
%   some default values if not provided by the user. 
%

%% Data tolerances
% Look into robust setting at a later day
%
fTolerance = 0.25; % bandpass length
lrTolerance = 0.01; % line-removal band

%% Warnings & Error messages
%
msg = []; 
wmsg = [];
% Keep the next two or else will get previous warns/errors on new runs
lastwarn(''); 
% lasterror('reset');

%% Check params.operations fields
%
if ~isfield( params, 'operations' ) || isempty( params.operations )
    wmsg = [ wmsg 'No params.operations fields set. No preprocessing done!\n'];
    ops.detrend = 0; ops.filter = 0; ops.lineRemoval = 0; params.operations = ops;
else
    ops = params.operations;
    if ~isfield( ops, 'detrend' ) || isempty( ops.detrend ), ops.detrend = 0; end
    if ~isfield( ops, 'filter' ) || isempty( ops.filter ), ops.filter = 0; end
    if ~isfield( ops, 'lineRemoval' ) || isempty( ops.lineRemoval ), ops.lineRemoval = 0; end
    params.operations = ops;
end

%% Check params.data fields
%
if (ops.filter == 1) || (ops.lineRemoval == 1) % Only two cases where data.sampfreq reqd
    if ~isfield( params, 'data' ) || isempty( params.data )
        msg = [ msg 'Field params.data missing!\n'];
    else
        d = params.data;
        % d.samplingFrequency
        %
        if ~isfield( d, 'samplingFrequency' )
            msg = [ msg 'Field params.data.samplingFrequency missing!\n'];
        else
            if isempty( d.samplingFrequency ) || ~isnumeric( d.samplingFrequency )
                msg = [ msg 'params.data.samplingFrequency is required and must be numeric.\n'];
            end
        end
    end
end

%% Check params.detrend fields
%
if (ops.detrend == 1) && isfield( params, 'detrend' )
    dt = params.detrend;
    % dt.type = 'linear' or 'constant'
    % 
    if ~isfield( dt, 'type' )
        msg = [ msg 'Field params.detrend.type missing!\n'];
    else
        if isempty(dt.type) || ~( strcmpi(dt.type, 'constant') || strcmpi(dt.type, 'linear') )
            msg = [ msg 'Field params.detrend.type empty or incorrect!\n'];
        end
    end
end

%% Minimum data length required
%
if size( data, 1 ) == 2, msg = [ msg 'Data length must be greater than 2. \n']; end

%% Check params.filter fields
%
if (ops.filter == 1) && isfield( params, 'filter' )
    f = params.filter;
    % f.type = 'rc-bw' or 'dpsspro'
    % 
    if ~isfield( f, 'type' )
        wmsg = [ wmsg 'Field params.filter.type missing! Assuming default to RC-high pass & BW-low pass.\n'];
        f.type = 'RC-BW';
    else
        if isempty(f.type) || ~( strcmpi(f.type, 'rc-bw') || strcmpi(f.type, 'dpsspro') )
            msg = [ msg 'Field params.filter.type empty or incorrect!\n'];
        end
    end
    % f.band: f.band(2) > f.band(1)
    %
    if ~isfield( f, 'band' )
        msg = [ msg 'Field params.filter.band missing!\n'];
    else
        if isempty( f.band ) || ~isnumeric( f.band )
            msg = [ msg 'params.filter.band is required and must be real.\n'];
        else
            if ( f.band(2) - f.band(1) < fTolerance )
                msg = [msg 'Field params.filter.band: band(2) - band(1) > ' num2str(fTolerance) '.\n'];
            end
        end
    end
    % Just in case any of the defaults were set here, transfer
    %
    params.filter = f;
end

%% Check params.lineRemoval fields
%
if (ops.lineRemoval == 1) && isfield( params, 'lineRemoval' )
    lr = params.lineRemoval;
    % lr.type = 'bwnotch' or 'thomson'
    % 
    if ~isfield( lr, 'type' )
        msg = [ msg 'Field params.lineRemoval.type missing!\n'];
    else
        if isempty(lr.type) || ~( strcmpi(lr.type, 'bwnotch') || strcmpi(lr.type, 'thomson') )
            msg = [ msg 'Field params.lineRemoval.type empty or incorrect!\n'];
        end
    end
    % lr.band: lr.band(2) > lr.band(1)
    %
    if ~isfield( lr, 'band' )
        msg = [ msg 'Field params.lineRemoval.band missing!\n'];
    else
        if isempty( lr.band ) || ~isnumeric( lr.band )
            msg = [ msg 'params.lineRemoval.band is required and must be real.\n'];
        else
            if ( lr.band(2) - lr.band(1) < lrTolerance )
                msg = [msg 'Field params.lineRemoval.band: band(2) - band(1) > ' num2str(lrTolerance) '.\n'];
            end
        end
    end
    
    % If bwNotch check parameters...
    %
    if strcmpi(lr.type, 'bwnotch')
        if ~isfield( lr, 'bwnotch' )
            msg = [ msg 'Field params.lineRemoval.bwnotch missing!\n'];
        else
            if ~isfield( lr.bwnotch, 'parms' )
                msg = [ msg 'params.lineRemoval.bwnotch.parms is required.\n'];
            else
                if ~isfield( lr.bwnotch.parms, 'order' ) || ...
                    isempty( lr.bwnotch.parms.order ) || ...
                   ~isnumeric( lr.bwnotch.parms.order )
                    msg = [ msg 'params.lineRemoval.bwnotch.parms.order is required and must be real.\n'];
                end
            end
        end
    end
            
    % If thomson check parameters...
    %
    if strcmpi(lr.type, 'thomson')
        if ~isfield( lr, 'thomson' )
            msg = [ msg 'Field params.lineRemoval.thomson missing!\n'];
        else
            % Don't provide extensive check here since chronux routine will
            % check it. 
            %
            if ~isfield( lr.thomson, 'parms' ) || isempty( lr.thomson.parms )
                msg = [ msg 'params.lineRemoval.thomson.parms is required.\n'];
            end
        end
    end
end

%% Offload the warnings accumulated to be displayed by caller
if ~isempty(wmsg), warning off checkInputs:Warnings; warning('checkInputs:Warnings', wmsg); end

%% If errors detected throw exception
if ~isempty(msg), error('checkInputs:Errors', msg); rethrow(lasterror); end


function [ y ] = rcHighPass(x, highPassFrequency, samplingFrequency)

% Implements a single-stage RC high pass digital filter
%
% x(time,channels/trials) a matrix as time x channels/trials
% highPassFrequency = Frequency point for high pass (Hertz). Scalar. 
% samplingFrequency = Frequency at which data was sampled (Hertz). Scalar.
% 

 try
    % time-constant is the inverse of the cut-off frequency;
    % timeConst = 1 / (2*pi*cutOff) = resistance(R) x capacitance(C) = RC
    % dt = delta_t = 1 / samplingFrequency
    const = samplingFrequency / (samplingFrequency + 2*pi*highPassFrequency);

    % data length
    t = size( x, 1 );
    c = size( x, 2 );
    y = zeros( size( x ) );

    % opt = 1 is twice as fast as opt=1 for large arrays
    opt = 1 ;
    
    % Recursion, so have to loop through
    % Arranged channel-wise to save memory for large arrays
    %
    if opt == 1
        yt = zeros(t,1);
        xd = zeros(t-1,1);
        for j=1:c %channel loop
            yt(1) = x(1,j);
            xd = x(2:end,j) - x(1:end-1,j);
            for i=2:t %time loop
                yt(i) = const*( yt(i-1) +  xd(i-1) );
            end
            y(:,j) = yt;
        end
    end
    
    % Direct Matlab 
    % direct form II transposed implementation of the standard difference equation
    %
    if opt == 2
        a(1) = 1.0; a(2) = -const;
        b(1) = const; b(2) = -const;

        %y = filter( b, a, x ); 
        y = filtfilt( b, a, x ); % zero phase distortion
            
    end
    
catch
    err = lasterror; err.identifier = 'rcHighPass:Error';
%     throwAsCaller(err);
    rethrow(err);
 end


function [ y ] = bwLowPass( x, order, cutoff, samplingFrequency )
% Implements an nth-order Butterworth low pass digital filter

% x(time,channels/trials) input data matrix as time x channels/trials
% order = Order of the Butterworth filter
% cutoff = Low-pass cut-off frequency in Hertz (Hz). 
% samplingFrequency = data sampling rate in Hertz (Hz).
%
% y(time,channels/trials) output filtered data as time x channels/trials
% 

try
    fc = 2*cutoff/samplingFrequency; % Normalized cut-off frequency
    y = zeros( size( x ) );
    
    [b,a] = butter( order, fc, 'low' );
    %y = filter( b, a, x ); 
    y = filtfilt( b, a, x ); % zero phase distortion
    
catch
    err = lasterror; err.identifier = 'bwLowPass:Error';
%     throwAsCaller(err);
    rethrow(err);
end


function [ y ] = bwBandStop( x, order, band, samplingFrequency )
% Implements an nth-order Butterworth band-stop digital filter

% x(time,channels/trials) input data matrix as time x channels/trials
% order = Order of the Butterworth filter
% band = [b1,b2] = bandstop frequency in Hertz (Hz). 
% samplingFrequency = data sampling rate in Hertz (Hz).
%
% y(time,channels/trials) output filtered data as time x channels/trials
% 

try
    fc = 2*band/samplingFrequency; % Normalized cut-off frequency
    y = zeros( size( x ) );
    
    [b,a] = butter( order, fc, 'stop' );
    %y = filter( b, a, x ); 
    y = filtfilt( b, a, x ); % zero phase distortion
catch
    err = lasterror; err.identifier = 'bwBandStop:Error';
%     throwAsCaller(err);
    rethrow(err);
end



