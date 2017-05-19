function [ cFig ] = getEEGHeadModel( cFig )
%
% This routine checks & helps configure the parameters required for an EEG
% 4-shell spherical model. Used for EEG leadfield calculations. The routine does the
% following: 
%
% (a) This routine will check the user's settings for the 4-shell EEG
% headmodel. In particular, it will check the radii & conductivity
% parameters. 
%
% (b) If the user specifies a head-model which isn't hard-coded here, it
% will then compute the c_n & a_n coefficients required for the EEG
% leadfield computations. 
%
% 
% Inputs: 
%
% A structure, cFig.*, with the following fields:
%
%   .radius(i), i=1,2,...,4. Optional. The radii of the regions delineating 
%   different conductivities. It must be ordered starting with the scalp radius,
%   radius(1), and descend towards the center. The radii should not be
%   normalized by the scalp radius (or any other). It is NOT necessary to
%   specify all the radii, but scalp radius, radius(1), must be specified.
%   If only the scalp radius is specified, a head model, *.model, must be
%   specified (see *.model below). Note, either, only scalp radius, or ALL
%   four radii must be specified; no partial specifications allowed. 
%
%   .conductivity(i), i=1,2...,4. i=1 a must, rest optional. The head 
%   model layer conductivities. They must be ordered starting with the
%   scalp conductivity, c(1), and descend towards the center. The 
%   conductivities should not be normalized by 
%   the scalp conductivity (or any other). It is NOT necessary to specify
%   all the conductivities, but scalp conductivities, c(1), must be specified.
%   If only the scalp conductivity is specified, a head model, *.model, must be
%   specified (see *.model below). Note, either, only scalp conductivity, or ALL
%   four conductivities must be specified; no partial specifications allowed. 
%
%   .model = 'rush-driscoll' or 'cuffin-cohen' or 'stok' or 'user'
%   If the .model field is not found, a user specified model is assumed 
%   internally.
%
%   .displayMsg = 0/1. Optional. Default = 0 = Don't display. 
%
% Notes: 
% 1. A 3-shell model can be simulated by merging 2 layers, ie specify 2 radii
% and conductivities the same. 
% 2. This routine will not check whether the sensors lie on the scalp
% radius. This issue is taken care of in the leadfield calling routine.
% 3. This routine calls anCoefficients.m. Parameters set in the calling
% program. Defaults set there too. Calling program is eegDipoleInSphere3d.m
%
%

rName = 'getEEGHeadModel:: '; 
spc   = '    ';
dm = 0; if isfield( cFig, 'displayMsg' ), dm = cFig.displayMsg; end
if dm, disp( rName ); end
defaultModel = 'rush-driscoll';

% Don't re-arrange IF-order below... will ruin logic...

%
% Check Radii
%
if ~isfield( cFig, 'radius' ) || isempty( cFig.radius )
    error( [rName 'At the very least, scalp (outer) radius, *.radius(1) must be provided.'] );
else
    nShells = length( cFig.radius );
    if nShells == 1
        if dm, disp( [spc 'Only 1 radius found, assuming outer scalp radius...'] ); end
        if ~isfield( cFig, 'model' ) || ~isfield( cFig.model, 'radii' ), cFig.model.radii = defaultModel; end
        if strcmp( cFig.model.radii, 'user' ), error( [rName 'Head model radii setting inconsistent with radii array input.']); end
    elseif nShells == 4
        cFig.model.radii = 'user';
        if dm, disp( [spc 'User specified 4-shell radii...'] ); end
        if ~isequal( cFig.radius, sort(cFig.radius, 'descend') )
            error( [rName 'Radii must be in descending order, r(1)>r(2)>r(3)>r(4). r(1)=scalp radius.' ] ); 
        end
    else
        error( [rName 'Must specify either outer (scalp) radius or 4-shell radii.'] )
    end;
end


%
% Check Conductivities
%
% if ~isfield( cFig, 'conductivity' ) || isempty( cFig.conductivity ) || length( cFig.conductivity ) < 4
%     if ~isfield( cFig, 'model' ) || ~isfield( cFig.model, 'conductivity' ), cFig.model.conductivity = defaultModel; end
%     if strcmp( cFig.model.conductivity, 'user' ), error( [rName 'Head model conductivity setting inconsistent with conductivity array input.']); end
% else    % Found a 4-shell conductivity profile
%     cFig.model.conductivity = 'user';
%     if dm, disp( [spc 'User specified conductivity profile...'] ); end
% end

% 31 Mar 2010
if ~isfield( cFig, 'conductivity' ) || isempty( cFig.conductivity )
    error( [rName 'At the very least, scalp (outer) conductivity, *.conductivity(1) must be provided.'] );
else
    nShells = length( cFig.conductivity );
    if nShells == 1
        if dm, disp( [spc 'Only 1 conductivity found, assuming outer scalp conductivity...'] ); end
        if ~isfield( cFig, 'model' ) || ~isfield( cFig.model, 'conductivity' ), cFig.model.conductivity = defaultModel; end
        if strcmp( cFig.model.conductivity, 'user' ), error( [rName 'Head model conductivity setting inconsistent with conductivity array input.']); end
    elseif nShells == 4
        cFig.model.conductivity = 'user';
        if dm, disp( [spc 'User specified 4-shell conductivity...'] ); end
    else
        error( [rName 'Must specify either outer (scalp) conductivity or 4-shell conductivity.'] )
    end;
end

iFlag = 0;
scalpRadius = cFig.radius( 1 ); % has to be present.
scalpConductivity = cFig.conductivity( 1 ); % has to be present.

cFig.model.radii = lower( deblank( cFig.model.radii ) );
cFig.model.conductivity = lower( deblank( cFig.model.conductivity ) );

% 
% Models
%
model{ 1 } = { 'rush', 'driscoll', 'rush-driscoll' };
model{ 2 } = { 'cuffin', 'cohen', 'cuffin-cohen' };
model{ 3 } = { 'stok' };
model{ 4 } = { 'user' };

%
% Set model radii
%
% Note: These are normalized r(*) values, Normalized by scalp radius.
% To get these dimensioned, need to multiply by scalp radius, which
% has a typical value of 0.10 meter.
%
switch cFig.model.radii
    case model{ 1 } % { 'rush', 'driscoll', 'rush-driscoll' };
        r(1) = 1.0; r(2) = 0.9467; r(3) = 0.8667; r(4) = 0.84;
        msg = 'rush-driscoll';
    case model{ 2 } %{ 'cuffin', 'cohen', 'cuffin-cohen' }
        r(1) = 1.0; r(2) = 0.95; r(3) = 0.9; r(4) = 0.75;
        msg = 'cuffin-cohen';
    case model{ 3 } % { 'stok' };
        r(1) = 1.0; r(2) = 0.9467; r(3) = 0.8667; r(4) = 0.84;
        msg = 'stok';
    case model{ 4 } % { 'user' };
        msg = 'user-specified';
        iFlag = 1;
        % If user provided radii, normalize them so that scalp radius = 1 
        r = cFig.radius / scalpRadius;
    otherwise
        error( [spc 'Specified head model for radii inappropriate...'] );        
end
if dm, disp( [spc 'EEG head model radii set to: ' msg] ); end
r=r(:);
cFig.radius = scalpRadius .* r;

%
% Set model conductivities
%
% Note: These are normalized c(*) values, Normalized by scalp conductivity.
% To get these dimensioned, need to multiply by scalp conductivity, which
% has a typical value of 0.33 (ohm-meter)^(-1).
%
switch cFig.model.conductivity
    case model{ 1 } % { 'rush', 'driscoll', 'rush-driscoll' };
        c(1) = 1.0; c(2) = 0.00125; c(3) = 3.0; c(4) = 1.0;
        a(1) = 2.292444243; a(2) = -0.169933547; a(3) = 0.014038780; a(4) = -0.000273030;
        msg = 'rush-driscoll';
    case model{ 2 } %{ 'cuffin', 'cohen', 'cuffin-cohen' }
        c(1) = 1.0; c(2) = 0.0125; c(3) = 3.0; c(4) = 1.0;
        a(1) = 3.719040299; a(2) = -0.278262604; a(3) = 0.013362279; a(4) = -0.000201385;
        msg = 'cuffin-cohen';
    case model{ 3 } % { 'stok' };
        c(1) = 1.0; c(2) = 0.0125; c(3) = 3.0; c(4) = 1.0;
        a(1) = 2.094413086; a(2) = -0.212973207; a(3) = 0.015140880; a(4) = -0.000311147;
        msg = 'stok';
    case model{ 4 } % { 'user' };
        msg = 'user-specified';
        iFlag = 1;
        % If user provided conductivity, normalize them so that scalp conductivity = 1 
        c = cFig.conductivity / scalpConductivity;
    otherwise
        error( [spc 'Specified head model for conductivity inappropriate...'] );
end
if dm, disp( [spc 'EEG head model conductivities set to: ' msg] ); end
c = c(:);
cFig.conductivity = scalpConductivity .* c;

%
% If either, radii or conductivities, are different from hard-coded models,
% recompute the coefficients, an().
%
if iFlag
    a = anCoefficients( cFig ); 
end

cFig.an = a(:);



        
