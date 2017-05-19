function [ leadfield ] = dipoleInSphere3d3( sensorLocations, dipoleLocations, doLocationCheck )
%
% Chronux:: LBEX - A Source Localization Toolbox 
%
% Purpose: Evaluate the leadfield matrix. 
% Model: Dipole within a sphere 
%
% Inputs:
% 
%   sensorLocations = The Cartesian coordinates of the sensor locations.
%   sensorLocations( m, k ); m = 1, ..., nSensors; k = 1,2,3 are the
%   Cartesian components.
%
%   dipoleLocations = The Cartesian coordinates of the dipole locations.
%   dipoleLocations( m, k ); m = 1, ..., nDipoles; k = 1,2,3 are the
%   Cartesian components. 
%
%   doLocationCheck = 0=No or 1=Yes
%                     Check that the dipole & sensor locations are within
%                     non-intersecting spheres. A rough check... 
%
%
% Outputs:
%
%   leadfields = Matrix of leadfield calculations. leadfields( i, j ); 
%   Storage scheme is: 
%   m = 1, ..., nDipole; n = 1, ..., nSensors; 
%   i = p + 3*(m-1); j = q + 3*(n-1);
%   p,q = 1,2,3.
%
% Note:
%   1. The { 1+3*(n-1), 2+3*(n-1), 3+3*(n-1) } columns of 
%      the leadfield corresponds to the magnetic field sensor distribution 
%      due to an orthogonal dipole (or other primary singularity) 
%      oriented along = {1,2,3}, at the nth voxel. For a given column, the rows are
%      arranged as: the sequential rows { 1+3*(m-1), 2+3*(m-1), 3+3*(m-1) }
%      are the 3 vector components of the magnetic field at the mth sensor.
%      
%   2. For a given sensor-dipole pair, labelled by (m,n), the leadfield is
%      a 3x3 matrix. Therefore the storage scheme are 3x3 such blocks, with a total 
%      of m*n blocks. It might be easier to think of extracting these 3x3
%      blocks from "leadfield(i,j)": given (m,n) sensor-dipole label, the
%      3x3 matrix is given by,
%      leadfield( 1+3*(m-1) : 3+3*(m-1), 1+3*(n-1) : 3+3*(n-1) )

% Don't switch the sequence below !!
if nargin < 3,  doSensorCheck = 1;   end
checkInputs( sensorLocations, dipoleLocations, doLocationCheck  )
nSensors = size( sensorLocations, 1 );  nDipoles = size( dipoleLocations, 1 );

rsV = zeros( nDipoles, 3 ); rs = zeros( nDipoles, 1 );
aV = zeros( nDipoles, 3 ); a2 = zeros( nDipoles, 1 ); a = zeros( nDipoles, 1 ); 
F = zeros( nDipoles, 1 ); FInv = zeros( nDipoles, 1 ); gradF = zeros( nDipoles, 3 ); 
pFinv = zeros(nDipoles, 3); pCdDS = zeros( nDipoles, 1 );
leadfield = zeros( 3*nDipoles, 3*nSensors );

rdV = dipoleLocations;    % nDipoles x 3
rd = sqrt( sum( rdV.*rdV, 2 ) );  % nDipoles x 1
rsArr = sqrt( sum( sensorLocations.*sensorLocations, 2 ) );  % nSensors x 1

p100crRdV = [ zeros(nDipoles,1), -rdV(:,3), rdV(:,2) ]; % p x rdV when p=[1,0,0]
p010crRdV = [ rdV(:,3), zeros(nDipoles,1), -rdV(:,1) ]; % p x rdV when p=[0,1,0]
p001crRdV = [ -rdV(:,2), rdV(:,1), zeros(nDipoles,1) ]; % p x rdV when p=[0,0,1]

for ns = 1 : nSensors
    rsV = repmat( sensorLocations( ns, : ), nDipoles, 1 ); % nDipoles x 3
    rs = repmat( rsArr(ns), nDipoles, 1 );    % nDipoles x 1
    
    aV = rsV-rdV; % nDipoles x 3
    a2 = sum( aV.*aV, 2 ); % nDipoles x 1
    a = sqrt( a2 ); % nDipoles x 1
    
    F = a .* ( rs .* a + sum( rsV .* aV, 2) ); % nDipoles x 1
    FInv = 1./F; % nDipoles x 1
    gradF = (aV .* repmat( ( a + rs + (F ./ a2)), 1, 3 ) + rsV .* repmat( ( a + (a2./rs) ), 1, 3 ) ); % nDipoles x 3
    
    pFinv = p100crRdV .* [ FInv FInv FInv ]; % nDipoles x 3
    pCdDS = sum( pFinv.*rsV, 2 ) .* FInv; 
    leadfield( [1:3:end], 1 + 3*(ns-1) : 3 + 3*(ns-1) ) = ( pFinv - [pCdDS pCdDS pCdDS] .* gradF );
    
    pFinv = p010crRdV .* [ FInv FInv FInv ]; % nDipoles x 3
    pCdDS = sum( pFinv.*rsV, 2 ) .* FInv; % nDipoles x 1
    leadfield( [2:3:end], 1 + 3*(ns-1) : 3 + 3*(ns-1) ) = ( pFinv - [pCdDS pCdDS pCdDS] .* gradF );
         
    pFinv = p001crRdV .* [ FInv FInv FInv ]; % nDipoles x 3
    pCdDS = sum( pFinv.*rsV, 2 ) .* FInv; 
    leadfield( [3:3:end], 1 + 3*(ns-1) : 3 + 3*(ns-1) ) = ( pFinv - [pCdDS pCdDS pCdDS] .* gradF );
end

% constant multiplier: mu_{0} / (4 * pi )
leadfield = 1e-7 * leadfield;	

if 0
    figure; imagesc( leadfield); colorbar('vert'); 
    title( ['dipoleInSphere3d:: (Vector) Leadfield Matrix'] );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to check inputs for dipoleInSphere3d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkInputs( sensorLocations, dipoleLocations, doLocationCheck )

if size( sensorLocations, 2 ) ~= 3, error( ' dipoleInSphere3d:: 3 Cartesian components of sensor locations required, ( x, y, z ) ' ); end
if size( dipoleLocations, 2 ) ~= 3, error( ' dipoleInSphere3d:: 3 Cartesian components of dipole locations required, ( x, y, z ) ' ); end

% This is a rough check...don't stop execution...
if doLocationCheck
    tolerance = 1.e-6;
    rs = sum( sensorLocations.^2, 2 );
    rd = sum( dipoleLocations.^2, 2 );
    rsMin = min( rs );
    rdMax = max( rd );
    
    % Check whether sensor locations are outside source space
    indx = find( ( rs - rdMax ) <= tolerance );
    if ~isempty( indx )
        disp( 'dipoleInSphere3d:: Following Sensors within Source Space Sphere... ' );
        disp( ['Sensors Indices:: ', num2str( indx' )] );
        %error( ' ' ); 
    end ;
    
    % Check whether dipole locations are inside sphere
    indx = find( ( rd - rsMin ) >= tolerance );
    if ~isempty( indx )
        disp( 'dipoleInSphere3d:: Following Dipoles outside Source Space Sphere... ' );
        disp( ['Dipoles Indices:: ', num2str( indx' )] );
        %error( ' ' ); 
    end ;
    
end

