function [ leadfield ] = eegDipoleInHomoSphere( sensorLocations, dipoleLocations, scalpRadius, conductivity, displayMsg )
%
% Chronux:: LBEX - A Source Localization Toolbox 
%
% Purpose: Evaluate the eeg leadfield matrix. 
% Model: Dipole within a homogenous sphere 
%
% Caution: It is assumed that the sensor & dipole locations are input
% relative to a coordinate system located at the center of the sphere. 
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
%   displayMsg = 0/1.  0 = Don't display. 
%
%
%%%%%%%%%%
% Outputs:
%%%%%%%%%%
%
%   leadfields = Matrix of leadfield calculations. leadfields( i, j ); 
%   Storage scheme is: 
%   m = 1, ..., nDipole; n = 1, ..., nSensors; 
%   i = p + 3*(m-1); j = n;
%   p = 1,2,3.
%
% Note:
%   1. The { 1+3*(m-1), 2+3*(m-1), 3+3*(m-1) } rows of 
%      the leadfield corresponds to the scalp potential distribution 
%      (ie sensor distribution) due to an orthogonal dipole  
%      oriented along = Cartesian {1,2,3}, at the nth voxel. 
%      

rName = 'eegDipoleInHomoSphere:: ';
spc   = '    ';
if size( sensorLocations, 2 ) ~= 3, error( [spc '3 Cartesian components of sensor locations required, ( x, y, z ) '] ); end
if size( dipoleLocations, 2 ) ~= 3, error( [spc '3 Cartesian components of dipole locations required, ( x, y, z ) '] ); end

dm = displayMsg;

% constant multiplier: 1 / (4 * pi * conductivity(1) * R^2)
% In paper it's sigma(4), but reverse of our indexing. It's the scalp
% conductivity.
const = 1 / ( 4*pi*conductivity);	

% 
% Currently, simply project all the EEG sensors onto the scalp sphere. 
% To be thought through & changed...
%
sensorLocations = eegSensorLayoutCheck( sensorLocations, scalpRadius ); 

if dm, disp( rName ); end

% Assume coords reported relative to sphere center...issue to be fixed?
nSensors = size( sensorLocations, 1 );  nDipoles = size( dipoleLocations, 1 );
rd2 = sum( dipoleLocations .* dipoleLocations, 2 );  % nDipoles x 1
rd = sqrt( rd2 );  % nDipoles x 1
%rs2 = sum( sensorLocations .* sensorLocations, 2 );  % nSensors x 1
%rs = sqrt( rs2 );  % nSensors x 1

r0 = dipoleLocations ./ repmat( rd, 1, 3 ); % unit radial, nDipoles x 3
f = rd / scalpRadius; % nDipoles x 1

rsVec=zeros(nDipoles,3); rs2=zeros(nDipoles,1); rs=zeros(nDipoles,1); 
rp2=zeros(nDipoles,1); rp=zeros(nDipoles,1); tmp = zeros( nDipoles,3); 
rdCosPhi=zeros(nDipoles,1); 
leadfield = zeros( 3*size(dipoleLocations,1), nSensors ); % Don't change this!!

for ns = 1 : nSensors
    rsVec = repmat( sensorLocations(ns, :), nDipoles, 1 ); % nDipoles x 3 = repetition
    rs2 = sum(rsVec .* rsVec, 2); 
    rs = sqrt( rs2.*rs2 ); % nDipoles x 1
    rp2 = sum( (rsVec - dipoleLocations) .* (rsVec - dipoleLocations), 2 );  % nDipoles x 1
    rp = sqrt( rp2 ); 
    rdCosPhi = sum( dipoleLocations .* rsVec, 2 ) ./ rs; % nDipoles x 1
    tmp = ( rsVec.*repmat( rdCosPhi, 1, 3 ) - dipoleLocations.*repmat(rs, 1, 3) )...
        ./ repmat( rs + rp - rdCosPhi, 1, 3 ); % nDipoles x 3
    leadfield( :, ns ) = reshape( (2*(rsVec-dipoleLocations)./repmat( rp.*rp2, 1, 3 )...
        + ( (rsVec - tmp) ./ repmat( rs2.*rp, 1, 3) ) )' , 3*nDipoles, 1 );
end

leadfield = const*leadfield;

if 0
    figure; imagesc( leadfield ); colorbar('vert'); 
    title( ['eegDipoleInHomoSphere:: Leadfield Matrix'] );
end

