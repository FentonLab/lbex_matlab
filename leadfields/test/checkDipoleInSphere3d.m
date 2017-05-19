function [ sensor, dipole ] = checkDipoleInSphere3d( maxSensorSpaceRadius, interfaceRadius, minDipoleSpaceRadius, sensorDim, dipoleDim, sensorSpace, dipoleSpace )

% Routine to check dipoleInSphere3d.m by generating cloud of sensor &
% dipole locations
%
% Sensor & dipoles locations are generated randomly in the radial
% direction. The sensor space is radially bound by:
%   interfaceRadius + tolerance  < sesnsor Space < maxSensorSpaceRadius 
% and the source (dipole) space is radially bound by
%   minDipoleSpaceRadius < dipole Space < interfaceRadius + tolerance  
%
% Inputs:
% 
%   maxSensorSpaceRadius, minDipoleSpaceRadius = Self evident. Latter can
%   be zero.
%
%   interfaceRadius = radius that separates the source & dipole spaces
%
%   sensorDim = [ Ts, Ps ] = Discretizations in theta & phi directions. 
%               Total number of sensor locations = Ts x Ps
%
%   dipoleDim = [ Td, Pd ] = Similar to sensorDim for dipole locations
%
%   sensorSpace = [ Theta_1, Theta_2 ; Phi_1, Phi_2 ] = limits of the
%                 sensor space projected on a unit sphere. Default is
%                 entire sphere. Optional.
%
%   dipoleSpace = [ Theta_1, Theta_2 ; Phi_1, Phi_2 ] = limits of the
%                 dipole space projected on a unit sphere. Optional. 
%                 Default is entire sphere.
%
%   Example:
%
%   For sensor & dipole spaces covering the entire sphere: 
%   checkDipoleInSphere3d( 1.25, 1.0, 0.1, [ 5 10 ], [ 10 20 ] );
%
%   Partial coverage:
%   checkDipoleInSphere3d( 1.25, 1.0, 0.1, [ 5 10 ], [ 10 20 ], [ 0, pi/6; 0, 2*pi], [0,pi; 0,2*pi] );

if nargin < 6 
    sensorSpace = [ 0, pi ; 0, 2*pi]; dipoleSpace = sensorSpace;
end

doRandom = 0;

tolerance = 1.e-2;

radius( 1 ) = minDipoleSpaceRadius; radius( 2 ) = interfaceRadius; radius( 3 ) = maxSensorSpaceRadius;

%%%%%%%%%%%%% Sensor Space %%%%%%%%%%%%%%%%%

% s = sensor; d = dipole
Ts = sensorDim( 1 ); Ps = sensorDim( 2 );
sTheta( 1 ) = sensorSpace( 1, 1 ); sTheta( 2 ) = sensorSpace( 1, 2 ); 
sPhi( 1 ) = sensorSpace( 2, 1 ); sPhi( 2 ) = sensorSpace( 2, 2 ); 

if doRandom
    theta = sTheta( 1 ) + ( sTheta( 2 ) - sTheta( 1 ) ) * rand( Ts, 1 );
    phi = sPhi( 1 ) + ( sPhi( 2 ) - sPhi( 1 ) ) * rand( Ps, 1 );
    r = radius( 2 ) + ( radius( 3 ) - radius( 2 ) ) * rand( Ts*Ps, 1 );
else
    theta = [ sTheta(1) : ( sTheta(2) - sTheta(1) ) / ( Ts - 1 ) : sTheta(2) ]' ;
    phi = [ sPhi(1) : ( sPhi(2) - sPhi(1) ) / ( Ps - 1 ) : sPhi(2) ]' ;
    r = repmat( radius( 3 ), Ts*Ps, 1 );
end
csTheta = cos( theta ) ; snTheta = sin( theta ); csPhi = cos( phi ); snPhi = sin( phi );

x = reshape( snTheta * csPhi' , Ts*Ps, 1 );
y = reshape( snTheta * snPhi' , Ts*Ps, 1 );
z = reshape( repmat( csTheta, 1, Ps ), Ts*Ps, 1 );
sensor.Location = repmat( r, 1, 3 ) .* [ x y z ];

% plot sensor cloud
figure; plot3( sensor.Location( :, 1 ), sensor.Location( :, 2 ), sensor.Location( :, 3 ), 'r.'  );
title( ['Sensor Locations'] );
%axis equal; 

%%%%%%%%%%%%% Dipole Space %%%%%%%%%%%%%%%%%

% s = sensor; d = dipole
Td = dipoleDim( 1 ); Pd = dipoleDim( 2 );
dTheta( 1 ) = dipoleSpace( 1, 1 ); dTheta( 2 ) = dipoleSpace( 1, 2 ); 
dPhi( 1 ) = dipoleSpace( 2, 1 ); dPhi( 2 ) = dipoleSpace( 2, 2 ); 

if doRandom
    theta = dTheta( 1 ) + ( dTheta( 2 ) - dTheta( 1 ) ) * rand( Td, 1 );
    phi = dPhi( 1 ) + ( dPhi( 2 ) - dPhi( 1 ) ) * rand( Pd, 1 );
    r = radius( 1 ) + ( radius( 2 ) - radius( 1 ) ) * rand( Td*Pd, 1 );
else
    theta = [ dTheta(1) : ( dTheta(2) - dTheta(1) ) / ( Td - 1 ) : dTheta(2) ]' ;
    phi = [ dPhi(1) : ( dPhi(2) - dPhi(1) ) / ( Pd - 1 ) : dPhi(2) ]' ;
    r = repmat( radius( 1 ), Td*Pd, 1 );
end
csTheta = cos( theta ) ; snTheta = sin( theta ); csPhi = cos( phi ); snPhi = sin( phi );

x = reshape( snTheta * csPhi' , Td*Pd, 1 );
y = reshape( snTheta * snPhi' , Td*Pd, 1 );
z = reshape( repmat( csTheta, 1, Pd ), Td*Pd, 1 );
dipole.Location = repmat( r, 1, 3 ) .* [ x y z ];

% plot dipole cloud
figure; plot3( dipole.Location( :, 1 ), dipole.Location( :, 2 ), dipole.Location( :, 3 ), 'b.'  );
title( ['Dipole Locations'] );
%axis equal; 

% Combined cloud
figure; plot3( [sensor.Location( :, 1 ); dipole.Location( :, 1 )], [sensor.Location( :, 2 ); dipole.Location( :, 2 )], [sensor.Location( :, 3 ); dipole.Location( :, 3 )], 'g.'  );
title( ['Sensor & Dipole Locations'] );
%axis equal; 

