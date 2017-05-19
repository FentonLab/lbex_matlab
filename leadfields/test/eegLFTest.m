clear all;

library = 'E:\hiren\lbex';
addpath( genpath( library ) );

r(1) = 1.0; r(2) = 0.9467; r(3) = 0.8667; r(4) = 0.84;
c(1) = 1.0; c(2) = 0.0125; c(3) = 3.0; c(4) = 5.0;
hdm.model.radii = 'user';
hdm.model.conductivity = 'user';
hdm.radius = 11*r;
hdm.conductivity = c;
hdm.displayMsg = 1;
params.eegHeadModel = hdm;

% Sensor distribution
Ts = 5; Ps = 20; radius = 1.00 * hdm.radius(1);
sTheta = [ 0, pi/3 ]; sPhi = [ 0, 2*pi ];

theta = [ sTheta(1) : ( sTheta(2) - sTheta(1) ) / ( Ts - 1 ) : sTheta(2) ]' ;
phi = [ sPhi(1) : ( sPhi(2) - sPhi(1) ) / ( Ps - 1 ) : sPhi(2) ]' ;
rad = repmat( radius, Ts*Ps, 1 );
csTheta = cos( theta ) ; snTheta = sin( theta ); csPhi = cos( phi ); snPhi = sin( phi );

x = reshape( snTheta * csPhi' , Ts*Ps, 1 );
y = reshape( snTheta * snPhi' , Ts*Ps, 1 );
z = reshape( repmat( csTheta, 1, Ps ), Ts*Ps, 1 );
sensor.Location = repmat( rad, 1, 3 ) .* [ x y z ];

% Dipole space
vxlSpace.gridLimits(:,1) = [-hdm.radius(1),hdm.radius(1)];
vxlSpace.gridLimits(:,2) = [-hdm.radius(1),hdm.radius(1)];
vxlSpace.gridLimits(:,3) = [-hdm.radius(1),hdm.radius(1)];
vxlSpace.radius = hdm.radius(1);
vxlSpace.origin = [0,0,0];
vxlSpace.units = 'cm';
vxlSpace.plot = 0; vxlSpace.display = 1;

vxlSpace.voxelSize = 7.5;
vxlSpace = voxelizeSphere( vxlSpace );


%[sen, dip]=checkDipoleInSphere3d( 11.0, 10.0, 0.0000001, [ 5 10 ], [ 20 40 ], [0 pi/3; 0 2*pi], [0 pi; 0 2*pi] );
tic, [lf, params] = eegDipoleInSphere3d3( sensor.Location, vxlSpace.grid(vxlSpace.inside, :), params ); toc

% Compare 
%
% To avoid ft's warning messages... project sensors
if 1
sensor.Location = eegSensorLayoutCheck( sensor.Location, params.eegHeadModel.radius(1) ); 
addpath( genpath( 'E:\hiren\modified_fieldtrip-20050522' ) );
% Fieldtrip arranges radius/conductivity from smallest to largest -
% important!!
vol.r = fliplr(params.eegHeadModel.radius); 
vol.c = fliplr(params.eegHeadModel.conductivity);
tic, lf2 = eeg_leadfield( vxlSpace.grid(vxlSpace.inside, :), sensor.Location, vol ); toc
lf2 = lf2';
end

f1 = lf ; f2 = lf2 ;
figure;
subplot(1,3,1), imagesc( f1 ); colorbar('vert'); 
    title( ['eegDipoleInSphere3d:: Leadfield Matrix'] );
subplot(1,3,2), imagesc( f2 );colorbar('vert'); 
    title( ['fieldtrip:: Leadfield Matrix'] );
subplot(1,3,3), imagesc( f1-f2 );colorbar('vert'); 
    title( ['Error'] );
