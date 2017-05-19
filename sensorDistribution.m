function sensorDistribution( varargin )
% Display the EEG/MEG Sensors
% Sensor locations   
% For top of coil use: coil( 1 ) or coilHC( 1 )
% For bottom of coil use: coil( 2 ) or coilHC( 2 )
% For coil positions in Dewar coordinates use: coil
% For coil positions in Head coordinates use: coilHC

if ~nargin || isempty( varargin ), error( 'sensorDistribution:: Data not specified' ); end

location = varargin{1};
x = location(:,1);
y = location(:,2);
z = location(:,3);
s = 9 * ones( size(x) );
%c = [ 0.0 * ones( size(x) ), 0.65 * ones( size(x) ), 0.0 * ones( size(x) ) ];

% Following designed to set centerline channels to black.
%
c=zeros(length(x),3); % first set all black
c(x>0, 1) = 1; % red
c(x<0, 3) = 1;  c(x<0, 2) = 0.75;  c(x<0, 1) = 0.5;  % bluish

if nargin == 1
    figure;
    scatter3( x, y, z, s, c, 'filled' );
    %scatter3( x, y, z, s, z, 'filled' );
    grid off; %title( '3D Sensor Distribution' );

    figure; scatter3( x, y, z, s, c, 'filled' ); view( 180, 90 ); grid off; title( 'Sensor Distribution: Top View' );
    figure; scatter3( x, y, z, s, c, 'filled' ); view( 180, 0 ); grid off; title( 'Sensor Distribution: Side View' );

elseif nargin == 2
    vxCentroid = varargin{2};
    
    figure;
    scatter3( x, y, z, s, c, 'filled' );
    %scatter3( x, y, z, s, z, 'filled' );
    grid off; %title( '3D Sensor Distribution' );
    hold on;
    
    x=vxCentroid(:,1);
    s = 8 * ones( size(x) );
    c=zeros(length(x),3); % first set all black
    c(x>0, 1) = 1; % red
    c(x<0, 3) = 1;  c(x<0, 2) = 0.75;  c(x<0, 1) = 0.5;  % bluish
    
    scatter3( vxCentroid(:,1), vxCentroid(:,2), vxCentroid(:,3), s, c );
    
end
