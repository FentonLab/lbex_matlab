function [ sensorLocations ] = eegSensorLayoutCheck( sensorLocations, scalpRadius )
% 
%   If sensor does not lie on the scalp surface, project it on the scalp. 
%   This needs more thought... currently project all sensors onto scalp sphere.  
%
%   sensorLocations = The Cartesian coordinates of the sensor locations.
%   sensorLocations( m, k ); m = 1, ..., nSensors; k = 1,2,3 are the
%   Cartesian components.

tolerance = 1.0e-5;
rs = sqrt( sum( sensorLocations .* sensorLocations, 2 ) );  % nSensors x 1
outliers = find( abs( rs - scalpRadius ) > tolerance*scalpRadius );
if ~isempty( outliers )
    disp( [ 'eegSensorLayoutCheck:: ' ] );
    disp( [ '    Following sensors not on specified scalp radius:' ] );
    disp( [ '    ' num2str( outliers' ) ] );
    disp( [ '    Outlier sensors projected onto scalp... ' ] );
    sensorLocations(outliers, :) = scalpRadius * ( sensorLocations(outliers,:) ./ repmat( rs(outliers), 1, 3 ) );
end
         