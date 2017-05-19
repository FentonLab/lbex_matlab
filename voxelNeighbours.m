function [ roiColumns, roiVoxelIndices ] = voxelNeighbours( roiVolume, roiCenteredOnVoxel, voxelCentroids, showDetails, showPlot )

%   Purpose: Given a voxel index, find the neighbouring voxel indices encompassed 
%            by a sphere center on the voxel.
%
%   Inputs:
%
%   roiVolume = Volume of the ROI expressed as a percentage of the total
%               source space volume. The roiRadius is determimed using
%               this.
%
%   voxelCentroids(i,k) = Centroid of all the voxels that constitute the source
%               space. i = 1,..,nVoxels; k=1,2,3. Cartesian coordimnates assumed.
%
%   roiCenteredOnVoxel = Index of the voxel used for centering the ROI.
%
%   Ouptuts:
%
%   roiVoxelIndices = Indices of the ROI defining voxels.
%

checkInputs( roiVolume, roiCenteredOnVoxel, voxelCentroids );

nVoxels = size( voxelCentroids, 1 );

centerVoxel = voxelCentroids( roiCenteredOnVoxel, : );  % 1 x 3
nROI = min( nVoxels, ceil(roiVolume * nVoxels) ); % round but max out to nVoxels
[distance2, sortedI] = sort( sum( ( voxelCentroids - repmat( centerVoxel, nVoxels, 1 ) ).^2, 2 ) );
roiVoxelIndices = sortedI(1:nROI);

if isempty( roiVoxelIndices ), error( 'voxelNeighbours:: Could not determine the ROI ... ' ); end ;

if length( roiVoxelIndices ) == 1, disp( ['voxelNeighbours:: Single voxel ROI selected...'] ); end

% Get the columns corresponding to the ROIs
roiColumns = []; 
for r = 1 : length( roiVoxelIndices )
    roiColumns = [ roiColumns 1+3*( roiVoxelIndices(r)-1 ) 2+3*( roiVoxelIndices(r)-1 ) 3+3*( roiVoxelIndices(r)-1 ) ]; 
end

% Plot the ROI cloud
if showPlot
    figure; 
    plot3( ...
    voxelCentroids( roiVoxelIndices, 1 ), voxelCentroids( roiVoxelIndices, 2 ), voxelCentroids( roiVoxelIndices, 3 ), 'k^', ...
    voxelCentroids( :, 1 ), voxelCentroids( :, 2 ), voxelCentroids( :, 3 ), 'y.');
    title( ['voxelNeighbours:: ROI Centroid Locations - See blue cluster'] );
    axis equal; 
end

% if Display
if showDetails
    disp( ['   ROI Centered at voxel: ' num2str( roiCenteredOnVoxel )  ...
           '   # of voxels in ROI: ' num2str( length( roiVoxelIndices ) )  ] );
%           '   ROI % Volume: ' num2str( length( roiVoxelIndices )*voxelLengthScale^3/volume )
end



function [ ] = checkInputs( roiVolume, roiCenteredOnVoxel, voxelCentroids )

if roiVolume > 1, disp( ['voxelNeighbours:: Warning...roiVolume must be =< 1' ] ); end

if roiCenteredOnVoxel == 0 || roiCenteredOnVoxel > size( voxelCentroids, 1)
    error( ['voxelNeighbours:: Centering voxel index specified incorrectly...' num2str(roiCenteredOnVoxel) ] );
end


