function [ power, freqs ] = localizePower_old( data, conc )

fprintf( 1, 'localizePower:: \n' ); 
tic;
if ~isfield( data, 'dataFile' ), error( ['localizePower:: Input data.dataFile not specified...' ] ); end
if ~isfield( data, 'filename' ), error( ['localizePower:: Output data.filename not specified...' ] ); end

% Open fequency data 
fidF = 0; 
[ fidF, message ] = fopen( [ data.dir filesep data.filename ], 'r'); 
if fidF == -1, error( 'localizePower:: Error opening FT data file...', message ); end

nTrials = fread( fidF, 1, 'uint32' ); 
tmp = fread( fidF, 3, 'uint32' ); 
nFreqs = tmp(1); nTapers = tmp(2); nChannels = tmp(3);

fprintf( 1, '   Found: Freqs = %g, Tapers = %g, Channels = %g, Trials = %g \n', nFreqs, nTapers, nChannels, nTrials );
fprintf( 1, '   Read from file: %s', data.filename ); 

blockSize = nFreqs * nTapers * nChannels;
fData = complex( zeros( nFreqs, nTapers, nChannels ) );

% Open concentration data 
fidC = 0; 
[ fidC, message ] = fopen( [ conc.dir filesep conc.filename ], 'r'); 
if fidC == -1, error( 'localizePower:: Error opening concentration data file...', message ); end

numberOfVoxels  = fread( fidC, 1, 'uint32' );
roiVolume       = fread( fidC, 1, 'double' );
nChannels       = fread( fidC, 1, 'uint32' );
S               = fread( fidC, nChannels, 'double' );
uDim            = fread( fidC, 2, 'uint32' );  
U               = fread( fidC, uDim', 'double' );

cutIndx = length( find( S  > S( 1 ) * conc.cutOff ) );
cutUS = U( :, 1 : cutIndx ) * diag( 1 ./ S( 1 : cutIndx ) );

power = zeros( nFreqs, numberOfVoxels, nTrials ); 
rFT = zeros( nFreqs, nTapers, nChannels ); 
iFT = zeros( nFreqs, nTapers, nChannels ); 

for vx = numberOfVoxels : -1 : 1
    v               = fread( fidC, 1, 'uint32' );
    rcDim           = fread( fidC, 1, 'uint32' ); 
    roiColumns      = fread( fidC, rcDim, 'uint32' );
    rvDim           = fread( fidC, 1, 'uint32' ); 
    roiVoxelIndices = fread( fidC, rvDim, 'uint32' );
    cDim            = fread( fidC, 1, 'uint32' );  
    concEigenvalues = fread( fidC, cDim, 'double' );
    vpDim           = fread( fidC, 2, 'uint32' );   
    Vp              = fread( fidC, vpDim', 'double' );
    
    if concEigenvalues( 1 ) >= conc.threshold
        mnIndx = min( nChannels, length( roiColumns ) );
        
        mStar = 0; mStar = length( find( concEigenvalues >= conc.threshold ) );
        concVectors = cutUS * Vp( 1:cutIndx, 1:mnIndx ) ;
        %concVectors = cutUS * Vp( 1:cutIndx, 1:mStar ) ;

        if vx==numberOfVoxels, figure; plot( [ 1: length(concEigenvalues)], concEigenvalues, 'k.' ); title( ['Conc Values '] ); drawnow; end
        %if vx==numberOfVoxels, figure; imagesc( concVectors' * kernel(:,roiColumns) ); colorbar( 'vert' ); drawnow; end
        %if vx==numberOfVoxels, figure; imagesc( concVectors' * kernel(:,roiColumns) * kernel(:,roiColumns)' * concVectors ); colorbar( 'vert' ); drawnow; end
        %if vx==numberOfVoxels, figure; imagesc( concVectors' * kernel * kernel' * concVectors ); colorbar( 'vert' ); drawnow; end
        %if vx==numberOfVoxels, figure; imagesc( V( roiColumns, : )'*V( roiColumns, : ) ); colorbar('vert'); title( ['V^t * V'] ); drawnow; end
        
        % Estimate power at voxel
        for tr = 1 : nTrials
            trIndx = fread( fidF, 1, 'uint32' );
            
            % fData = frequency x tapers x channels 
            rFT = fread( fidF, blockSize, 'real*8' );
            iFT = fread( fidF, blockSize, 'real*8' );
                        
            %rFT = reshape( rFT, [ nFreqs, nTapers, nChannels ] );
            %iFT = reshape( iFT, [ nFreqs, nTapers, nChannels ] );
            %fData = complex( rFT, iFT );
            fData = reshape( complex( rFT, iFT ), [ nFreqs, nTapers, nChannels ] );
            
            for nf = 1 : nFreqs
                % taperedBCf =  mStar x nTapers = (mStar x nChannels) x (nChannels, nTapers)
                cTaperedBf = concVectors( :, 1:mStar )' * ( squeeze( fData( nf, :, : ) ).') ; 
                % Current power estimate from voxel only
                power( nf, vx, tr ) = sum( sum( repmat( concEigenvalues( 1:mStar ), 1, nTapers ) .* ( cTaperedBf .* conj( cTaperedBf ) ) ) ) / nTapers;
            end
        end
        freqs = fread( fidF, nFreqs, 'real*8' );
        frewind( fidF );
        tmp = fread( fidF, 4, 'uint32' ); %read & discard
    else
        disp( ['Truncation: ' num2str( 10*log10( 1 - concEigenvalues(1) ) ) '  \lamda_1: ' num2str( concEigenvalues(1) ) ] );
        power( :, vx, : ) = reshape( repmat( realmin, [ nFreqs, nTrials ] ), [ nFreqs, 1, nTrials] );
    end
    
end

status = fclose( fidF ); if status == -1, error( 'localizePower:: Error closing file ...' ); end
status = fclose( fidC ); if status == -1, error( 'localizePower:: Error closing file ...' ); end

fprintf( 1, '   Done...%gs\n', toc );
