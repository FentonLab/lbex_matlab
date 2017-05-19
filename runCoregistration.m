clear all;

dataSet     = 'marmar_Niehoo1_20030903_01.ds';
sessionDir  = 'E:\nienke\NienkeScripts\specialIEEE_Calcs\';
dataCD      = 'E:\nienke\raw\cd_10a\' ;
mriFile = 'E:\nienke\fieldtripCompare\leadfields\leadfields\mri\martens_m.mri';


foi = 72.0848;
%foi = 70;
regularization = 0;


undscr = find( dataSet == '_' ); bckslsh = find( dataSet == filesep );
if isempty( bckslsh ), bckslsh = 0; end
postFix = [ dataSet( bckslsh+1 : undscr(1)-1 ) '_' ];
timeData  = [ sessionDir filesep 'data_' postFix ];
powData   = [ sessionDir filesep 'pow_' postFix ];
csdData   = [ sessionDir filesep 'csd_' postFix ];

% Get leadfield
% See getLeadfield.m for details
% What's requied for Fieldtrips beamformer is in "leadfield", a varaible
% called "ldf"
saveIn  = [ sessionDir filesep 'leadfield' ];
load( saveIn );  
clear trueLead;

% Get CSD files

% Pre condition
load( [ csdData 'bas' ], 'freq' );
cfg = [];
cfg.grid = source2sparse( ldf );
cfg.hdmfile = ldf.cfg.hdmfile;
cfg.frequency = foi;
cfg.method = 'power';
cfg.projectnoise = 'yes';
cfg.lambda = regularization;
cfg.feedback = 'gui';
sourcePre = sourceanalysis( cfg, freq );

% Post condition
load( [ csdData 'act' ], 'freq' ); 
cfg = [];
cfg.grid = source2sparse( ldf );
cfg.hdmfile = ldf.cfg.hdmfile;
cfg.frequency = foi;
cfg.method = 'power';
cfg.projectnoise = 'yes';
cfg.lambda = regularization;
cfg.feedback = 'gui';
sourcePost = sourceanalysis( cfg, freq );

clear freq;

% Plot on MRI using "sliceinterp"
%mriFile = 'E:\nienke\fieldtripCompare\leadfields\leadfields\mri\geijn_h_vd.mri';
cfg=[];
cfg.downsample = 2;
cfg.funparameter = 'pow';
cfg.feedback = 'gui';
%cfg.interactive = 'yes';
cfg.nslices = 28;   % default = 20, set in multiples of 4 
cfg.dim = 3;    % default = 3 = horizontal slice; 1 = saggital slice; 2 = coronal slice
%cfg.rotate = 3; % Set to 3 when cfg.dim = 1;
%cfg.alpha = 0.5;    % Controls opacity; choice = 0-1

sourceRelative = sourcePost; 
if 0
    pre = sourcePre.avg.pow;
    post = sourcePost.avg.pow;
else
    sourceRelative.avg.noise = zeros( size(sourceRelative.avg.noise) );
    sourceRelative.avg.csd = zeros( size(sourceRelative.avg.csd) );

    if 1
        load( ['lbex_leak']); 
        load( ['lbex_act']); post = trAvgPower; clear trAvgPower;
        load( ['lbex_bas']); pre = trAvgPower; clear trAvgPower;
    end
    
    if 0
        load( ['lbexFRatio_act'] );
    end
    
    if 0
        load( ['mnlsRes'] );
    end
    
    if 0
        load( ['mnls_act']); post = trAvgPower; clear trAvgPower;
        load( ['mnls_bas']); pre = trAvgPower; clear trAvgPower;    
    end
    
end
%ratio = 10*log10( trAvgFRatio );
%ratio = trAvgFRatio;
%ratio = leakage ;
%ratio = post;
%ratio = pre;
%ratio = post ./ pre;
%ratio = ratio/max(ratio);
%ratio = spread2;
ratio=repmat( realmin, size(post));
ratio(1) = 1.0;

clear sourcePre; clear sourcePost;
sourceRelative.avg.pow = ratio;

figure;
sourceRelF = source2full( sourceRelative );
sourceInterp = sourceinterpolate( cfg, sourceRelF, mriFile );
sliceinterp( cfg, sourceInterp );


% Output for ploting on MRI using "MRIcro"
if 0
    % write MRI anatomy in "analyze" format
    cfg=[];
    cfg.parameter = 'anatomy';
    cfg.fileprefix = strcat( [ csdData ], '_ana' );
    cfg.datatype = 'uint8';
    sourcewrite( cfg, sourceInterp );
    
    cfg=[];
    cfg.parameter = 'pow';
    cfg.fileprefix = strcat( [ csdData ], '_dat' );
    cfg.datatype = 'float';
    sourcewrite( cfg, sourceInterp );
    
end

disp( ['Frequency being displayed is: ' num2str( foi )] );
disp( ['Regularization Parameter: ' num2str( regularization )] );