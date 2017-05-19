function [ usrlabels ] = lmx( config )
% This routine reads from a usr provided xml file providing which channel 
% numbers are to be relabelled, with the original label and the new label 
% The xml schema is as used by Niko's group, the original channels are 
% very likely Xltek specific. 
%
% This routine uses the open source xml package called "xmltree"

rname = [mfilename ':: ']; 
fprintf( 1, [rname '\n'] ); 
if ~isfield( config, 'filename' ), error( [rname 'Input config.filename not specified...' ] ); end
if ~isfield( config, 'dir' ), config.dir = pwd; fprintf( 1, '   Directory not specified... set as: %s\n', config.dir ); end

tree = xmltree( [ config.dir filesep config.filename ] );

% Below, the paths are xml schema specific and will need to be changed if
% the xml schema is changed. 

% Note: 
% usrlabels.original(1) returns a cell type
% usrlabels.original{1} returns a char type
% Holds true for *.channelNumbers and *.new too.

uid = find(tree,'/ReLabelDetails/chngLbl/ndx');
usrlabels.channelNumbers = get(tree,children(tree,uid),'value') ;

uid = find(tree,'/ReLabelDetails/chngLbl/orig');
usrlabels.original = upper( get(tree,children(tree,uid),'value') );

uid = find(tree,'/ReLabelDetails/chngLbl/chng');
usrlabels.new = upper( get(tree,children(tree,uid),'value') );
