% This routine loads crx data. Snippets are treated as trials. No meta
% data has been returned. Data is of raw channels.
function varargout = getCrxData(varargin)  
% get file name and load
[crxFileName] = getFileToBeLoaded(varargin{:});
if isempty(crxFileName)
    disp('try again...');
    return;
end
try
    load(crxFileName, '-mat', 'Info');
    numSavedSnippets = Info.Count;
catch
    disp('Unable to read file.');
    return
end
oData = cell(numSavedSnippets,1);
% loop over num of snippets
for nSnippetCount = 1 : numSavedSnippets
    strSnippetName = sprintf('Snippet_%d', nSnippetCount);
    SnippetData = load(crxFileName, '-mat', strSnippetName);
    %packaging for output
    oData{nSnippetCount} = SnippetData.(strSnippetName).Data.modData{1};
end
varargout{1} = oData;

%% get file name to be loaded if not provided as command line argument.
function [crxFileName] = getFileToBeLoaded(varargin)
try
    crxFileName = [];
    if ~nargin
        [filename pathname] = uigetfile('*.crx','Select .crx file');
        if ~isequal(filename,0)
            crxFileName = fullfile(pathname, filename);
        end
    else
        crxFileName = varargin{1};
    end
catch
    disp('Unable to get the file.');
end