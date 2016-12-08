% pop_bemobil_signal_decomposition() - Prepare data for decomposition (AMICA)
%
% Usage:
%   >>  [ EEG ] = pop_bemobil_signal_decomposition( EEG )
%
% Inputs:
%   EEG                   - eeglab EEG struct%    
%
% Outputs:
%   EEG     - EEGLAB EEG structure with ICA weight and sphere matrix
%
% See also: 
%   RUNAMICA15, POP_RUNAMICA EEGLAB

function [ EEG ] = pop_bemobil_signal_decomposition( EEG )

if nargin < 1
	help pop_bemobil_signal_decomposition;
	return;
end;

% add preprocess GUI
uigeom_cols = 2;
uigeom_rows = 4;

uigom = {
    {uigeom_cols uigeom_rows [0 0] [1 1]} 
    {uigeom_cols uigeom_rows [1 0] [1 1]}

    {uigeom_cols uigeom_rows [0 1] [1 1]} 
    {uigeom_cols uigeom_rows [1 1] [1 1]}

    {uigeom_cols uigeom_rows [0 2] [1 1]} 
    {uigeom_cols uigeom_rows [1 2] [1 1]}
    
    {uigeom_cols uigeom_rows [0 3] [1 1]} 
    {uigeom_cols uigeom_rows [1 3] [1 1]}
    
    %{uigeom_cols uigeom_rows [0 4] [1 1]} 
    %{uigeom_cols uigeom_rows [1 4] [1 1]}
    };

uilist = {
    { 'style', 'text', 'string', 'Iteration'}
    { 'style', 'edit', 'string', '1', 'tag', 'iteration'}

    { 'style', 'text', 'string', 'AMICA yes/no'}
    { 'style', 'checkbox', 'value', 1, 'tag', 'amica'}
    
    { 'style', 'text', 'string', 'AMICA: number of models'}
    { 'style', 'edit', 'string', '1', 'tag', 'num_models'}
    
    { 'style', 'text', 'string', 'AMICA: max_threads'}
    { 'style', 'edit', 'string', '4', 'tag', 'max_threads'}

    %{ 'style', 'text', 'string', 'will add others here later'}
    %{ 'style', 'edit', 'string', 'other decomp algorithm'}
    };

[ tmp1 tmp2 strhalt structout ] = inputgui('geom', uigom, 'uilist', uilist, ...
    'title', 'Signal decomposition', 'helpcom', 'pophelp(''pop_bemobil_signal_decomposition'');');
    
if isempty(strhalt)
    return; 
end;

% close menu figure after OK
close(gcbf);

% prepare values from input gui
%structout.iteration = str2double(structout.iteration);
structout.num_models = str2double(structout.num_models);
structout.max_threads = str2double(structout.max_threads);

% run processing function with values from gui input
EEG = bemobil_signal_decomposition(EEG, structout.iteration, structout.amica,...
    structout.num_models, structout.max_threads); %structout.other);
end

