function fout = tilefigure(varargin)
%TILEFIGURE Creates tiled figures
%   fh = tilefigure(tileSize,key,value,...) starts creating tiled figures
%   from top-left corner. The size of figure is determined by tileSize =
%   [number of columns, number of rows] (default value is tileSize = [4,
%   3]). The tile size can be followed by the key-value pairs supported by
%   figure command. The tiling is restarted when different tileSize is
%   given
%
%   Copyright (C) 2022  Ivo Puncochar
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

persistent tileSize numberOfDisplayedFigures

% Get number of optional arguments
numberOfVarArgs = length(varargin);

if numberOfVarArgs>0
    % If some input arguments are given thre are two posibilities 1) first
    % argument is vector of tile size possibly folowed by key-value pairs
    % 2) only key-value pairs are given
    if isnumeric(varargin{1}) && length(varargin{1})==2 && mod(numberOfVarArgs-1,2)==0 && iscellstr(varargin(2:2:end))
        
        % Reset tiles if different tile size is given
        if ~isequal(tileSize(:),varargin{1}(:))
            tileSize = varargin{1};
            numberOfDisplayedFigures = 0;
        end
        
        % Remove the first argument from varargin
        varargin(1) = [];
    elseif mod(numberOfVarArgs,2)~=0 || ~iscellstr(varargin(1:2:end))
        error('Incorrect input arguments')
    end
        
end

% If the function is called for the first time and tile size input argument
% is not given,set default tile size and set the number of displayed
% figures to zero
if isempty(tileSize)
    tileSize = [4 3];
    numberOfDisplayedFigures = 0;
end


% Compute maximum number of figures
maxFigures = tileSize(1)*tileSize(2);


if numberOfDisplayedFigures == maxFigures
    error('Cannot display more than %i figures.',maxFigures)
end

% Compute normalized width and height of figure
normWH = 1./tileSize;

% Increase the number of diplayed figures
numberOfDisplayedFigures = numberOfDisplayedFigures + 1;

% Compute indices of tile for current figure
[col,row] = ind2sub(tileSize,numberOfDisplayedFigures);

% Note the screen is indexed bottom-up, therefore tileSize(2)-row is used
% to get top-down ordering

% Display figure
f = figure('Units','normalized',...
    'OuterPosition',[normWH(1)*(col-1) normWH(2)*(tileSize(2)-row) normWH],...
    varargin{:});

% Set Units of figure to default 'pixels'
set(f,'Units','pixels');


% Assign output argument if required
if nargout==1
    fout = f;
end


end