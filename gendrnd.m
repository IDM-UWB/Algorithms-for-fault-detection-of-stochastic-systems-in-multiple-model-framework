function x = gendrnd(p,varargin)
% GENDRND  Generate random array using user defined discrete distribution.
%
% X = GENDRND(p,m,n,values,cdfTol) returns an M-by-N matrix containing
% pseudorandom values drawn from the set VALUES (default set {1, 2, 3, ...,
% K-1, K}) with user defined probability ditribution P over the set {1, 2,
% 3, ..., K-1, K}.
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

sizp = size(p);
if numel(sizp)>2 || sizp(2)>1 || ~isreal(p)
    exception = MException('MyToolbox:gendrnd:wrongReqInArgs',...
        'Argument P has to be a column vector of real numbers.');
    throw(exception);
end

% Check the maximum number of optional input arguments
maxnumvarargs = 4;
numvarargs = length(varargin);
if numvarargs > maxnumvarargs
    exception = MException('MyToolbox:gendrnd:tooManyOptInArgs',...
        'The function accepts at most %d optional input arguments.',maxnumvarargs);
    throw(exception);
end

% Skip any new inputs if they are empty
newVals = cellfun(@(x) isnotempty(x), varargin);

% Set defaults for optional input arguments
optargs = {1,1,(1:length(p))',sqrt(eps)};

% Now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
optargs(newVals) = varargin(newVals);

% Place optional input arguments in memorable variable names
[m,n,values,cdfTol] = optargs{:};

% Compute cumulative distribution
cdfp = [0;cumsum(p)];

% Check that the conditions for p to be a probability vector are satisfied
if any(p < 0) || abs(cdfp(end)-1) > cdfTol
    exception = MException('MyToolbox:gendrnd:notValidDistribution',...
        'Argument P has to contain nonnegative real numbers that sum up to 1 to be a valid pmf.');
    throw(exception);
end

% Generates random numbers from the open interval (0,1)
randomNumbers = rand(m*n,1);

% Find the indices
%[~,~,idx] = histcounts(randomNumbers,cdfp);
idx = discretize(randomNumbers,cdfp); % Discretize seems to be slightly faster than histcounts

% Index the values and reshape them
x = values(idx);
x = reshape(x,m,n);


end