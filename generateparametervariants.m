function [indices,combinationsTransProbabilities] = generateparametervariants(paramTransitionProbs)
%GENERATEPARAMETERVARIANTS Generates all possible combibation of multiple
%parameter variants
%   [indices,combinationsTransProbabilities] = generateParameterVariants(paramTransitionProbs)
%   Input:  structure with n atributes containing matices of transitional 
%           probabilities of individual parameter variants
%   Output: matrix containing all indces of all variants combination 
%           matrix of combined transitional probabilitis of all parameters
%           variants (considering than than only one index can change by
%           the value of 1)
%
%   Copyright (C) 2022  Miroslav Flidr, Ivo Puncochar
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

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

% creation of an index matrix for combinations indexing
paramTransitionProbsMatrices = struct2cell(paramTransitionProbs);                       % "extraction" of transition probabilities matrices into one cell array
dims = cellfun(@length,paramTransitionProbsMatrices)';                                  % determine dimension of transition probabilities matrices 
numberOfParameters = length(dims);                                                      % determine the number of considered parameters
indexVectors = arrayfun(@linspace,ones(1,numberOfParameters),dims,dims,'UniformOutput',false);  % generation of all posible index combinations
[indexVectors{:}]=ndgrid(indexVectors{:});                                              %
indices = reshape(cat(numberOfParameters+1,indexVectors{:}),[],numberOfParameters);     % generation of all unique index combinations

% determine transitional probabilities matrix of model combination
combinationsTransProbabilities = paramTransitionProbsMatrices{1};
for parameter = 2:numberOfParameters
    combinationsTransProbabilities = kron(combinationsTransProbabilities,paramTransitionProbsMatrices{parameter});
end

% reduction of possible model transitions
[i,j,val] = find(combinationsTransProbabilities);
mask = sum(abs(indices(i,:)-indices(j,:)),2)>1;
i(mask) = [];
j(mask) = [];
val(mask) = [];
combinationsTransProbabilities = sparse(i,j,val);


columnSum = sum(combinationsTransProbabilities, 1);                    % evaluace sum of all columns
combinationsTransProbabilities = bsxfun(@rdivide, combinationsTransProbabilities, columnSum); % sum of transition probabilities in columns has to be one
