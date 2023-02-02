function [estimatevs,flagUnlikely,flagAdjacent] = vsimmadaptpre(estimatevs,estimatevsPrevious,model,thvs,uPrevious,yCurrent)
%VSIMMADAPTPRE Perform set adaptation of the VS IMM - first part
%   [estimatevs,flagUnlikely,flagAdjacent] =
%   vsimmadaptpre(estimatevs,estimatevsPrevious,model,thvs,uPrevious,yCurrent)
%   uses the estimate form the previsous time step 'estimatevsPrevious' and
%   current time step 'estimatevs' to perform adaptation of the set of
%   active models. It assumes that the model of the system is given in
%   'model', the thresholds 'thvs' are used for partitioning the set of
%   active models into three groups: unlikely (probability<thvs(1)),
%   principal  (probability>thvs(2)) and the rest. The input from the
%   previous time step should be provided in 'uPrevious'. The output from
%   the curren time step should be provided in 'yCurrent'. If new models
%   are added to the active set the filtering estimates are computed and
%   added to the 'estimatesvs'. The models that were flag unlikely are
%   returned in bolean vector 'flagUnlikely' and models flag as adjacent
%   are returned in bolean vector 'flagDajacent'
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

% Get the total number of models
nModel = length(model.M);

% Convert indices of currently active models into flag of currently active models
flagActiveCurrent = false(nModel,1);
flagActiveCurrent(estimatevs.idxActiveCurrent) = true;

% Flag principal models within currently active models
flagPrincipal = false(nModel,1);
flagPrincipal(estimatevs.idxActiveCurrent) = estimatevs.pmuf(estimatevs.idxActiveCurrent)>thvs(2);

% Flag unlikely models within currently active models
flagUnlikely = false(nModel,1);
flagUnlikely(estimatevs.idxActiveCurrent) = estimatevs.pmuf(estimatevs.idxActiveCurrent)<thvs(1);

% Flag models adjacent to principal models
flagAdjacent = any(model.P(:,flagPrincipal)>0,2); % 'any' takes care for empty matrix as well

% Flag new models (adjacent but not active)
flagNew = flagAdjacent;
flagNew(estimatevs.idxActiveCurrent) = false;

% If there are new models perform the whole prediction and
% filtering steps for new models and fuse estimates with estimates
% of original active models
if any(flagNew)
    % Perform prediction and filtering step for new models in the
    % set only
    estimatevsNewPrevious = estimatevsPrevious;
    estimatevsNewPrevious.idxActiveNext = find(flagNew');
    estimatevsNewCurrent = vsimmkfp(uPrevious,model,estimatevsNewPrevious);
    estimatevsNewCurrent = vsimmkff(yCurrent,model,estimatevsNewCurrent);
    
    % Add the filtering statistics (mean and covariance matrix) of
    % new models
    estimatevs.xf(:,estimatevsNewCurrent.idxActiveCurrent) =...
        estimatevsNewCurrent.xf(:,estimatevsNewCurrent.idxActiveCurrent);
    estimatevs.Pxxf(:,:,estimatevsNewCurrent.idxActiveCurrent) =...
        estimatevsNewCurrent.Pxxf(:,:,estimatevsNewCurrent.idxActiveCurrent);
    estimatevs.expyp(estimatevsNewCurrent.idxActiveCurrent) =...
        estimatevsNewCurrent.expyp(estimatevsNewCurrent.idxActiveCurrent);
    estimatevs.srdetPyp(estimatevsNewCurrent.idxActiveCurrent) =...
        estimatevsNewCurrent.srdetPyp(estimatevsNewCurrent.idxActiveCurrent);
    estimatevs.pmup(estimatevsNewCurrent.idxActiveCurrent)=...
        estimatevsNewCurrent.pmup(estimatevsNewCurrent.idxActiveCurrent);
    estimatevs.idxActiveCurrent = find(flagActiveCurrent | flagNew)';
    estimatevs.idxActiveNext = estimatevs.idxActiveCurrent;
    
    % Normalize probabilities
    estimatevs.pmuf(estimatevs.idxActiveCurrent) =...
        normalizeweightexp(estimatevs.expyp(estimatevs.idxActiveCurrent),...
        estimatevs.srdetPyp(estimatevs.idxActiveCurrent).*estimatevs.pmup(estimatevs.idxActiveCurrent)')';
end


end