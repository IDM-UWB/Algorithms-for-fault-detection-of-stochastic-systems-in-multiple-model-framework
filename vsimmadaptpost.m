function estimatevs = vsimmadaptpost(estimatevs,flagUnlikely,flagAdjacent,nModel,discardingStrategy,nMinModel)
%VSIMMADAPTPOST Perform set adaptation of the VS IMM - second part
%   estimatevs =
%   vsimmadaptpost(estimatevs,flagUnlikely,flagAdjacent,nModel,discardingStrategy,nMinModel)
%   finishes the adaptation of the active set model for the VS IMM by
%   discarding estimates in 'estimatevs' for models that are deemed
%   unnecessary based on the unlikely statuses given in 'flagUnlikely' and
%   'flagAdjacent'. The parameters that controls the discarding is the
%   total number of possible models 'nModels', the minimum number of models
%   in the active set ' nMinModel' and discarding stratyegy
%   'discardingStrategy' (1 - simple strategy that does not respect the
%   minimum number of models in the active set, 2 - discarding strategy
%   that keeps at minimum the required number)
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

% Flag discardable
flagDiscard = flagUnlikely & ~flagAdjacent;

% Double check that all models to be potentialy discarded has
% probabilities different from nans
if any(isnan(estimatevs.pmuf(flagDiscard)))
    error('Some models to be discared has nan as probabilities')
end

% Convert indices of next active to flag next active
flagActiveNext = false(nModel,1);
flagActiveNext(estimatevs.idxActiveNext) = true;



switch discardingStrategy
    case 1
        % Discard all models in discard set
        flagActiveNext = flagActiveNext & ~flagDiscard;
        
    case 2
        % Discard maximum number of models with lowest probabilities such
        % that resulting set has at least minimum number of models
        
        % Convert flag to idx for discarded models
        idxDiscard = find(flagDiscard);
        
        % Sort discard model from lowest to highest probability
        [~,idx] = sort(estimatevs.pmuf(flagDiscard));
        idxDiscardSorted = idxDiscard(idx);
        
        % Compute the number of models to be discarded
        nToDiscard = min(length(estimatevs.idxActiveNext)-nMinModel,length(idxDiscard));
        
        % Discard models if possible
        if nToDiscard>0
            flagDiscard = false(nModel,1);
            flagDiscard(idxDiscardSorted(1:nToDiscard)) = true;
            flagActiveNext = flagActiveNext & ~flagDiscard;
        end
end
estimatevs.idxActiveNext = find(flagActiveNext)';


end