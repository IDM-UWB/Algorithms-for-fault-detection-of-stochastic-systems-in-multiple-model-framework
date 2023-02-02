function [pmuseqnew,idxMaxProbability] = ps2pfs(pmuseq,nSequencesNew)
%PS2PFS Computes probabilities of terminal subsequences
%   [pmuseqnew,idxMaxProbability] = ps2pfs(pmuseq,nSequencesNew) computes probabilities of
%   terminal
%   subsequences 'pmuseqterminal' given the probabilities of individual
%   sequences 'pmuseq' and the number of required terminal sequences 'n'.
%   Note that the algorithm do not check if nSequencesNew is an integer
%   power of the number of models.
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

% Number of model sequences
nSequences = length(pmuseq);

% Check that the number of required sequences is the same or lower than
% current number of sequnces
if nSequencesNew < nSequences
    % Reshape the column vector of probabilities to get matrix P, where
    % each column of P contains probabilities of model sequences with the
    % same final subsequence
    pmuseqReshaped = reshape(pmuseq,nSequencesNew,[]);
    
    % The probabilities of model sequences that have the same final sequence
    pmuseqnew = sum(pmuseqReshaped,2);
    
    % Index of model sequence that has maximum probability within group of sequences with common final model sequence
    [~,idxMaxProbability] = max(pmuseqReshaped,[],2);
    idxMaxProbability = (idxMaxProbability-1)*nSequencesNew + (1:nSequencesNew)';
else
    pmuseqnew = pmuseq;
    [~,idxMaxProbability] = max(pmuseq,[],2);
end


end