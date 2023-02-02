function [pmuseqmerged,xseqmerged,Pxxseqmerged] = mmmerge(pmuseq,xseq,Pxxseq,nSequencesMerged)
%MMMERGE Computes statistics of terminal subsequences
%   [pmuseqmerged,xseqmerged,Pxxseqmerged] =
%   mmmerge(pmuseq,xseq,Pxxseq,nSequencesMerged) computes the mean values,
%   covariance matrices and probabilities by merging terminal subsequnces
%   of multiple hypothesis tree
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

% The dimension of the state  and the number of model sequences
[nx,nSequences] = size(xseq);

if nSequencesMerged < nSequences
    % Compute probabilities of model sequences after merging
    pmuseqmerged = ps2pfs(pmuseq,nSequencesMerged);
    
    % Compute the mean values of the merged sequences
    xseqmerged = reshape(sum(reshape(bsxfun(@times,xseq,pmuseq'),nSequencesMerged*nx,[]),2),nx,nSequencesMerged);
    xseqmerged = bsxfun(@rdivide,xseqmerged,pmuseqmerged');
    
    % Replace nan that results from zero division by zeros
    xseqmerged(isnan(xseqmerged)) = 0;
    
    
    % Compute the covariance matrices of the merged sequences
    % Preallocate arays
    Pw = zeros(nx,nx,nSequences);
    Pxxseqmerged = zeros(nx,nx,nSequencesMerged);
    
    dxseq = xseq - repmat(xseqmerged,1,nSequences/nSequencesMerged);
    for i = 1:nSequences
        Pw(:,:,i) = (Pxxseq(:,:,i) + dxseq(:,i)*dxseq(:,i)')*pmuseq(i);
    end
    
    for i = 1:nSequencesMerged
        if pmuseqmerged(i) > 0
            Pxxseqmerged(:,:,i) = sum(Pw(:,:,i:nSequencesMerged:end),3)/pmuseqmerged(i);
        else
            % Prevent getting nan as a result of division by zero, set up
            % to the identity matrix. The matrix can be chosen arbitrary as
            % the probability of the sequqnce is zero anyway
            Pxxseqmerged(:,:,i) = eye(nx,nx);
        end
    end
else
    pmuseqmerged = pmuseq;
    xseqmerged = xseq;
    Pxxseqmerged = Pxxseq;
end


end