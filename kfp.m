function [xp,Pxxp] = kfp(u,xf,Pxxf,A,B,Sigmaw,meanw)
%KFP Predictive step of the Kalman filter
%   [xp,Pxxp] = kfp(u,xf,Pxxf,A,B,Sigmaw,meanw) computes predictive mean
%   value 'xp' and predictive covariance matrix 'Pxxp', based on the input
%   'u', filtering mean value 'xf', filtering covariance matrix 'Pxxf',
%   state matrix 'A', input matrix 'B', noise covariance matrix 'Sigmaw',
%   and noise mean value 'meanw'
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

% If mean value of state noise is not given, zero is assumed
if nargin<7
    nx = size(xf,1);
    meanw = zeros(nx,1);
end

isemptyU = isempty(u);
isemptyB = isempty(B);
if (isemptyB && ~isemptyU) || (~isemptyB && isemptyU)
    error('Input matrix B and input u have to be both either empty or nonempty')
end

if isemptyU
    xp = A*xf + meanw;
else
    xp = A*xf + B*u + meanw;
end
Pxxp = A*Pxxf*A' + Sigmaw;


end