function x = normrndm(mx,Px,m)
%NORMRNDM Random sample from normal distribution
%   X = NORMRNDM(MX,PX,M) generates random vectors from normal distribution
%   with mean value MX and positive definite covariance matrix PX.
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

if nargin < 3
    m = 1;
end

% Dimension of x
nx = length(mx);

% Compute Cholesky decomposition
sPx = chol(Px,'lower');

% Generate sample vectors
x = sPx*randn(nx,m) + mx(:,ones(1,m));