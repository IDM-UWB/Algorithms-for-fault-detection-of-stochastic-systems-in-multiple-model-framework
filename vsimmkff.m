function estimate = vsimmkff(y,model,estimate)
%VSIMMKFF Filtering step of variable structure interating multiple model Kalman filter
%   estimate = vsimmkff(y,model,estimate) performs the filtering step
%   of the variable structure interacting multiple model Kalman filter
%   using the output 'y', model 'model' and estimate 'estimate'
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

% Get dimension of state and number of models
[nx,nModel] = size(estimate.xp);

% Initialize filtering statistics to nans
estimate.xf = nan(nx,nModel);
estimate.Pxxf = nan(nx,nx,nModel);
estimate.expyp = nan(1,nModel);
estimate.srdetPyp = nan(1,nModel);
estimate.pmuf = nan(nModel,1);

% Perform the filtering step for curently active models only
for i = estimate.idxActiveCurrent
    
    % Compute filtering estimates using Kalman filter
    [estimate.xf(:,i),estimate.Pxxf(:,:,i),~,Pyyp,~,e] = ...
        kff(y,estimate.xp(:,i),estimate.Pxxp(:,:,i),model.M(i).C,model.M(i).Sigmav,model.M(i).r);
    estimate.expyp(i) = -0.5*e'*(Pyyp\e);
    estimate.srdetPyp(i) = 1/sqrt(det(Pyyp)); % Note that constant 2*pi is not included as it gets cancelled during weight computation
end

% Normalize filtering probabilities
estimate.pmuf(estimate.idxActiveCurrent) = normalizeweightexp(estimate.expyp(estimate.idxActiveCurrent),estimate.srdetPyp(estimate.idxActiveCurrent).*estimate.pmup(estimate.idxActiveCurrent)')';


end