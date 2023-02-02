function estimateNew = vsimmkfp(u,model,estimate)
%VSIMMKFP Predictive step of variable structure interacting multiple model Kalman filter
%   estimateNew = vsimmkfp(u,model,estimate) performs the predictive step
%   of the variable structure interacting multiple model Kalman filter
%   using the input 'u', model 'model' and estimate 'estimate'
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
[nx,nModel] = size(estimate.xf);

% Initialize the estimate at next time step
estimateNew.xp = nan(nx,nModel);
estimateNew.Pxxp = nan(nx,nx,nModel);
estimateNew.pmup = nan(nModel,1);
estimateNew.Pxxf = [];
estimateNew.xf = [];
estimateNew.expyp = [];
estimateNew.srdetPyp = [];
estimateNew.pmuf = [];

% Compute joint probability of two step (current, future) model sequences
pmuj = model.P(estimate.idxActiveNext,estimate.idxActiveCurrent)'.*estimate.pmuf(estimate.idxActiveCurrent);

% Compute predictive probability of models at future time step
pmup = sum(pmuj,1);

% Compute smoothed probabilites
pmustmp  = pmuj./pmup;

% Set nans smoothed probabilities to zeros
flagNans = isnan(pmustmp);
pmustmp(flagNans)=0;

% Arrange mixing probabilities into a full matrix -- just due to
% implementation
pmus = nan(nModel,nModel);
pmus(estimate.idxActiveCurrent,estimate.idxActiveNext) = pmustmp;

% Compute mixing means
xs = nan(nx,nModel);
xs(:,estimate.idxActiveNext) = estimate.xf(:,estimate.idxActiveCurrent)*pmus(estimate.idxActiveCurrent,estimate.idxActiveNext);

% Compute mixing covariances, predictive means and predictive covariances
Pxxs = nan(nx,nx,nModel);
for i = estimate.idxActiveNext
    % Compute mixing covariance matrices
    Pxxs(:,:,i) = zeros(nx,nx);
    for j = estimate.idxActiveCurrent
        Pxxs(:,:,i) = Pxxs(:,:,i) + (estimate.Pxxf(:,:,j) + (xs(:,i)-estimate.xf(:,j))*(xs(:,i)-estimate.xf(:,j))')*pmus(j,i);
    end
    
    % Compute predictive means and covariance matrices
    [estimateNew.xp(:,i),estimateNew.Pxxp(:,:,i)] = kfp(u,xs(:,i),Pxxs(:,:,i),model.M(i).A,model.M(i).B,model.M(i).Sigmaw,model.M(i).q);
     % estimateNew.xp(:,i) = model.M(i).A*xs(:,i) + model.M(i).B*u + model.M(i).q;
     %estimateNew.Pxxp(:,:,i) = model.M(i).A*Pxxs(:,:,i)*model.M(i).A' + model.M(i).Sigmaw;
end

% Store predictive probabilites
estimateNew.pmup(estimate.idxActiveNext,1) = pmup';

% Set active sets for the next time step
estimateNew.idxActiveCurrent = estimate.idxActiveNext;
estimateNew.idxActiveNext = estimate.idxActiveNext; % The active set is changed outside this function, implicitly sctive set is not changing


end