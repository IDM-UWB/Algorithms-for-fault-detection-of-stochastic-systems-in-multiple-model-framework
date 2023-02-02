%   example_maneuveringtarget.m illustration of variable structure interacting
%   multiple model estimator for tracking of a manuevering target
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

close all
clear all
clc

% Vector of accelerations
% aLevel = {[-40 -20 0 20 40],[-40 -20 0 20 40]};
% naLevel = cellfun(@length,aLevel);
% n = prod(naLevel);
% deaggregationidx(aLevel,linidx2subidx(naLevel,1:n))

a = [0 20 0 -20 0 20 -20 -20 20 40 0 -40 0
    0 0 20 0 -20 20 20 -20 -20 0 40 0 -40];
nModel = size(a,2);

% figure
% plot(a(1,:),a(2,:))
% grid on
% hold on
% axis equal


% Dimension of state
nx = 6;

% Dimension of input
nu = 1;

% Dimension of output
ny = 2;

% Sampling period
Ts = 0.1;

for i = 1:nModel
    model.M(i).A = [1 0 Ts 0 0 0
        0 1 0 Ts 0 0
        0 0 1 0 0 0
        0 0 0 1 0 0
        0 0 0 0 1 0
        0 0 0 0 0 1];
    
    model.M(i).B = zeros(nx,1);
    
    Bq = [Ts^2/2 0
        0 Ts^2/2
        Ts 0
        0 Ts
        0 0
        0 0];
    model.M(i).q = Bq*a(:,i);
    
    model.M(i).G = [Ts^2/2 0
        0 Ts^2/2
        Ts 0
        0 Ts
        1 0
        0 1];
    model.M(i).Sigmaw = model.M(i).G*model.M(i).G';
    
    model.M(i).C = [1 0 0 0 0 0
        0 1 0 0 0 0];
    model.M(i).r = [0
        0];
    model.M(i).H = 0.01*eye(2);
    model.M(i).Sigmav = model.M(i).H*model.M(i).H';
end
model.P = [116/120 0.02 0.02 0.02 0.02 0 0 0 0 0 0 0 0
    1/120 0.95 0 0 0 1/30 0 0 1/30 0.1 0 0 0
    1/120 0 0.95 0 0 1/30 1/30 0 0 0 0.1 0 0
    1/120 0 0 0.95 0 0 1/30 1/30 0 0 0 0.1 0
    1/120 0 0 0 0.95 0 0 1/30 1/30 0 0 0 0.1
    0 0.01 0.01 0 0 28/30  0 0 0 0 0 0 0
    0 0 0.01 0.01 0 0 28/30 0 0 0 0 0 0
    0 0 0 0.01 0.01 0 0 28/30 0 0 0 0 0
    0 0.01 0 0 0.01 0 0 0 28/30 0 0 0 0
    0 0.01 0 0 0 0 0 0 0 0.9 0 0 0
    0 0 0.01 0.01 0 0 0 0 0 0 0.9 0 0
    0 0 0 0 0.01 0 0 0 0 0 0 0.9 0
    0 0 0 0 0 0 0 0 0 0 0 0 0.9];

model.xp0 = [100;100;0;0;0;0];
model.Pxxp0 = eye(6);
model.pmup0 = [1;zeros(nModel-1,1)];

% Define system to be the same as the model
sys.S = model.M;
sys.P = model.P;
sys.xp0 = model.xp0;
sys.Pxxp0 = model.Pxxp0;
sys.pmup0 = model.pmup0;

% Select discarding strategy of VS IMM
discardingStrategy = 2;

% Threshold for the variable structure set optimizations
thvs = [0.1 0.9];

% Minimum number of models of the variable structure IMM
nMinModel = 5;

% Number of times steps
F = 500;
Fp1 = F + 1;

% Preallocate arrays
x = nan(nx,Fp1);
mu = nan(1,Fp1);
u = nan(nu,F);
y = nan(ny,Fp1);
%estimatefs = nan(1,Fp1);
%estimatevs = nan(1,Fp1);
timefs = zeros(1,Fp1);
timevs = zeros(1,Fp1);

% Initial estimate of the FS IMM
estimatefs(1).idxActiveCurrent = [1 2 3 4 5 6 7 8 9 10 11 12 13];
estimatefs(1).idxActiveNext = [1 2 3 4 5 6 7 8 9 10 11 12 13];
estimatefs(1).xp = nan(nx,nModel);
estimatefs(1).Pxxp = nan(nx,nx,nModel);
estimatefs(1).pmup = nan(nModel,1);
estimatefs(1).xf = [];
estimatefs(1).Pxxf = [];
estimatefs(1).expyp = [];
estimatefs(1).srdetPyp = [];
estimatefs(1).pmuf = [];
for i = estimatefs(1).idxActiveCurrent
    estimatefs(1).xp(:,i) = model.xp0;
    estimatefs(1).Pxxp(:,:,i) = model.Pxxp0;
end
estimatefs(1).pmup(estimatefs(1).idxActiveCurrent) = model.pmup0(estimatefs(1).idxActiveCurrent);
% Check consistency
if abs(sum(estimatefs(1).pmup(estimatefs(1).idxActiveCurrent))-1)>1e-6
    error('active set is wong')
end

% Initialize estimate of the VS IMM
estimatevs(1).idxActiveCurrent = [1 2 3 4 5]; % This should be consistent with the minimum number of models of VSIMM
estimatevs(1).idxActiveNext = [1 2 3 4 5]; % This should be consistent with the minimum number of models of VSIMM
estimatevs(1).xp = nan(nx,nModel);
estimatevs(1).Pxxp = nan(nx,nx,nModel);
estimatevs(1).pmup = nan(nModel,1);
estimatevs(1).xf = [];
estimatevs(1).Pxxf = [];
estimatevs(1).expyp = [];
estimatevs(1).srdetPyp = [];
estimatevs(1).pmuf = [];
for i = estimatevs(1).idxActiveCurrent
    estimatevs(1).xp(:,i) = model.xp0;
    estimatevs(1).Pxxp(:,:,i) = model.Pxxp0;
end
estimatevs(1).pmup(estimatevs(1).idxActiveCurrent) = model.pmup0(estimatevs(1).idxActiveCurrent);
% Check consistency
if abs(sum(estimatevs(1).pmup(estimatevs(1).idxActiveCurrent))-1)>1e-6
    error('active set is wong')
end

% Generate initial condition and noises
x(:,1) = sys.xp0 + chol(sys.Pxxp0,'lower')*randn(nx,1);
mu(1) = gendrnd(sys.pmup0);
v = randn(2,Fp1);
w = randn(2,F);


for k = 1:Fp1
    y(:,k) = sys.S(mu(k)).C*x(:,k) + sys.S(mu(k)).r + sys.S(mu(k)).H*v(:,k);
    
    % Perform filtering step of the FS IMM estimator
    tic
    estimatefs(k) = vsimmkff(y(:,k),model,estimatefs(k));
    timefs(k) = timefs(k) + toc;
    
    % VS IMM
    
    % Perform filtering step of the VS IMM estimator
    tic
    estimatevs(k) = vsimmkff(y(:,k),model,estimatevs(k));
    
    
    if k>1
        % Perform model set adaptation of the VS IMM estimator (it cannot be performed at initial time)
        [estimatevs(k),flagUnlikely,flagAdjacent] = vsimmadaptpre(estimatevs(k),estimatevs(k-1),model,thvs,u(:,k-1),y(:,k));
    end
    
    % Compute global estimate of the VS IMM estimator
    flagActiveNonzero = false(1,nModel);
    flagActiveNonzero(estimatevs(k).idxActiveCurrent) =  estimatevs(k).pmuf(estimatevs(k).idxActiveCurrent)~=0;
    xfvsimm(:,k) = estimatevs(k).xf(:,flagActiveNonzero)*estimatevs(k).pmuf(flagActiveNonzero);
    pmufvsimm(:,k) = estimatevs(k).pmuf;
    [~,idxtmp] = max(estimatevs(k).pmuf(estimatevs(k).idxActiveCurrent));
    dvsimm(k) = estimatevs(k).idxActiveCurrent(idxtmp);
    nModelvs(k) = length(estimatevs(k).idxActiveCurrent);
    

    if k>1
        % Finish the model set addaptation of the VS IMM estimator (it cannot be performed at initial time)
        estimatevs(k) = vsimmadaptpost(estimatevs(k),flagUnlikely,flagAdjacent,nModel,discardingStrategy,nMinModel);
    end
    timevs(k) = timevs(k) + toc;
    
    % Compute global estimate of the FS IMM estimator
    tic
    flagActiveNonzero = false(1,nModel);
    flagActiveNonzero(estimatefs(k).idxActiveCurrent) =  estimatefs(k).pmuf(estimatefs(k).idxActiveCurrent)~=0;
    xffsimm(:,k) = estimatefs(k).xf(:,flagActiveNonzero)*estimatefs(k).pmuf(flagActiveNonzero);
    pmuffsimm(:,k) = estimatefs(k).pmuf;
    [~,dfsimm(k)] = max(estimatefs(k).pmuf(estimatefs(k).idxActiveCurrent));
    timefs(k) = timefs(k) + toc;
    
    
    %  Generate control using filtering estimates
    u(:,k) = 0;
    
    
    if k<Fp1
        % Simulate the system
        mu(k+1) = gendrnd(sys.P(:,mu(k)));
        x(:,k+1) = sys.S(mu(k+1)).A*x(:,k) + sys.S(mu(k+1)).B*u(:,k) + sys.S(mu(k+1)).q + sys.S(mu(k+1)).G*w(:,k);
        
        % Perform the prediction step of the FS IMM
        tic
        estimatefs(k+1) = vsimmkfp(u(:,k),model,estimatefs(k));
        timefs(k) = timefs(k) + toc;
        
        % Perform the prediction step of the VS IMM
        tic
        estimatevs(k+1) = vsimmkfp(u(:,k),model,estimatevs(k));
        timevs(k) = timevs(k) + toc;
    end
end

xffsimmerr = x-xffsimm;
xfvsimmerr = x-xfvsimm;

msefsimm = sum(sum(xffsimmerr.^2,1))/Fp1;
msevsimm = sum(sum(xfvsimmerr.^2,1))/Fp1;



tilefigure
plot(mu,'.')
hold on
grid on
plot(dfsimm,'s')
plot(dvsimm,'o')
legend({'true model','model estimate (FS)','model estimate (VS)'})
title('Mode')

tilefigure
plot(pmuffsimm')
hold on
grid on
plot(pmufvsimm',':')
title('Mode probabilities')

tilefigure
plot((pmuffsimm-pmufvsimm)')
grid on
title('Difference in model prob of FS IMM and VS IMM')


tilefigure
plot(x','.-')
grid on
hold on
plot(xffsimm','-s')
plot(xfvsimm','-o')
legend({'true state','state estimate (FS)','state estimate (VS)'})


tilefigure
plot(x(1,:),x(2,:),'.-')
grid on
hold on
plot(xffsimm(1,:),xffsimm(2,:),'-s')
plot(xfvsimm(1,:),xfvsimm(2,:),'-o')
legend({'true position','position estimate (FS)','position estimate (VS)'})

tilefigure
plot(xffsimmerr')
grid on
title(sprintf('FS IMM MSE = %f',msefsimm))

tilefigure
plot(xfvsimmerr')
grid on
title(sprintf('VS IMM MSE = %f',msevsimm))

tilefigure
plot(nModelvs)
grid on
title('Number of models used by VS IMM')

tilefigure
plot(timefs)
grid on
hold on
plot(timevs)
title('Time requirements of FS and VS IMM')
legend({'FS IMM','VS IMM'})