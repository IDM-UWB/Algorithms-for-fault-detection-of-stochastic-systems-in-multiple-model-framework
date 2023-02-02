%   example_faultdetection.m illustration of variable structure interacting
%   multiple model estimator for multiple fault detection
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

% Model
% Sampling period
Ts = 1;

% Dynamic matrix
A = [1 Ts
    0 1];

% Measurement matrix
C = [1 0
    2 0
    3 0
    4 0];

% Intensity of the continuous-time acceleration process noise
q = 0.1;

% Discrete-time noise matrices
G = sqrt(q)*[sqrt(Ts^3/3) 0
    sqrt(3*Ts)/2 sqrt(Ts)/2];
H = sqrt(0.01)*eye(4);

% Covariance matrices of state and measurement noise
SigmaGw = G*G';
SigmaHv = H*H';


% Get dimensions
[ny,nx] = size(C);
nw = size(G,2);
nv = size(H,2);

% Initial condition
xp0 = [5
    0.02];
Pxxp0 = diag([10,0.01]);

% Desired probability of false alarm
Pfa = 1e-3;

% The simulation horizon
N = 500;
Np1 = N + 1;


% Simulate the measurement fault
idxC = 1;
f = zeros(ny,Np1);
f(1,50:149) = 3;
f(2,100:199) = (0:99)*0.04;
f(3,120:219) = 2;
%f(3,170:269) = 4*sin(0.1*(0:99));
%f(:,300:399) = C(4:-1:1,idxC*ones(1,100));
f(:,300:309) = [1;0;0;0]*ones(1,10);
f(:,310:319) = [1;1;0;0]*ones(1,10);
f(:,320:329) = [1;2;0;0]*ones(1,10);
f(:,330:339) = [1;2;1;0]*ones(1,10);
f(:,340:349) = [1;2;2;0]*ones(1,10);
f(:,350:359) = [1;2;3;0]*ones(1,10);
f(:,360:369) = [1;2;3;1]*ones(1,10);
f(:,370:379) = [1;2;3;2]*ones(1,10);
f(:,380:389) = [1;2;3;3]*ones(1,10);
f(:,390:419) = [1;2;3;4]*ones(1,30);


% Residual generator design using the measurement equation only
WstaticAux = null(C')';
tmp = chol(WstaticAux*SigmaHv*WstaticAux','lower');
Wstatic = inv(tmp)*WstaticAux;

% Compute threshold for statistic
nrstatic = size(Wstatic,1);
SstaticTh = chi2inv(1-Pfa,nrstatic);

% Residual genearator design using the whole model
Nwindow = 1;
O = [];
T = zeros((Nwindow+1)*nx,Nwindow*nw);
for i = 0:Nwindow
    O = [O;C*A^i];
    if i>0
        T(i*nx+1:(i+1)*nx,:) = A*T((i-1)*nx+1:i*nx,:);
        T(i*nx+1:(i+1)*nx,(i-1)*nw+1:i*nw) = G;
    end
end
T = kron(eye(Nwindow+1),C)*T;
WdynamicAux = null(O')';
SigmaY = T*kron(eye(Nwindow),SigmaGw)*T' + kron(eye(Nwindow+1),SigmaHv);
tmp = chol(WdynamicAux*SigmaY*WdynamicAux','lower');
Wdynamic = inv(tmp)*WdynamicAux;

% Compute threshold for statistic
nrdynamic = size(Wdynamic,1);
SdynamicTh = chi2inv(1-Pfa,nrdynamic);

% Direct fault estimation inspired by Zhong
WstaticAux = null(C')';
EstaticUnc = WstaticAux'*inv(WstaticAux*(SigmaHv+eye(ny))*WstaticAux')*WstaticAux;


% Define augmented model for Kalman filter
Aaug = [A zeros(nx,ny)
    zeros(ny,nx) eye(ny)];
Caug = [C eye(ny)];
Gaug = [G zeros(nx,ny)
    zeros(ny,nx) 1*eye(ny)];
Haug = H;
SigmaHvaug = Haug*Haug';
SigmaGwaug = Gaug*Gaug';

% Compute theoretical detectability
detectabilityIndex = rank(obsv(Aaug,Caug)); % in terms of augmented matrices
if detectabilityIndex<(nx+ny)
    fprintf('Not detectable according to Caglayan criterion\n')
end
%Ktmp = -place(A',C',[0.2 0.3])'; % in terms of original matrices (should be the same)
%eig(A+Ktmp*C)
%rank(C*inv(eye(nx)-(A+Ktmp*C))*Ktmp*eye(ny)+eye(ny))


% Define models for the multiple model estimator
scenario = 2;
switch scenario
    case 1
        model.M(1).A = A;
        model.M(1).B = [];
        model.M(1).G = G;
        model.M(1).q = zeros(nx,1);
        model.M(1).C = C;
        model.M(1).H = H;
        model.M(1).r = zeros(ny,1);
        
        model.M(2).A = A;
        model.M(2).B = [];
        model.M(2).G = G;
        model.M(2).q = zeros(nx,1);
        model.M(2).C = C;
        model.M(2).H = H;
        model.M(2).r = C(:,idxC);
        
        
        model.P = [0.5 0.5
            0.5 0.5];
        model.pmup0 = [0.5;0.5];
        model.xp0 = xp0;
        model.Pxxp0 = Pxxp0;
        
        estimate.xpseq = model.xp0;
        estimate.Pxxpseq = model.Pxxp0;
        estimate.pmupseq = model.pmup0;
        estimate.ypseq = [];
        estimate.Pyypseq = [];
        estimate.xfseq = [];
        estimate.Pxxfseq = [];
        estimate.pmufseq = [];
    case 2
        model = createmultiplemodels(A,[],G,zeros(nx,1),C,H,[0:4 -4:-1],0.2);
        model.pmup0 = zeros(size(model.P,1),1);
        model.pmup0(1) = 1;
        model.xp0 = xp0;
        model.Pxxp0 = Pxxp0;
        
        
        % Select discarding strategy of VS IMM
        discardingStrategy = 2;
        
        % Threshold for the variable structure set optimizations
        thvs = [0.1 0.9];
        
        % Minimum number of models of the variable structure IMM
        nMinModel = 9;
        
        % Initialize estimate of the VS IMM
        nModel = length(model.M);
        estimatevs(1).idxActiveCurrent = find(model.P(:,1))'; % This should be consistent with the minimum number of models of VSIMM
        estimatevs(1).idxActiveNext = find(model.P(:,1))'; % This should be consistent with the minimum number of models of VSIMM
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
        
        
    otherwise
        error('unknown scenario')
end


% Preallocate arrays
x = nan(nx,Np1);
y = nan(ny,Np1);

rstatic = nan(nrstatic,Np1);
Sstatic = nan(1,Np1);
dstatic = nan(1,Np1);
festStaticCon = zeros(ny,Np1); % fault estimate with assumption maximum two faults
festStaticUnc = zeros(ny,Np1); % fault estimate without assumption of maximum two faults

rdynamic = nan(nrdynamic,Np1);
Sdynamic = nan(1,Np1);
ddynamic = nan(1,Np1);

xaugp = nan(nx+ny,Np1);
xaugf = nan(nx+ny,Np1);
Pxxaugp = nan(nx+ny,nx+ny,Np1);
Pxxaugf = nan(nx+ny,nx+ny,Np1);


% Initialze system
x(:,1) = normrndm(xp0,Pxxp0);
w = normrndm(zeros(nx,1),SigmaGw,N);
v = normrndm(zeros(ny,1),SigmaHv,Np1);


% Initialize Kalman filter for augmented model
xaugp(:,1) = [xp0
    zeros(ny,1)];
Pxxaugp(:,:,1) = [Pxxp0 zeros(nx,ny)
    zeros(ny,nx) 100*eye(ny)];

for k = 1:Np1
    fprintf('Time instant %i/%i\n',k,Np1)
    
    % Measurement simulation
    y(:,k) = C*x(:,k) + f(:,k) + v(:,k);
    
    % Compute direct estimate of the fault
    festStaticUnc(:,k) = EstaticUnc*y(:,k);
    
    % Compute static residual
    rstatic(:,k) = Wstatic*y(:,k);
    
    % Compute static statistics
    Sstatic(k) = rstatic(:,k)'*rstatic(:,k);
    
    if Sstatic(k)<=SstaticTh
        dstatic(k) = 0;
        
    else
        dstatic(k) = 1;
                % Try to identify the measurements with faults assume at most two
                % faults
                % Compute distances to all planes
                indices = [1 2
                    1 3
                    1 4
                    2 3
                    2 4
                    3 4];
        
                for i = 1:size(indices,1)
                    Wtmp = Wstatic(:,indices(i,:));
        
                    pd(i) = rstatic(:,k)'*(eye(nrstatic)-Wtmp*inv(Wtmp'*Wtmp)*Wtmp')*rstatic(:,k);
                end
                [~,imin] = min(pd);
                fidx(:,k) = indices(imin,:)';
                Wtmp = Wstatic(:,indices(imin,:));
                festStaticCon(indices(imin,:),k) = inv(Wtmp'*Wtmp)*Wtmp'*rstatic(:,k);
    end
    
    % Dynamic fault detector
    if k>Nwindow
        tmp = y(:,k-Nwindow:k);
        rdynamic(:,k) = Wdynamic*tmp(:);
        Sdynamic(k) = rdynamic(:,k)'*rdynamic(:,k);
        if Sdynamic(k)<=SdynamicTh
            ddynamic(k) = 0;
        else
            ddynamic(k) = 1;
        end
        
    end
    
    % Filtering step of Kalman filter for augmented model
    [xaugf(:,k),Pxxaugf(:,:,k)] = kff(y(:,k),xaugp(:,k),Pxxaugp(:,:,k),Caug,SigmaHvaug,zeros(ny,1));
    
    
    % Filtering step of multiple model state estimator
    switch scenario
        case 1
            
            estimate(k) = smmkff(y(:,k),estimate(k),model,false);
            [estimate(k).pmufseq,estimate(k).xfseq,estimate(k).Pxxfseq] =...
                mmmerge(estimate(k).pmufseq,estimate(k).xfseq,estimate(k).Pxxfseq,length(model.M));
            prst1(k) = estimate(k).pmufseq(1);
            
        case 2
            estimatevs(k) = vsimmkff(y(:,k),model,estimatevs(k));
            
            
            if k>1
                % Perform model set adaptation of the VS IMM estimator (it cannot be performed at initial time)
                [estimatevs(k),flagUnlikely,flagAdjacent] = vsimmadaptpre(estimatevs(k),estimatevs(k-1),model,thvs,[],y(:,k));
            end
            
            % Compute global estimate of the VS IMM estmator
            flagActiveNonzero = false(1,nModel);
            flagActiveNonzero(estimatevs(k).idxActiveCurrent) =  estimatevs(k).pmuf(estimatevs(k).idxActiveCurrent)~=0;
            xfvsimm(:,k) = estimatevs(k).xf(:,flagActiveNonzero)*estimatevs(k).pmuf(flagActiveNonzero);
            pmufvsimm(:,k) = estimatevs(k).pmuf;
            [~,idxtmp] = max(estimatevs(k).pmuf(estimatevs(k).idxActiveCurrent));
            dvsimm(k) = estimatevs(k).idxActiveCurrent(idxtmp);
            nModelvs(k) = length(estimatevs(k).idxActiveCurrent);
            
            % Compute the estimate of measuremetn fault
            festvsimm(:,k) = model.M(dvsimm(k)).r;
            
            if k>1
                % Continue with the VS IMM estimator
                
                estimatevs(k) = vsimmadaptpost(estimatevs(k),flagUnlikely,flagAdjacent,nModel,discardingStrategy,nMinModel);
            end
            
    end
    
    
    
    
    if k<Np1
        % Simulate system state
        x(:,k+1) = A*x(:,k) + w(:,k);
        
        % Predictive step of the Kalman filter for the augmented model
        [xaugp(:,k+1),Pxxaugp(:,:,k+1)] = kfp([],xaugf(:,k),Pxxaugf(:,:,k),Aaug,[],SigmaGwaug,zeros(nx+ny,1));
        
        switch scenario
            case 1
                % Predictive step of the multiple model state estimator
                %estimate(k+1) = smmkfp([],estimate(k),model);
            case 2
                estimatevs(k+1) = vsimmkfp([],model,estimatevs(k));
            otherwise
                error('unknown scenario')
        end
    end
end


tk = 0:N;

figure
plot(tk,x)
grid minor
xlabel('Time instant $k$ [-]','Interpreter','latex')
ylabel('State $\mathbf{x}_{k}$','Interpreter','latex')


figure
plot(tk,x)
hold on
set(gca,'ColorOrderIndex',1)
plot(tk,xaugf(1:2,:),'.-')
grid minor
xlabel('Time instant $k$ [-]','Interpreter','latex')
ylabel('State $\mathbf{x}_{k}$','Interpreter','latex')



figure
plot(tk,y)
grid minor
xlabel('Time instant $k$ [-]','Interpreter','latex')
ylabel('Measurement $\mathbf{y}_{k}$','Interpreter','latex')

figure('Name','fault-modelnumber')
subplot(1,2,1)
plot(tk,f,'Linewidth',1.5)
xlim([0,N])
grid minor
xlabel('Time instant $k$','Interpreter','latex')
ylabel('Measurement fault $\mathbf{f}_{k}$','Interpreter','latex')
legend({'fault $\mathbf{f}_{k,1}$','fault $\mathbf{f}_{k,2}$','fault $\mathbf{f}_{k,3}$','fault $\mathbf{f}_{k,4}$'},'Interpreter','latex')
subplot(1,2,2)
plot(tk,nModelvs,'LineWidth',1.5,'Marker','x','LineStyle','none')
ylim([8 17])
%set(gca,'YTick', 8:17)
grid on
xlabel('Time instant $k$','Interpreter','latex')
ylabel('Number of active models of VSIMM','Interpreter','latex')
%myprint

figure('Name','faultestimation')
subplot(1,2,1)
plot(tk,f,'LineWidth',1.5,'LineStyle',':')
grid minor
hold on
set(gca,'ColorOrderIndex',1)
plot(tk,festStaticCon,'LineWidth',1.5)
xlim([0,500])
xlabel('Time instant $k$','Interpreter','latex')
ylabel('Measurement fault','Interpreter','latex')
% subplot(2,2,3)
% plot(tk,f)
% grid minor
% hold on
% set(gca,'ColorOrderIndex',1)
% plot(tk,festStaticUnc,'.:')
% xlabel('Time instant $k$ [-]','Interpreter','latex')
% ylabel('Measurement fault','Interpreter','latex')
subplot(1,2,2)
plot(tk,f,'LineWidth',1.5,'LineStyle',':')
grid minor
hold on
set(gca,'ColorOrderIndex',1)
plot(tk,festvsimm,'LineWidth',1.5)
xlim([0,500])
xlabel('Time instant $k$','Interpreter','latex')
ylabel('Measuremetn fault','Interpreter','latex')
%myprint

figure
plot(tk,f)
grid minor
hold on
plot(tk,festStaticCon)
xlabel('Time instant $k$ [-]','Interpreter','latex')
ylabel('Measurement fault','Interpreter','latex')
title('Measurement fault estimate (static residual with constraint)')


figure
plot(tk,f)
grid minor
hold on
plot(tk,festStaticUnc)
xlabel('Time instant $k$ [-]','Interpreter','latex')
ylabel('Measurement fault','Interpreter','latex')
title('Measurement fault estimate (static residual without constraint)')

figure
plot(tk,f)
grid minor
hold on
plot(tk,xaugf(nx+1:nx+ny,:))
xlabel('Time instant $k$ [-]','Interpreter','latex')
ylabel('Measurement fault','Interpreter','latex')
title('Measurement fault estimate (augmented model)')

figure('Name','ex2_staticresidual')
plot(tk,rstatic)
grid minor
xlabel('Time instant $k$ [-]','Interpreter','latex')
ylabel('Static residual $\mathbf{r}_{k}$','Interpreter','latex')


figure('Name','ex2_dynamicresidual')
plot(tk,rdynamic)
grid minor
xlabel('Time instant $k$ [-]','Interpreter','latex')
ylabel('Dynamic residual $\mathbf{r}_{k}$','Interpreter','latex')



figure
semilogy(tk,Sstatic)
hold on
plot([tk(1), tk(end)],[SstaticTh SstaticTh],'k')
grid minor
xlabel('Time instant $k$ [-]','Interpreter','latex')
ylabel('Statistics $S_{k}$','Interpreter','latex')

figure
semilogy(tk,Sdynamic)
hold on
plot([tk(1), tk(end)],[SdynamicTh SdynamicTh],'k')
grid minor
xlabel('Time instant $k$ [-]','Interpreter','latex')
ylabel('Statistics $S_{k}$','Interpreter','latex')


switch scenario
    case 1
        figure
        plot(tk,prst1)
        grid minor
        xlabel('Time instant $k$ [-]','Interpreter','latex')
        ylabel('Probability of fault-free model','Interpreter','latex')
    case 2
        figure
        plot(tk,festvsimm)
        grid minor
        xlabel('Time instant $k$ [-]','Interpreter','latex')
        ylabel('Measurement fault estimate','Interpreter','latex')
    otherwise
        error('unknown scenario')
end

figure('Name','nmodels-vsimm')
plot(tk,nModelvs)
grid minor
xlabel('Time instant $k$ [-]','Interpreter','latex')
ylabel('Number of models VS IMM','Interpreter','latex')
%myprint