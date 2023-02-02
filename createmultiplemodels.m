function model = createmultiplemodels(A,B,G,q,C,H,f_i,p)
%CREATEMULTIPLEMODELS Creates mutliple model structure for
%example_faultdetection.m
%   model = createmultiplemodels(A,B,G,q,C,H,f_i,p) creates structure model
%   containing models with additive faults that are constructed using
%   dynamical matrix 'A', input marix 'B', state noise matrix 'G', state
%   noise mean value 'q', measurement matrix 'C', measurement noise matrix
%   'H', levels of aditive faults 'f_i', and single probability of
%   transition 'p'
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

parameterNames = {'f1', 'f2', 'f3','f4'};
parameterNumber = length(parameterNames);


nf = length(f_i);
for i = 1:parameterNumber
    % Values of fault component
    parameterValues.(parameterNames{i}) = num2cell(f_i); % parameter f_i
    % Probability of transition to neighboring fault magnitude
    parameterTransitionProbs.(parameterNames{i}) = sparse(diag([ 1-p (1-2*p)*ones(1,nf-2) 1-p]) + diag(p*ones(1,nf-1),1) + diag(p*ones(1,nf-1),-1));
end

[indices,combinationsTransProbabilities] = generateparametervariants(parameterTransitionProbs);
model.P = combinationsTransProbabilities;

modelParameterNames = {'A','B','G','Sigmaw','q','C','H','Sigmav','r'};

modelParameterNumber = length(modelParameterNames);

modelParameterFunctions = {...
    @(parameters)A;
    @(parameters)B;
    @(parameters)G;
    @(parameters)G*G';
    @(parameters)q;
    @(parameters)C;
    @(parameters)H;
    @(parameters)H*H';
    @(parameters)[parameters{1};parameters{2};parameters{3};parameters{4}]
    };

numberOfIndices = size(indices,1);

for modelNumber = 1:numberOfIndices
    % for every model variant select a unique parameter combination
    % specified by the indices matrix
    parameters = cell(parameterNumber,1);
    for parameter = 1:parameterNumber
        parameters{parameter} = parameterValues.(parameterNames{parameter}){indices(modelNumber,parameter)};
    end
    % for the given parameter set evaluate the model parameters
    for modelParameter = 1:modelParameterNumber
        model.M(modelNumber).(modelParameterNames{modelParameter}) = modelParameterFunctions{modelParameter}(parameters);
    end
end


end