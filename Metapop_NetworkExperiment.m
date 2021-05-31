%     Metapop_NetworkExperiment.m runs simulations on networks with different densities and rewiring probabilities.
%     Copyright (C) 2021 Kathyrn R Fair
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

clear all;
close all;

sympref('HeavisideAtOrigin',0);

rng('shuffle','twister');

N=100; %%% Specify number of patches to split IC between
neighpick = 26; %(one of 17, 26, corresponding to network density of 0.343, 0.525 resp.)

eperc=16;
global eportion
eportion = eperc/100; %assign value for proportion of food available as export (mu)


%Read parameters in from fitting exercise
params = csvread('Ky_35_params_NOurban_Feb5.csv');

%Set global parameters
peak.loc=params(2);
peak.height=params(3);
global delta
delta = params(4);
global rho
rho = params(1);
global sigma
sigma = 1/((rho*peak.loc^2) - peak.loc);
global a0
a0=(peak.height*peak.loc*exp(1/(peak.loc*rho - 1)))/(peak.loc*rho - 1);
global K_y
K_y = params(12);
global r_y
r_y = params(11);
global h
h = 0.06;
global kappa
kappa = params(5);
global zeta
zeta = params(8);
global food
food = params(10);
global spoil
spoil = params(9);
global gamma_A
gamma_A = params(6);
global beta_A
beta_A = params(7);
global gamma_I
gamma_I = params(6);
global beta_I
beta_I = params(7);

global net

tspan = linspace(0, 100, 101); %specify points to evaluate at

%for Y, P, A, F respectively (assume 100 patches for now)
initdat = csvread('Meta_FittingTS_Global_Feb1.csv',1); %the 1 indicates we skip the header row introduced by writing in R
global T
T = initdat(55, 5)-initdat(55, 7);

init = [params(17); params(14)/N; params(15)/N; params(16)/N;];
for i = 2:N
    init = [init; params(17); params(14)/N; params(15)/N; params(16)/N;];
end

global T_patch
T_patch=T/N;

% ode45 calls incomplete/completnetCall with (t,x) which calls incomplete/completenetfcn with (t,x,N)
% Do this so we can specify # patches, in addition to our IC for the system
incompletenetCall = @(t,x) incompletenetfcn(t,x,N);

%Ensure all pop and resource levels are geq 0
options = odeset('NonNegative',1:4*N);
%, 'RelTol',1e-8,'AbsTol',1e-10

%Read in list of networks to use
fstruct = dir(sprintf('WSnet_1_%iN_%i_*.csv',N, neighpick));

nsims=length(fstruct);
storage_total=zeros(length(tspan)*N*nsims,15);

oops=0;
for loop1=1:nsims

net = csvread(fstruct(loop1).name,1); %the 1 indicates we skip the header row introduced by writing in R
file.string=strsplit(fstruct(loop1).name, '_');
neigh =str2double(file.string(4));
rewire=str2double(file.string(5));
realization=str2double(file.string(7));

    gamma_I = 7.5;
    beta_I = 0.25;

    gamma_A=gamma_I;
    beta_A=beta_I;

%Simulate ODE with incomplete network specified by net matrix
[Tic, Ric] = ode45(incompletenetCall, tspan, init, options);

if sum(sum(isnan(Ric)))==0 %excludes runs where variables go nan

    %save data on imports
Y = Ric(:,1:4:4*N);
P = Ric(:,2:4:4*N);
A = Ric(:,3:4:4*N);
F = Ric(:,4:4:4*N);

%recreate demand function values
b_A=zeros(length(tspan),N);
b_I=zeros(length(tspan),N);

for i=1:length(Tic)
    for j=1:N
    b_A(i,j) = 1/(1+ exp(gamma_A*(beta_A- ((P(i,j))/(F(i,j))))));
    b_I(i,j) = 1/(1+ exp(gamma_I*(beta_I- ((P(i,j))/(F(i,j))))));
    end
end

%recreate lambda values
lambda=zeros(length(tspan),N);

for i=1:length(Tic)
    for j=1:N

       if (sum(net(j,:)))
        lambda(i,j) = 1/sum(net(j,:)*transpose(b_I(i,:)));
    else
        lambda(i,j) = 0;
       end

    end
end

%recreate import values
imports=zeros(length(tspan),N);

for i=1:length(Tic)
    for j=1:N

    imports(i,j) = b_I(i,j)*(sum(food*(1-spoil)*(eportion).*transpose(net(:,j)).*Y(i,:).*A(i,:).*lambda(i,:)));

    end
end

node = 1;
%preallocate array size
storage_ts=zeros(length(tspan),15);
storage_run=zeros(length(tspan)*N,15);

counter_run=1;
for loop2=1:4:4*N

    for loop3=1:length(Tic)
            storage_ts(loop3,:) = [length(net(:,1)) neigh rewire gamma_I beta_I gamma_A beta_A node Tic(loop3) Ric(loop3,loop2) Ric(loop3,loop2+1) Ric(loop3,loop2+2) Ric(loop3,loop2+3) realization imports(loop3,node)];
    end
        rowadd=counter_run*length(tspan);
        storage_run((rowadd-(length(tspan)-1)):rowadd,:)=storage_ts;
        counter_run=counter_run+1;
        node=node+1;

end
        rowadd_tot=loop1*length(tspan)*N;
        storage_total((rowadd_tot-(length(tspan)*N-1)):rowadd_tot,:)=storage_run;

else
  oops=oops+1
end

if rem(loop1,10)==0
    loop1
end

end

dlmwrite(sprintf('PPlane_NOurban_betaA_gammaA_sWnet_N%i_neigh%i_WSnetexp_lowyield_export%i.csv',N, neigh, eperc),storage_total, 'delimiter', ',', 'precision', 5);

function f = incompletenetfcn(t, x, N)

global net
epsilon = 10^-4;

global a0
global delta
global sigma
global rho
global K_y
global r_y
global h
global kappa
global zeta
global food
global spoil
global gamma_A
global beta_A
global gamma_I
global beta_I
global T_patch
global eportion


Y = x(1:4:4*N);
P = x(2:4:4*N);
A = x(3:4:4*N);
F = x(4:4:4*N);

b_A(1) = 1/(1+ exp(gamma_A*(beta_A- ((P(1))/(F(1)))))); %agricultural land expansion demand sigmoid
b_I(1) = 1/(1+ exp(gamma_I*(beta_I- ((P(1))/(F(1)))))); %food import demand sigmoid

for j = 2:N
    b_A(j) = 1/(1+ exp(gamma_A*(beta_A- ((P(j))/(F(j))))));
    b_I(j) = 1/(1+ exp(gamma_I*(beta_I- ((P(j))/(F(j))))));
end

if sum(sum(net))>1 %for networked system
%Network is (i,j) (from,to)

for j = 1:N
    if (sum(net(j,:)))
        lambda(j) = 1/sum(net(j,:)*b_I(:));
    else
        lambda(j) = 0;
    end
end

else %for isolated nodes
    for j = 1:N
    lambda(j) = 0;
    end
end


f=[ Y(1)*r_y*(1-(Y(1)/K_y)); %yield DE
    P(1)*(a0*exp(-sigma*((F(1))/P(1)))*(rho-(P(1)/(F(1)))) - delta); % Population DE
    kappa*b_A(1)*(T_patch-A(1))-zeta*A(1); %Agricultural land area DE
    food*(1-spoil)*Y(1)*A(1)*(1 - lambda(1)*(eportion)*(net(1,:)*b_I(:))) + b_I(1)*(sum(food*(1-spoil)*eportion.*net(:,1).*Y(:).*A(:).*lambda(:))) - F(1);]; %Food supply DE

for k = 2:N
    f = [f; Y(k)*r_y*(1-(Y(k)/K_y));
        P(k)*(a0*exp(-sigma*((F(k))/P(k)))*(rho-(P(k)/(F(k)))) - delta);
        kappa*b_A(k)*(T_patch-A(k))-zeta*A(k);
        food*(1-spoil)*Y(k)*A(k)*(1 - lambda(k)*(eportion)*(net(k,:)*b_I(:))) + b_I(k)*(sum(food*(1-spoil)*eportion.*net(:,k).*Y(:).*A(:).*lambda(:))) - F(k);];

end

t;

end
