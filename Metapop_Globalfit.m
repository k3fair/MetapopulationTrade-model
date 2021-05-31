%     Metapop_Globalfit.m runs global paramter fitting under low & high yield scenarios
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


function de_fit

sympref('HeavisideAtOrigin',0);
paramsfull = csvread('Meta_FittingTS_Global_Feb1.csv',1); %the 1 indicates we skip the header row introduced by writing in R

%Yield fitting, done independent of fitting for remaining variables

skips_yield=[];

    function C=basedynamics_fit(alpha, t)

        c0=paramsfull(1,9); %initial value of agricultural yield

        [T,Cv]=ode45(@DifEq,t,c0);

        function dC=DifEq(t,c)

            dcdt=zeros(1,1);
            dcdt(1)=alpha(1)*c(1)*(1-((c(1))/(alpha(2)))); %agricultural yield DE
            dC=dcdt;
        end
        C=Cv;

        for i = 1:length(skips_yield)
        C(skips_yield(i,2), skips_yield(i,1))=0;
        end

    end

t=paramsfull(1:53, 1);
c=paramsfull(1:53,9);
cmod=c;

       for i = 1:length(skips_yield)
        cmod(skips_yield(i,2), skips_yield(i,1))=0;
       end

%High yield scenario fitting for yield parameters
param0=[0.03,3.5]; %r, K initial guess values
lb = [0,2]; %r, K lower bounds
ub = [1,7]; %r, K upper bounds

[alpha,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat]=lsqcurvefit(@basedynamics_fit,param0,t,cmod,lb,ub); %cmod contains the yield data we're fitting to

r_yH=alpha(1); %stores these values separately for later use w. high yield fitting
K_yH=alpha(2);

%Low yield scenario fitting for yield parameters
param0=[0.03,3.5]; %r, K initial guess values
lb = [0,2]; %r, K lower bounds
ub = [1,3.5]; %r, K upper bounds

[alpha,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat]=lsqcurvefit(@basedynamics_fit,param0,t,cmod,lb,ub); %cmod contains the yield data we're fitting to

%For plotting results of yield fitting

    function C=basedynamics_display(alpha, t)

        c0=paramsfull(1,9); %initial value of agricultural yield

        [T,Cv]=ode45(@DifEq,t,c0);

        function dC=DifEq(t,c)

            dcdt=zeros(1,1);
            dcdt(1)=alpha(1)*c(1)*(1-((c(1))/(alpha(2)))); %agricultural yield DE
            dC=dcdt;
        end
        C=Cv;
    end

tv = linspace(min(t), max(t), max(t)-min(t)+1);
Cfit = basedynamics_display(alpha, tv);

set(groot,'defaultAxesColorOrder',[27 158 119; 217 95 2;]/255)

% %%%%%%%%%%%%%%%%%%%%%

%Fitting for remaining (non-yield) parameters under low yield scenario
close all

skips_main=[2 3 4];
skiptimes_main_start_ag=58;
skiptimes_main_start=54;

    function C=dynamics_fit(theta, t)

        c0=[paramsfull(1,2); paramsfull(1,6);
            paramsfull(1,3); paramsfull(1,9);]; %initial values for population, agricultural land area, food supply, and agricultural yield

        tau = paramsfull(55, 5)-paramsfull(55, 7);

        options = odeset('NonNegative', 1:4);
        [T,Cv]=ode45(@DifEq,t,c0, options);

        function dC=DifEq(t,c)

            rho = theta(1);
            sigma = 1/((theta(1)*theta(2)^2) - theta(2));
            alpha0=(theta(3)*theta(2)*exp(1/(theta(2)*theta(1) - 1)))/(theta(2)*theta(1) - 1);

            dcdt=zeros(4,1);

            dcdt(1)= c(1)*(alpha0*exp(-sigma*((c(3))/c(1)))*(rho-(c(1)/(c(3)))) - theta(4)); %Population DE
            dcdt(2)= theta(5)*((1/(1+ exp(theta(6)*(theta(7)- ((c(1))/(c(3))))))))*(tau - c(2)) - theta(8)*(c(2)); %Agricultural land DE
            dcdt(3)= (1-theta(9))*theta(10)*c(4)*c(2) - c(3); %Food supply DE
            dcdt(4) = theta(11)*c(4)*(1-((c(4))/(theta(12)))); %Agricultural yield DE

            dC=dcdt;

        end

        C=Cv;

        for i = 1:length(skips_main)
            if skips_main(i)==2 %exclude from fitting all agricultural land past 2017  (no data available)
                for j = skiptimes_main_start_ag:length(C(:,1))
                    C(j, skips_main(i))=0;
                end
            else %exclude from fitting all food and yield past 2013  (no data available)
                for j = skiptimes_main_start:length(C(:,1))
                    C(j, skips_main(i))=0;
                end
            end
        end

    end

t=paramsfull(1:end, 1);
c=horzcat(paramsfull(1:end,2),paramsfull(1:end,6),paramsfull(1:end,3), paramsfull(1:end,9));
cmod=c;

for i = 1:length(skips_main)
    if skips_main(i)==2 %exclude from fitting all agricultural land past 2017 (no data available)
        for j = (skiptimes_main_start_ag-1):length(cmod(:,1))
            cmod(j, skips_main(i))=0;
        end
    else %exclude from fitting all food and yield past 2013  (no data available)
        for j = (skiptimes_main_start-1):length(cmod(:,1))
            cmod(j, skips_main(i))=0;
        end
    end
end

%%%Current fit
p0 = [6.25,0.3673048,0.04556713,0.0113,0.003,10,0,0.001,0.069,0.75,alpha(1),alpha(2)]; %initial guesses for parameter values
lb = [6.25,0.3673048,0.04556713,0.0113,0,0,0,0.001,0.02,0.745,alpha(1),alpha(2)]; %lower bounds on parameters
ub = [6.25,0.3673048,0.04556713,0.0113,0.1,10,10,0.1,0.2,0.755,alpha(1),alpha(2)]; %upper bounds on parameters

options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',100000, 'MaxIterations',10000, 'FunctionTolerance', 1e-6);

[theta,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat]=lsqcurvefit(@dynamics_fit,p0,t,cmod,lb,ub,options); %cmod contains the population, agricultural land area, food supply, and yield data we're fitting to

fprintf(1,'\t Global parameters, low yield scenario:\n')
fprintf(1, '\t\t a0=%3.4f, sigma=%3.4f, delta=%3.4f, rho=%3.4f, \n', (theta(3)*theta(2)*exp(1/(theta(2)*theta(1) - 1)))/(theta(2)*theta(1) - 1), 1/((theta(1)*theta(2)^2) - theta(2)), theta(4), theta(1))
fprintf(1, '\t\t r_y=%3.4f, K_y=%3.4f \n', theta(11), theta(12))
fprintf(1, '\t\t f=%3.4f, gamma_A=%3.4f, beta_A=%3.4f, s=%3.4f, kappa=%3.4f, zeta=%3.4f  \n', theta(10), theta(6), theta(7), theta(9), theta(5), theta(8))

%To plot results from fitting

    function C=dynamics_display(theta, t)

               c0=[paramsfull(1,2); paramsfull(1,6);
            paramsfull(1,3); paramsfull(1,9);];  %initial values for population, agricultural land area, food supply, and agricultural yield

        tau = paramsfull(55, 5)-paramsfull(55, 7);

        options = odeset('NonNegative', 1:4);
        [T,Cv]=ode45(@DifEq,t,c0, options);

        function dC=DifEq(t,c)

            rho = theta(1);
            sigma = 1/((theta(1)*theta(2)^2) - theta(2));
            alpha0=(theta(3)*theta(2)*exp(1/(theta(2)*theta(1) - 1)))/(theta(2)*theta(1) - 1);

            dcdt=zeros(4,1);

            dcdt(1)= c(1)*(alpha0*exp(-sigma*((c(3))/c(1)))*(rho-(c(1)/(c(3)))) - theta(4)); %Population DE
            dcdt(2)= theta(5)*((1/(1+ exp(theta(6)*(theta(7)- ((c(1))/(c(3))))))))*(tau - c(2)) - theta(8)*(c(2)); %Agricultural land DE
            dcdt(3)= (1-theta(9))*theta(10)*c(4)*c(2) - c(3); %Food supply DE
            dcdt(4) = theta(11)*c(4)*(1-((c(4))/(theta(12)))); %Agricultural land yield DE

            dC=dcdt;

        end

        C=Cv;

    end

tv = linspace(min(t), 2100, 2100-min(t)+1);
Cfit = dynamics_display(theta, tv);

dlmwrite('Ky_35_params_NOurban_Feb5.csv', [theta tv(53) Cfit(53,:)] , 'delimiter', ',', 'precision', 10); %write file containing parameters for low yield scenario

cplot=c;              % make a copy of the data specifically for plotting
cplot(cplot==0)=nan;

%Save timeseries data for plotting low yield scenario
dlmwrite('GlobalTrajectoryMODEL_lowyield.csv',[transpose(tv), Cfit], 'delimiter', ',', 'precision', 5);
dlmwrite('GlobalTrajectoryFAO_lowyield.csv',[t, cplot], 'delimiter', ',', 'precision', 5);


fig1=figure(1);
subplot(2,3,1)
plot(t, cplot(:,1) , '.')
hold on
hlp = plot(tv, Cfit(:,1));
hold off
grid
xlim([min(tv) max(tv)])
xlabel('Year')
ylabel('P(t)')
title('Low yield scenario')

subplot(2,3,2)
plot(t, cplot(:,2) , '.');
hold on
hlp = plot(tv, Cfit(:,2));
plot(2030, cplot(40,2)+0.125, 'g.')
plot(2030, cplot(40,2)+0.416, 'g.')
plot(2030, cplot(40,2)+0.277, 'b.')
plot(2030, cplot(40,2)+0.168, 'r.')
plot(1998, 5.072, 'b.')
plot(2030, 5.349, 'b.')
plot(2050, cplot(50,2)+0.593, 'c.')
hold off
grid
xlim([min(tv) max(tv)])
xlabel('Year')
ylabel('A(t)')

subplot(2,3,3)
plot(t, cplot(:,3) , '.');
hold on
hlp = plot(tv, Cfit(:,3));
hold off
grid
xlim([min(tv) max(tv)])
xlabel('Year')
ylabel('F(t)')


subplot(2,3,4)
plot(t, cplot(:,4) , '.');
hold on
hlp = plot(tv, Cfit(:,4));
hold off
grid
xlim([min(tv) max(tv)])
xlabel('Year')
ylabel('y(t)')

x=0:0.01:2;

subplot(2,3,5)
hold on
plot(x, (1./(1+ exp(theta(6).*(theta(7)- (x))))))
hold off
grid
xlabel('P/F')
ylabel('b')

subplot(2,3,6)
hold on
plot(t, c(:,1)./c(:,3),'.')
plot(tv, (Cfit(:,1)./(Cfit(:,3))))
hold off
grid
xlim([min(tv) max(tv)])
ylabel('P/F')
xlabel('Year')

rho = theta(1);
sigma = 1/((theta(1)*theta(2)^2) - theta(2));
alpha0=(theta(3)*theta(2)*exp(1/(theta(2)*theta(1) - 1)))/(theta(2)*theta(1) - 1);

%%%% High yield fit

%Fitting for beta_A, gamma_A under high yield scenario

skips_main=[2 3 4];
skiptimes_main_start_ag=58;
skiptimes_main_start=54;

    function C=dynamics_fit_high(omega, t)

        c0=[paramsfull(1,2); paramsfull(1,6);
            paramsfull(1,3); paramsfull(1,9);];  %initial values for population, agricultural land area, food supply, and agricultural yield

        tau = paramsfull(55, 5)-paramsfull(55, 7);

        options = odeset('NonNegative', 1:4);
        [T,Cv]=ode45(@DifEq,t,c0, options);

        function dC=DifEq(t,c)

            rho = omega(1);
            sigma = 1/((omega(1)*omega(2)^2) - omega(2));
            alpha0=(omega(3)*omega(2)*exp(1/(omega(2)*omega(1) - 1)))/(omega(2)*omega(1) - 1);

            dcdt=zeros(4,1);

            dcdt(1)= c(1)*(alpha0*exp(-sigma*((c(3))/c(1)))*(rho-(c(1)/(c(3)))) - omega(4)); %Populaton DE
            dcdt(2)= omega(5)*((1/(1+ exp(omega(6)*(omega(7)- ((c(1))/(c(3))))))))*(tau - c(2)) - omega(8)*(c(2)); %Agricultural land area DE
            dcdt(3)= (1-omega(9))*omega(10)*c(4)*c(2) - c(3); %Food supply DE
            dcdt(4) = omega(11)*c(4)*(1-((c(4))/(omega(12)))); %Agricultural land yield DE

            dC=dcdt;

        end

        C=Cv;

        for i = 1:length(skips_main)
            if skips_main(i)==2 %exclude from fitting all agricultural land past 2017  (no data available)
                for j = skiptimes_main_start_ag:length(C(:,1))
                    C(j, skips_main(i))=0;
                end
            else %exclude from fitting all food and yield past 2013  (no data available)
                for j = skiptimes_main_start:length(C(:,1))
                    C(j, skips_main(i))=0;
                end
            end
        end

    end

t=paramsfull(1:end, 1);
c=horzcat(paramsfull(1:end,2),paramsfull(1:end,6),paramsfull(1:end,3), paramsfull(1:end,9));
cmod=c;

for i = 1:length(skips_main)
    if skips_main(i)==2 %exclude from fitting all agricultural land past 2017  (no data available)
        for j = (skiptimes_main_start_ag-1):length(cmod(:,1))
            cmod(j, skips_main(i))=0;
        end
    else %exclude from fitting all food and yield past 2013  (no data available)
        for j = (skiptimes_main_start-1):length(cmod(:,1))
            cmod(j, skips_main(i))=0;
        end
    end
end

%high yield fit for beta, gamma
p0 = [6.25,0.3673048,0.04556713,0.0113,theta(5),5,0.8,theta(8),theta(9),theta(10),r_yH,K_yH];  %initial guesses for parameter values
lb = [6.25,0.3673048,0.04556713,0.0113,theta(5),0,0,theta(8),theta(9),theta(10),r_yH,K_yH]; %lower bounds on parameter values
ub = [6.25,0.3673048,0.04556713,0.0113,theta(5),10,10,theta(8),theta(9),theta(10),r_yH,K_yH]; %upper bounds on parameter values

options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',100000, 'MaxIterations',10000, 'FunctionTolerance', 1e-6);

[omega,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat]=lsqcurvefit(@dynamics_fit_high,p0,t,cmod,lb,ub,options); %cmod contains the population, agricultural land area, food supply, and yield data we're fitting to


fprintf(1,'\t Global parameters, high yield scenario:\n')
fprintf(1, '\t\t a0=%3.4f, sigma=%3.4f, delta=%3.4f, rho=%3.4f, \n', (omega(3)*omega(2)*exp(1/(omega(2)*omega(1) - 1)))/(omega(2)*omega(1) - 1), 1/((omega(1)*omega(2)^2) - omega(2)), omega(4), omega(1))
fprintf(1, '\t\t r_y=%3.4f, K_y=%3.4f \n', omega(11), omega(12))
fprintf(1, '\t\t f=%3.4f, gamma_A=%3.4f, beta_A=%3.4f, s=%3.4f, kappa=%3.4f, zeta=%3.4f  \n', omega(10), omega(6), omega(7), omega(9), omega(5), omega(8))

%To plot results from fitting

    function C=dynamics_display_high(omega, t)

               c0=[paramsfull(1,2); paramsfull(1,6);
            paramsfull(1,3); paramsfull(1,9);];

        tau = paramsfull(55, 5)-paramsfull(55, 7);

        options = odeset('NonNegative', 1:4);
        [T,Cv]=ode45(@DifEq,t,c0, options);

        function dC=DifEq(t,c)

            rho = omega(1);
            sigma = 1/((omega(1)*omega(2)^2) - omega(2));
            alpha0=(omega(3)*omega(2)*exp(1/(omega(2)*omega(1) - 1)))/(omega(2)*omega(1) - 1);

            dcdt=zeros(4,1);

            dcdt(1)= c(1)*(alpha0*exp(-sigma*((c(3))/c(1)))*(rho-(c(1)/(c(3)))) - omega(4)); %population DE
            dcdt(2)= omega(5)*((1/(1+ exp(omega(6)*(omega(7)- ((c(1))/(c(3))))))))*(tau - c(2)) - omega(8)*(c(2)); %agricultural land area DE
            dcdt(3)= (1-omega(9))*omega(10)*c(4)*c(2) - c(3); %food supply DE
            dcdt(4) = omega(11)*c(4)*(1-((c(4))/(omega(12)))); %agricultural yield DE

            dC=dcdt;

        end

        C=Cv;

    end

%tv = linspace(min(t), max(t));
tv = linspace(min(t), 2100, 2100-min(t)+1);
Cfit = dynamics_display_high(omega, tv);

dlmwrite('Ky_70_params_NOurban_Feb5.csv', [omega tv(53) Cfit(53,:)] , 'delimiter', ',', 'precision', 10); %write file containing parameters for high yield scenario

cplot=c;              % make a copy of the data specifically for plotting
cplot(cplot==0)=nan;

%Save timeseries data for plotting high yield scenario
dlmwrite('GlobalTrajectoryMODEL_highyield.csv',[transpose(tv), Cfit], 'delimiter', ',', 'precision', 5);

fig2=figure(2);
subplot(2,3,1)
plot(t, cplot(:,1) , '.')
hold on
hlp = plot(tv, Cfit(:,1));
hold off
grid
xlim([min(tv) max(tv)])
xlabel('Year')
ylabel('P(t)')
title('High yield scenario')

subplot(2,3,2)
plot(t, cplot(:,2) , '.');
hold on
hlp = plot(tv, Cfit(:,2));
plot(2030, cplot(40,2)+0.125, 'g.')
plot(2030, cplot(40,2)+0.416, 'g.')
plot(2030, cplot(40,2)+0.277, 'b.')
plot(2030, cplot(40,2)+0.168, 'r.')
plot(1998, 5.072, 'b.')
plot(2030, 5.349, 'b.')
plot(2050, cplot(50,2)+0.593, 'c.')
hold off
grid
xlim([min(tv) max(tv)])
xlabel('Year')
ylabel('A(t)')

subplot(2,3,3)
plot(t, cplot(:,3) , '.');
hold on
hlp = plot(tv, Cfit(:,3));
hold off
grid
xlim([min(tv) max(tv)])
xlabel('Year')
ylabel('F(t)')


subplot(2,3,4)
plot(t, cplot(:,4) , '.');
hold on
hlp = plot(tv, Cfit(:,4));
hold off
grid
xlim([min(tv) max(tv)])
xlabel('Year')
ylabel('y(t)')

x=0:0.01:2;

subplot(2,3,5)
hold on
plot(x, (1./(1+ exp(omega(6).*(omega(7)- (x))))))
hold off
grid
xlabel('P/F')
ylabel('b')

subplot(2,3,6)
hold on
plot(t, c(:,1)./c(:,3),'.')
plot(tv, (Cfit(:,1)./(Cfit(:,3))))
hold off
grid
xlim([min(tv) max(tv)])
ylabel('P/F')
xlabel('Year')


end
