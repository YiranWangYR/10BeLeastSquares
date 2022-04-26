% A basic code which generates deviated depth profile and compares linear 
% regression (Wang and Oskin, 202x) and a Bayesian approach for 10Be depth 
% profiles.

% This code generate hypothetical profiles based on given parameters.
% Imposed deviation and uncertainties are applied to mimic deviations in 
% realistic cases.

% Wang and Oskin, 202x

clear all;

%============'True' profile==========
Pn_0=10;
Pm1_0=0.014;
Pm2_0=0.0506;
density=2;
decay=4.997E-07;
Ln=160; %attenuation
Lm1=1500;
Lm2=4320;
La=[Ln,Lm1,Lm2];
P0=[Pn_0,Pm1_0,Pm2_0];

z_true=[25; 50; 70; 110; 150; 200]; 


inh_true=100000;  % atoms/g
D_true=1*Ln/density;  %cm
Time=200000; %yr
r_true=D_true/Time;

Pzn_true=Pn_0*exp(-density*z_true/Ln);
Pzm1_true=Pm1_0*exp(-density*z_true/Lm1);
Pzm2_true=Pm2_0*exp(-density*z_true/Lm2);

Rn_true=density*r_true/Ln+decay;
Rm1_true=density*r_true/Lm1+decay;
Rm2_true=density*r_true/Lm2+decay;

Ten_true=(1-exp(-Rn_true*Time))/Rn_true;
Tem1_true=(1-exp(-Rm1_true*Time))/Rm1_true;
Tem2_true=(1-exp(-Rm2_true*Time))/Rm2_true;

% true concentration
Cn_true=Ten_true*Pzn_true;
Cm_true=Tem1_true*Pzm1_true+Tem2_true*Pzm2_true;
C_true=Cn_true+Cm_true+inh_true;


%=============Introduce uncertainty+noise===========
C_sd=0.02; %analytical uncertainty of sample concentration; percentage
C_err=0.05; %imposed deviation; percentage
z_sd=5;  %sample depth uncertainty

D_est=D_true; % Estimated denudation
D_sd=10;    % uncertainty of estimated denudation


N=length(C_true); % sample size
C_mean=round(C_true.*(1+C_err*randn(N,1)));  % sample concentration (deviated)


%==============Linear inversion (Denudation-depth approach)===========
P=50000; %numbers of repetition for the Monte Carlo simulation of linear inversion 

%---------MC sampling
% Concentration
for i=1:P
        y(:,i)=C_mean.*(1+(C_sd.*randn(N,1)));
end
% Sample depth
for i=1:P
        z_measure(:,i)=z_true-z_sd+2*z_sd.*rand(N,1);
end

% denudation
D_measure(:,1)=D_est+D_sd*randn(P,1);


%----------inversion
G2=zeros(N,2);
for i=1:P
    Pzn=Pn_0*exp(-density*z_measure(:,i)/Ln);
    Pzm1=Pm1_0*exp(-density*z_measure(:,i)/Lm1);
    Pzm2=Pm2_0*exp(-density*z_measure(:,i)/Lm2);
    if D_measure(i)==0  % no erosin case
        fm1=1;
        fm2=1;
    else     % erosion exist
        fm1=exp(-0.5*(density*D_measure(i)/Lm1-density*D_measure(i)/Ln)+(1/24)*((density*D_measure(i)/Lm1)^2-(density*D_measure(i)/Ln)^2));
        fm2=exp(-0.5*(density*D_measure(i)/Lm2-density*D_measure(i)/Ln)+(1/24)*((density*D_measure(i)/Lm2)^2-(density*D_measure(i)/Ln)^2));
    end
    x=Pzn+Pzm1*fm1+Pzm2*fm2;
    G2(:,1)=1;
    G2(:,2)=x;
    [Up, Lp, Vp] = svd( G2);
    B = Vp/Lp*Up';
    mest2=B*y(:,i);
    inh_est2(i,1)=mest2(1,1);
    Ten_est2(i,1)=mest2(2,1);
    if D_measure==0
         t_est2(i,1)=-log(1-Ten_est2(i,1)*decay)/decay;      % actual exposure age when there's no erosion, yr
    else
        t_est2(i,1)=Be10Newton(Ten_est2(i,1),Ten_est2(i,1),D_measure(i),100,Ln,density,decay); 
    end
end

t_est2_ci(1,:)=prctile(t_est2(1:P,1),[2.5 97.5]);
t_est2_mean=mean(t_est2(1:P,1));
Cinh_est2_ci(1,:)=prctile(inh_est2(1:P,1),[2.5 97.5]);
Cinh_est2_mean=mean(inh_est2(1:P,1));



%========bayesian (MCMC scheme)================
%prior
a=20000;    %lower boundary of age
b=1000000;   %upper boundary of age
c=000000;  %lower boundary of inh
d=260000;   %upper boundary of inh


t0=150000;
D0=D_est;
Cinh0=80000;
z0=z_true;

if (t0>=a) && (t0<=b)
    p_t=1/(b-a);
else
    p_t=0;
end

p_D=1/(D_sd*sqrt(2*pi))*exp(-0.5*((D0-D_est)/D_sd)^2);

if (Cinh0>=c)&&(Cinh0<=d)
    p_inh=1/(d-c);
else
    p_inh=0;
end

for i=1:N
    if (z0(i,1)>=(z_true(i,1)-z_sd))&&(z0(i,1)<=(z_true(i,1)+z_sd))
        p_z(i,1)=1/(2*z_sd);
    else
        p_z(i,1)=0;
    end
end
PZ=prod(p_z);

p0=log(p_t*p_D*p_inh*PZ);
f0=modelconcentration(P0,z0,La,decay,density,D0/t0,t0,Cinh0);
L0=loglike(C_mean,f0,C_mean*C_sd);



% number of steps
k=300000;

% step size
t_step=10000;
D_step=1;
C_step=10000;
z_step=1;

t=zeros(k,1);
D=zeros(k,1);
Cinh=zeros(k,1);
z=zeros(k,N);
p=zeros(k,1);
L=zeros(k,1);
t(1,1)=t0;
D(1,1)=D0;
Cinh(1,1)=Cinh0;
z(1,:)=z0';
p(1,1)=p0;
L(1,1)=L0;
t_rej=zeros(k,1);
D_rej=zeros(k,1);
Cinh_rej=zeros(k,1);
z_rej=zeros(k,N);
acceptance=0;



% M-H MCMC
for i=2:k
    T1=t(i-1,1)+t_step*(rand-0.5);
    D1=D(i-1,1)+D_step*(rand-0.5);
    Cin1=Cinh(i-1,1)+C_step*(rand-0.5);
    Z1=z(i-1,:)+z_step*(rand(1,N)-0.5);
    Z1=Z1';
    if (T1>=a) && (T1<=b)
        p_t1=1/(b-a);
    else
        p_t1=0;
    end
    
    p_D1=1/(D_sd*sqrt(2*pi))*exp(-0.5*((D1-D_est)/D_sd)^2);
    
    for j=1:N
        if (Z1(j,1)>=(z_true(j,1)-z_sd))&&(Z1(j,1)<=(z_true(j,1)+z_sd))
        p_z1(j,1)=1/(2*z_sd);
        else
        p_z1(j,1)=0;
        end
    end
    PZ1=prod(p_z1);
    
    if (Cin1>=c)&&(Cin1<=d)
        p_inh1=1/(d-c);
    else
        p_inh1=0;
    end
    
    p1=log(p_t1*p_D1*p_inh1*PZ1);
    
    f1=modelconcentration(P0,Z1,La,decay,density,D1/T1,T1,Cin1);
    L1=loglike(C_mean,f1,C_mean*C_sd);
    mu=rand(1);
    if exp(p1+L1-p(i-1,1)-L(i-1,1))>mu
        t(i,1)=T1;
        D(i,1)=D1;
        Cinh(i,1)=Cin1;
        z(i,:)=Z1';
        p(i,1)=p1;
        L(i,1)=L1;
        acceptance=acceptance+1;
    else
        t(i,1)=t(i-1,1);
        D(i,1)=D(i-1,1);
        Cinh(i,1)=Cinh(i-1,1);
        z(i,:)=z(i-1,:);
        p(i,1)=p(i-1,1);
        L(i,1)=L(i-1,1);
        t_rej(i,1)=T1;
        D_rej(i,1)=D1;
        Cinh_rej(i,1)=Cin1;
        z_rej(i,:)=Z1';
    end
end

acc_ratio=acceptance/k;

t_est1_ci(1,:)=prctile(t(5000:k,1),[2.5 97.5]);
t_est1_mean=mean(t(5000:k,1));
Cinh_est1_ci(1,:)=prctile(Cinh(5000:k,1),[2.5 97.5]);
Cinh_est1_mean=mean(Cinh(5000:k,1));

%{
%--------------------figure---------------
figure(1)
subplot(2,2,1);
histogram(t(5000:k,1));
title('Age-Bayesian');

subplot(2,2,2);
histogram(t_est2);
title('Age-Linear regression');

subplot(2,2,3)
histogram(Cinh(5000:k,1));
title('Inheritance-Bayesian');

subplot(2,2,4)
histogram(inh_est2);
title('Inheritance-Linear regression');


% depth profile fit
C_min=C_mean-C_sd;
C_max=C_mean+C_sd;

K=500;

zmodel=(1:max(z_true)+10)';
figure (2)

% LR
% best fit lines. P vs C
x1=zeros(P,1);
x2=Pn_0*ones(P,1);
y1=inh_est2;
TeN=Ten_est2;
y2=TeN.*x2+y1;

% error of production rate at depth
Pz_err=sqrt(0+(density*z_sd*ones(N,1)/Ln).^2);

subplot(2,2,1)
%axis([0 Pn_0 0 max(y2)]);
xlabel('Production rate at depth (Pz; atome*g^{-1}*yr^{-1})')
ylabel('^{10}Be concentration (atoms/g)')
hold on
%============line fit================
for j=1:K
    i=round((P-1)*rand(1))+1;
   plot([x1(i),x2(i)],[y1(i),y2(i)],'color',[0.5 0.5 0.5])
end

%data
errorbar(x,C_mean,C_sd*C_mean,C_sd*C_mean,Pz_err,Pz_err,'bo','linewidth',1.5);%vertical errorbar
plot(Pzn_true,C_true,'rx','linewidth',1.5);

set(gca,'FontSize',20)
set(gca,'linewidth',2)
title('Linear Regression', 'Fontsize',20,'fontweight','bold');
hold off


% depth profile
for j=1:K
    i=round((P-1)*rand(1))+1;
    Pzn=Pn_0*exp(-density*zmodel/Ln);
    Pzm1=Pm1_0*exp(-density*zmodel/Lm1);
    Pzm2=Pm2_0*exp(-density*zmodel/Lm2);
    B=density*D_measure(i)./[Ln, Lm1, Lm2];
    Te_sim=(1-exp(-B-decay*t_est2(i)))./(B/(t_est2(i))+decay);
    fmodel(:,j)=Pzn*Te_sim(1)+Pzm1*Te_sim(2)+Pzm2*Te_sim(3)+Cinh(i);
end
subplot(2,2,2)
%axis([0 round(max(C_mean)*1.1, -4) -(max(z_true)+10) 0])
xlabel('^{10}Be concentration (atoms/g)')
ylabel('Depth (m)')
hold on
for i=1:K
    plot(fmodel(:,i),-zmodel,'Color', [0.5 0.5 0.5]);
end

errorbar(C_mean,-z_true,z_sd*ones(N,1),z_sd*ones(N,1),C_sd*C_mean,C_sd*C_mean,'bo','linewidth',1.5); %vertical errorbar
plot(C_true,-z_true,'rx','linewidth',1.5);
set(gca,'FontSize',20)
set(gca,'linewidth',2)
title('Linear Regression', 'Fontsize',20,'fontweight','bold');
hold off


%Bayesian
% best fit lines. P vs C
x1=zeros(k,1);
x2=Pn_0*ones(k,1);
y1=Cinh;
Dn=density*D/Ln;
TeN_B=(1-exp(-Dn-decay*t))./(Dn./t+decay);
y2=TeN_B.*x2+y1;

% error of production rate at depth
Pz_err=sqrt(0+(density*z_sd*ones(N,1)/Ln).^2);

subplot(2,2,3)
%axis([0 Pn_0 0 max(y2)]);
xlabel('Production rate at depth (Pz; atome*g^{-1}*yr^{-1})')
ylabel('^{10}Be concentration (atoms/g)')
hold on
%============line fit================
for j=1:K
    i=round((k-5001)*rand(1))+5001;
   plot([x1(i),x2(i)],[y1(i),y2(i)],'color',[0.5 0.5 0.5])
end

%data
errorbar(x,C_mean,C_sd*C_mean,C_sd*C_mean,Pz_err,Pz_err,'bo','linewidth',1.5);%vertical errorbar
plot(Pzn_true,C_true,'rx','linewidth',1.5);

set(gca,'FontSize',20)
set(gca,'linewidth',2)
title('Bayesian', 'Fontsize',20,'fontweight','bold');
hold off


% depth profile


for j=1:K
    i=round((k-5001)*rand(1))+5001;
    Pzn=Pn_0*exp(-density*zmodel/Ln);
    Pzm1=Pm1_0*exp(-density*zmodel/Lm1);
    Pzm2=Pm2_0*exp(-density*zmodel/Lm2);
    B=density*D(i)./[Ln, Lm1, Lm2];
    Te_sim=(1-exp(-B-decay*t(i)))./(B/(t(i))+decay);
    fmodel(:,j)=Pzn*Te_sim(1)+Pzm1*Te_sim(2)+Pzm2*Te_sim(3)+Cinh(i);
end

subplot(2,2,4)
%axis([0 round(max(C_mean)*1.1, -4) -(max(z_true)+10) 0])
xlabel('^{10}Be concentration (atoms/g)')
ylabel('Depth (m)')
hold on
for i=1:K
    plot(fmodel(:,i),-zmodel,'Color', [0.5 0.5 0.5]);
end

errorbar(C_mean,-z_true,z_sd*ones(N,1),z_sd*ones(N,1),C_sd*C_mean,C_sd*C_mean,'bo','linewidth',1.5); %vertical errorbar
plot(C_true,-z_true,'rx','linewidth',1.5);

set(gca,'FontSize',20)
set(gca,'linewidth',2)
title('Bayesian', 'Fontsize',20,'fontweight','bold');

hold off


% display
fprintf('t-Bayesian \n');
fprintf('95CI: [%d, %d] \n', round(t_est1_ci(1,1)), round(t_est1_ci(1,2)));
fprintf('mean: %d \n', round(t_est1_mean));

fprintf('Inheritance-Bayesian \n');
fprintf('95CI: [%d, %d]\n', round(Cinh_est1_ci(1,1)), round(Cinh_est1_ci(1,2)));
fprintf('mean: %d \n', round(Cinh_est1_mean));

fprintf('t-LS \n');
fprintf('95CI: [%d, %d] \n', round(t_est2_ci(1,1)), round(t_est2_ci(1,2)));
fprintf('mean: %d \n', round(t_est2_mean));

fprintf('Inheritance-LS \n');
fprintf('95CI: [%d, %d]\n', round(Cinh_est2_ci(1,1)), round(Cinh_est2_ci(1,2)));
fprintf('mean: %d \n', round(Cinh_est2_mean));
%}

Re_LS=[t_est2_mean,t_est2_ci,Cinh_est2_mean,Cinh_est2_ci];
Re_Bayesian=[t_est1_mean,t_est1_ci,Cinh_est1_mean,Cinh_est1_ci];


function f=modelconcentration(P0,depth,La,decay,density,erosion,t,Cinh)
Pzn=P0(1)*exp(-density.*depth/La(1));
Pzm1=P0(2)*exp(-density.*depth/La(2));
Pzm2=P0(3)*exp(-density.*depth/La(3));
B=density*erosion./La+decay;
Te=(1-exp(-B.*t))./B;
f=Pzn.*Te(:,1)+Pzm1.*Te(:,2)+Pzm2.*Te(:,3)+Cinh;
end

function L=loglike(d,f,sig)
% natural log of the likelihood of a certain group of concentrations
% different sd for each sample
% d: measured data
% f: modeled 
esquare=(d-f).^2;
k=length(d);
for i=1:k
    li(i,1)=log((1/(sqrt(2*pi)*sig(i,1)))*exp(-(1/(2*sig(i,1)^2))*esquare(i,1)));
end
L=sum(li);
end
