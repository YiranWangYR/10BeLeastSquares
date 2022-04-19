% version 2.0
% Code to apply a least squares linear regression method of exposure age 
% estimation from cosmogenic nuclide depth-profiles. 
% Detailed explaination can be find in: Wang and Oskin (202x)
% This code estimate age with eroded thickness(denudation depth) D, muogenic production included.
% A Monte Carlo simulation is applied to find the distribution of age and
% inheritance. Erosion may be zero.

% Update notes on 2.0: 
% 1. This version is tested on matlab2020, and may not work on older versions.
% Version 1.0 is compatible with matlab 2014, no essential difference
% between these two versions.
% 2. The Newton method function has been integrated into this script.
% 3. Small revisions on the figures. PDFs were added to the plots.

% If you encounter any problem, please contact: yrwwang@ucdavis.edu


% Please edit the values betwenn line 24-68, following the instructions
% below

clear all;

%=============================Inputs==============================
% -------Main variables---------
% !! Format !!
% variable=[x,y,z]
% x=the mean of the variable;
% y=the standard diviation of the variable (normal distribution), or the
% error (normal distribution),y may equal to 0;
% z=type of distribution (1=normal distribution; 0=uniform distribution).
P0nD=[20,0,1]; %surface production rate for nucleon spallation, atoms/(g*yr)
P0m1D=[0.31,0,1]; %surface production rate for negative muons, atoms/(g*yr)
P0m2D=[0.13,0,1]; %surface production rate for fast muons, atoms/(g*yr)
densityD=[2.2,0.2,1]; %g/cm3, sediment density
LanD=[160,0,1]; % nucleon spallation attenuation length; g/cm2
Lam1D=[1500,0,1]; %negative muon attenuation length; g/cm2
Lam2D=[4300,0,1]; %fast muon attenuation length
ErodedD=[40, 10, 1]; % cm, eroded thickness. No erosion when the first and second terms ==0
% When there is no erosion, this code returns results with muogenic
% production fully incorperated, and no approximation. (Eq.3-5 of Wang and Oskin (202x)
%-----------

% load sample data. Change the file name accordingly
D1=load('../Be10 code/data/T2_Fin.txt');
% Data format:
% 1st column: sample depth, cm;
% 2nd column: standard deviation of sample depth, cm;
% 3rd column: sample concentration, atoms/g
% 4th column: standard deviation of C, atoms/g (!!!put in actual value here, 
% do not use the percentage format)

% distribution type of samples
z_dis=1;    %distribution type of sample depth. 1=normal distribution; 0=uniform distribution
y_dis=1;    %distribution type of sample concentration. 1=normal distribution; 0=uniform distribution

% decay constant
decay=0.0000004997;    %decay constant

P=5000; % Number of iterations for the Monte Carlo simulation

Np=100;  % Number of iterations for the Newton's method to find exposure age. Default is 100.


% output excel file name
filename=(['10Be_thickness.xlsx']);



%=====================Input ends=====================
%===============Do not change following codes=====================

%==============Data preperation===============
z_mean=D1(:,1); % cm, sample depth (from gravel top)
z_sd=D1(:,2); % cm, standard deviation of z
y_mean=D1(:,3); % sample concentration, atoms/g
y_sd=D1(:,4); % standard deviation of C

N=length(y_mean); % number of samples


%============Stope if sd of C is in percentage format===========
if y_sd(1)<1
    fprintf('Please use actual value of the standard deviation of concentration. Do not use the percentage format.');
    return;
end

%-------------Sampling data (Monte Carlo)--------------
%P0n
if P0nD(3)==1
    P0n=P0nD(1)+P0nD(2)*randn(P,1);
else P0n=P0nD(1)+2*P0nD(2)*(rand(P,1)-0.5);
end

%P0m1
if P0m1D(3)==1
    P0m1=P0m1D(1)+P0m1D(2)*randn(P,1);
else P0m1=P0m1D(1)+2*P0m1D(2)*(rand(P,1)-0.5);
end

%P0m2
if P0m2D(3)==1
    P0m2=P0m2D(1)+P0m2D(2)*randn(P,1);
else P0m2=P0m2D(1)+2*P0m2D(2)*(rand(P,1)-0.5);
end

%density
if densityD(3)==1
    density=densityD(1)+densityD(2)*randn(P,1);
else density=densityD(1)+2*densityD(2)*(rand(P,1)-0.5);
end

%attenuation
if LanD(3)==1
    Lan=LanD(1)+LanD(2)*randn(P,1);
else Lan=LanD(1)+2*LanD(2)*(rand(P,1)-0.5);
end

if Lam1D(3)==1
    Lam1=Lam1D(1)+Lam1D(2)*randn(P,1);
else Lam1=Lam1D(1)+2*Lam1D(2)*(rand(P,1)-0.5);
end

if Lam2D(3)==1
    Lam2=Lam2D(1)+Lam2D(2)*randn(P,1);
else Lam2=Lam2D(1)+2*Lam2D(2)*(rand(P,1)-0.5);
end


%eroded thickness
if ErodedD(3)==1
    D=ErodedD(1)+ErodedD(2)*randn(P,1);
else D=ErodedD(1)+2*ErodedD(2)*(rand(P,1)-0.5);
end

%sample concentration
if y_dis==1
    for i=1:P
        y(:,i)=y_mean+round(y_sd.*randn(N,1));
    end
else
    for i=1:P
        y(:,i)=y_mean+2*y_sd.*(rand(N,1)-0.5);
    end
end

%sample depth
if z_dis==1
    for i=1:P
        z(:,i)=z_mean+z_sd.*randn(N,1);
    end
else
    for i=1:P
        z(:,i)=z_mean+2*z_sd.*(rand(N,1)-0.5);
    end
end

%---------Least Squares inversion---------
M=2; % number of unknowns
G=zeros(N,M);

Re=zeros(P,6);  % results to be saved here

for i=1:P
    La=[Lan(i),Lam1(i),Lam2(i)];
    % calculate production rate at each sample depth
    Pzn=P0n(i)*exp(-density(i)*z(:,i)/La(1));
    Pzm1=P0m1(i)*exp(-density(i)*z(:,i)/La(2));
    Pzm2=P0m2(i)*exp(-density(i)*z(:,i)/La(3));
    if D(i)==0  % no erosion case
        gm1=1;
        gm2=1;
    else     % erosion exists, calculate g based on eq. 9
        gm1=exp(-0.5*(density(i)*D(i)/La(2)-density(i)*D(i)/La(1))+(1/24)*((density(i)*D(i)/La(2))^2-(density(i)*D(i)/La(1))^2));
        gm2=exp(-0.5*(density(i)*D(i)/La(3)-density(i)*D(i)/La(1))+(1/24)*((density(i)*D(i)/La(3))^2-(density(i)*D(i)/La(1))^2));
    end
    x=Pzn+Pzm1*gm1+Pzm2*gm2;    % effective production rate at sample depth
    % build matrix for the linear equation. x is production rate at depth,
    % y is sample concentration; eq. 4 and 10
    G(:,1)=1;
    G(:,2)=x;
    [Up, Lp, Vp] = svd(G);  % singular value decomposition
    Inv = Vp/Lp*Up';
    mest = Inv*y(:,i);  % inversion results
    Te=mest(2,1);  % effective exposure age of nucleon (age without decay)
    C_inh=mest(1,1);  % inherited concentration
    if D==0
         t=-log(1-Te*decay)/decay;      % actual exposure age when there's no erosion, yr; eq. 5
    else
        t=Be10NewtonN(Te,Te,D(i),Np,La(1),density(i),decay);    % actual exposure age when erosion exist; eqs. 11-13
    % 1st variable: Te, effective age
    % 2nd variable: x, first guess of the age, default as Te
    % 3rd variable: D, eroded thickness
    % 4th variable: n, number of iterations, default as 100
    % 5th variable: attenuation length of nucleon
    % 6th variable: density, sample sediment density
    % 7th variable: decay constant
    end
    Cinh_o = C_inh/exp(-decay*t);  % the origional inheritance corrected for decay
    r=D(i)/t;      % erosion rate
    Re(i,:)=[Te, C_inh, t, Cinh_o, r, D(i)];
end

%==========================Save results====================
%xlswrite(filename,Re); %results saved as excel, please change file name.
%first column: Te value (effective exposure age,yr)
%second column: Cinh (inhereted concentration in samples)
%third column: t(exposure age before loess),
%forth column: Cinh_o(original inhereted concentration, corrected for decay),
%fifth column: r(erosion rate)
%sixth column: D(total erosion)



%==========================Display results=================
Te=Re(:,1)/1000; %effective exposure age in kyr
C_inh=Re(:,2); % atoms/g
t=Re(:,3)/1000; %exposure age in kyr
Cinh_o=Re(:,4); % atoms/g
r=Re(:,5)*1000;  %erosion rate in cm/kyr
De=Re(:,6); %eroded thickness, cm

%t
t_mean=mean(t);
t_median=median(t);
t_ci(1,:)=prctile(t,[2.5 97.5]);
t_ci(2,:)=prctile(t,[18 82]);
%Te
Te_mean=mean(Te);
Te_median=median(Te);
Te_ci(1,:)=prctile(Te,[2.5 97.5]);
Te_ci(2,:)=prctile(Te,[18 82]);
%C_inh
Cinh_mean=mean(C_inh);
Cinh_median=median(C_inh);
Cinh_ci(1,:)=prctile(C_inh,[2.5 97.5]);
Cinh_ci(2,:)=prctile(C_inh,[18 82]);
%r
r_mean=mean(r);
r_median=median(r);
r_ci=prctile(r,[2.5 97.5]);
%De
De_mean=mean(De);
De_median=median(De);
De_ci=prctile(De,[2.5 97.5]);


%==============display============
fprintf('Exposure age (kyr) \n');
fprintf('2-sigma CI: [%.3f, %.3f] \n', round(t_ci(1,1),3), round(t_ci(1,2),3));
fprintf('1-sigma CI: [%.3f, %.3f] \n', round(t_ci(2,1),3), round(t_ci(2,2),3));
fprintf('Mean: [%.3f] \n', round(t_mean,3));
fprintf('Median: [%.3f] \n', round(t_median,3));

fprintf('Inheritance (atoms/g) \n');
fprintf('2-sigma CI: [%d, %d]\n', round(Cinh_ci(1,1)), round(Cinh_ci(1,2)));
fprintf('1-sigma CI: [%d, %d]\n', round(Cinh_ci(2,1)), round(Cinh_ci(2,2)));
fprintf('Mean: [%d] \n', round(Cinh_mean));
fprintf('Median: [%d] \n', round(Cinh_median));



%---------------Figure 1 inverted results vs data----------------
K=500;   %number of simulation

% best fit lines
x1=zeros(P,1);
x2=P0nD(1)*ones(P,1);
y1=Re(:,2);
TeN=Re(:,1);
y2=TeN.*x2+y1;

% error of production rate at depth
Pz_err=sqrt(P0nD(2)^2+(densityD(1)*z_sd/LanD(1)).^2+(densityD(2)*z_mean/LanD(1)).^2+(LanD(2)*densityD(1)*z_mean/LanD(1)^2).^2).*(P0nD(1)*exp(-densityD(1)*z_mean/LanD(1)));


figure(1)
subplot(2,2,1)
axis([0 P0nD(1) 0 max(y2)]);
xlabel('Production rate at depth (Pz; atome*g^{-1}*yr^{-1})')
ylabel('^{10}Be concentration (atoms/g)')
hold on
%============line fit================
for j=1:K
    i=round((P-1)*rand(1))+1;
   plot([x1(i),x2(i)],[y1(i),y2(i)],'color',[0.5 0.5 0.5])
end

%data
errorbar(x,y_mean,y_sd,y_sd,Pz_err,Pz_err,'bo','linewidth',1.5);%vertical errorbar
%{
y_lim=get(gca,'ylim');
y_l=y_lim(1,2)/100;

%this plot the horizontal errorbar. Newer versions of matlab can simplify
%this step
for i=1:6
    plot([x(i)-Pz_err(i),x(i)+Pz_err(i)],[y_mean(i),y_mean(i)],'b','linewidth',1.5);
    plot([x(i)-Pz_err(i),x(i)-Pz_err(i)],[y_mean(i)-y_l,y_mean(i)+y_l],'b','linewidth',1.5);
    plot([x(i)+Pz_err(i),x(i)+Pz_err(i)],[y_mean(i)-y_l,y_mean(i)+y_l],'b','linewidth',1.5);
end
%}
set(gca,'FontSize',20)
set(gca,'linewidth',2)
hold off



% depth profile fit
y_min=y_mean-y_sd;
y_max=y_mean+y_sd;

zmodel=(1:max(z_mean)+10)';
for j=1:K
    i=round((P-1)*rand(1))+1;
    Pzn=P0n(i)*exp(-density(i)*zmodel/Lan(i));
    Pzm1=P0m1(i)*exp(-density(i)*zmodel/Lam1(i));
    Pzm2=P0m2(i)*exp(-density(i)*zmodel/Lam2(i));
    B=density(i)*De(i)./[Lan(i), Lam1(i), Lam2(i)];
    Te_sim=(1-exp(-B-decay*t(i)*1000))./(B/(t(i)*1000)+decay);
    fmodel(:,j)=Pzn*Te_sim(1)+Pzm1*Te_sim(2)+Pzm2*Te_sim(3)+C_inh(i);
end

subplot(2,2,2)
axis([0 round(max(y_mean)*1.1, -4) -(max(z_mean)+10) 0])
xlabel('^{10}Be concentration (atoms/g)')
ylabel('Depth (m)')
hold on
for i=1:K
    plot(fmodel(:,i),-zmodel,'Color', [0.5 0.5 0.5]);
end
%{
%this plot the horizontal errorbar. Newer versions of matlab can simplify
%this step
for i=1:N
    plot([y_min(i,1),y_max(i,1)],[-z_mean(i,1),-z_mean(i,1)],'-b','linewidth',1.5);
    plot([y_min(i,1),y_min(i,1)],[-z_mean(i,1)-2,-z_mean(i,1)+2],'-b','linewidth',1.5);
    plot([y_max(i,1),y_max(i,1)],[-z_mean(i,1)-2,-z_mean(i,1)+2],'-b','linewidth',1.5);
end
%}
errorbar(y_mean,-z_mean,z_sd,z_sd,y_sd,y_sd,'bo','linewidth',1.5); %vertical errorbar
set(gca,'FontSize',20)
set(gca,'linewidth',2)

hold off


%histogram of Te
subplot(2,2,3);
histogram(Te,'Normalization','pdf');
xlabel('T_e value (kyr)');
ylabel('Probability');
%Te_y=get(gca,'ylim');
hold on
[f_Te,xi_Te]=ksdensity(Te);
plot(xi_Te,f_Te,'m','LineWidth',2);
xline(Te_ci(1,1),'r','LineWidth',2);
xline(Te_ci(1,2),'r','LineWidth',2);
xline(Te_mean,'b','LineWidth',2);
xline(Te_median,'g','LineWidth',2);
set(gca,'box','off')
set(gca,'FontSize',20)
set(gca,'linewidth',2)
hold off

%histogram of C_inh
subplot(2,2,4);
histogram(C_inh,'Normalization','pdf');
hold on
[f_inh,xi_inh]=ksdensity(C_inh);
plot(xi_inh,f_inh,'m','LineWidth',2);
xlabel({'Inherited concentration (atoms/g)'});
ylabel('Probability');
%Cinh_y=get(gca,'ylim');
xline(Cinh_ci(1,1),'r','LineWidth',2);
xline(Cinh_ci(1,2),'r','LineWidth',2);
xline(Cinh_mean,'b','LineWidth',2);
xline(Cinh_median,'g','LineWidth',2);
set(gca,'box','off')
set(gca,'FontSize',20)
set(gca,'linewidth',2)
hold off


%==============plot data===========

figure(2)
subplot(2,2,1:2)
% histogram of t1
histogram(t,'Normalization','pdf'); 
xlabel('Exposure age (kyr)');
ylabel('Probability');
%t1_y=get(gca,'ylim');
hold on
[f_t,xi_t]=ksdensity(t);
plot(xi_t,f_t,'m','LineWidth',2);
xline(t_ci(1,1),'r','LineWidth',2);
xline(t_ci(1,2),'r','LineWidth',2);
xline(t_mean,'b','LineWidth',2);
xline(t_median,'g','LineWidth',2);
set(gca,'box','off')
set(gca,'FontSize',20)
set(gca,'linewidth',2)
hold off


%histogram of r
subplot(2,2,3);
histogram(r,'Normalization','pdf');
xlabel('Denudation rate (cm/kyr)');
ylabel('Probability');
%r_y=get(gca,'ylim');
hold on
[f_r,xi_r]=ksdensity(r);
plot(xi_r,f_r,'m','LineWidth',2);
xline(r_ci(1,1),'r','LineWidth',2);
xline(r_ci(1,2),'r','LineWidth',2);
xline(r_mean,'b','LineWidth',2);
xline(r_median,'g','LineWidth',2);
set(gca,'box','off')
set(gca,'FontSize',20)
set(gca,'linewidth',2)
hold off

%histogram of De
subplot(2,2,4);
histogram(De,'Normalization','pdf');
xlabel('Denudation length (cm)');
ylabel('Probability');
%De_y=get(gca,'ylim');
hold on
[f_De,xi_De]=ksdensity(De);
plot(xi_De,f_De,'m','LineWidth',2);
xline(De_ci(1,1),'r','LineWidth',2);
xline(De_ci(1,2),'r','LineWidth',2);
xline(De_mean,'b','LineWidth',2);
xline(De_median,'g','LineWidth',2);
set(gca,'box','off')
set(gca,'FontSize',20)
set(gca,'linewidth',2)
hold off



function t=Be10NewtonN(Te,x,D,n,attenuation,density,decay)
% Newton methods to find t using eroded thickness (D). 
% Te: effective age
% x: first guess of the age
% D: eroded thickness
% n: number of iterations
% attenuation: attenuation length of nucleon
% density: sample sediment density
% decay: decay constant
B=-density*D/attenuation;
f=@(x)exp(B)*exp(-decay*x)-B*Te/x+(decay*Te-1);
df=@(x)-decay*exp(B-decay*x)+B*Te/(x^2);
for i=1:n
    x=x-f(x)/df(x);
end
t=x;
end

