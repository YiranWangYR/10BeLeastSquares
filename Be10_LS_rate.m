% version 1.0
% Code to apply a least squares method exposure age estimation from
% cosmogenic nuclide depth-profiles. Detailed explaination can be find in:
% Wang and Oskin (202x)
% This code estimate age with erosion/denudation rate r, muogenic production ignored.
% A Monte Carlo simulation is applied to find the distribution of age and
% inheritance. Erosion rate can be set to zero.

% If you encounter any problem, please contact: yrwwang@ucdavis.edu

% How to use this code:
% Please edit the values betwenn line 25-57, following the instructions
% below

clear all;

%=============================Inputs==============================
% -----Main variables-------
% !! Format !!
% variable=[x,y,z]
% x=the mean of the variable;
% y=the standard diviation of the variable (normal distribution), or the
% error (normal distribution),y may equal to 0;
% z=type of distribution (1=normal distribution; 0=uniform distribution).
P0nR=[20,0,1]; %surface production rate for nucleon spallation, atoms/(g*yr)
P0m1R=[0.31,0,1]; %surface production rate for negative muons, atoms/(g*yr)
P0m2R=[0.13,0,1]; %surface production rate for fast muons, atoms/(g*yr)
densityR=[2,0.2,1]; %g/cm3, sediment density
LanR=[160,0,1]; % nucleon spallation attenuation length; g/cm2
Lam1R=[1500,0,1]; %negative muon attenuation length; g/cm2
Lam2R=[4300,0,1]; %fast muon attenuation length
rateR=[0.0003, 0.00005, 1]; % cm, erosion rate. No erosion when the first and second terms ==0
% When there is no erosion, this code returns results with muogenic
% production fully incorperated. (Eq.3-5 of Wang and Oskin (202x)
%--------------

% load sample data. Change the file name accordingly
D1=load('../Be10 code/data/T2_Fin.txt');
% Data format:
% 1st column: sample depth, cm
% 2nd column: standard deviation of sample depth, cm
% 3rd column: sample concentration, atoms/g
% 4th column: standard deviation of C, atoms/g (!!!put in actual value here, 
% do not use the percentage form)

% distribution type of samples
z_dis=1;    %distribution type of sample depth. 1=normal distribution; 0=uniform distribution
y_dis=1;    %distribution type of sample concentration. 1=normal distribution; 0=uniform distribution

% decay constant
decay=0.0000004997;    %decay constant

P=5000; % Number of iterations for the Monte Carlo simulation


% output excel file name
filename=(['10Be_rate.xlsx']);

%======================Input ends=====================

%==============Data preperation===============
z_mean=D1(:,1); % cm, sample depth (from top)
z_sd=D1(:,2); % cm, standard deviation of z
y_mean=D1(:,3); % sample concentration, atoms/g
y_sd=D1(:,4); % standard deviation of C1

N=length(y_mean);  % number of samples


%============Stope if sd of C is in percentage format===========
if y_sd(1)<1
    fprintf('Please use actual value of the standard deviation of concentration. Do not use the percentage format.');
    return;
end

%-------------Sampling data (Monte Carlo)--------------
%P0n
if P0nR(3)==1
    P0n=P0nR(1)+P0nR(2)*randn(P,1);
else P0n=P0nR(1)+2*P0nR(2)*(rand(P,1)-0.5);
end

%P0m1
if P0m1R(3)==1
    P0m1=P0m1R(1)+P0m1R(2)*randn(P,1);
else P0m1=P0m1R(1)+2*P0m1R(2)*(rand(P,1)-0.5);
end

%P0m2
if P0m2R(3)==1
    P0m2=P0m2R(1)+P0m2R(2)*randn(P,1);
else P0m2=P0m2R(1)+2*P0m2R(2)*(rand(P,1)-0.5);
end

%density
if densityR(3)==1
    density=densityR(1)+densityR(2)*randn(P,1);
else density=densityR(1)+2*densityR(2)*(rand(P,1)-0.5);
end

%attenuation
if LanR(3)==1
    Lan=LanR(1)+LanR(2)*randn(P,1);
else Lan=LanR(1)+2*LanR(2)*(rand(P,1)-0.5);
end

if Lam1R(3)==1
    Lam1=Lam1R(1)+Lam1R(2)*randn(P,1);
else Lam1=Lam1R(1)+2*Lam1R(2)*(rand(P,1)-0.5);
end

if Lam2R(3)==1
    Lam2=Lam2R(1)+Lam2R(2)*randn(P,1);
else Lam2=Lam2R(1)+2*Lam2R(2)*(rand(P,1)-0.5);
end


% erosion rate
if rateR(3)==1
    r=rateR(1)+rateR(2)*randn(P,1);
else r=rateR(1)+2*rateR(2)*(rand(P,1)-0.5);
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


%---------Least Squares---------
M=2;  % number of unknowns
G=zeros(N,M);

Re=zeros(P,6);  % results to be saved here


for i=1:P
    La=[Lan(i),Lam1(i),Lam2(i)];
    % calculate production rate at each sample depth
    Pzn=P0n(i)*exp(-density(i)*z(:,i)/La(1));
    Pzm1=P0m1(i)*exp(-density(i)*z(:,i)/La(2));
    Pzm2=P0m2(i)*exp(-density(i)*z(:,i)/La(3));
    if r(i)==0   % no erosion case
        x=Pzn+Pzm1+Pzm2;    % production rate at sample depth
    else x=Pzn;     % erosion exists, approximate by ignoring muons
    end
    % build matrix for the linear equation. x is production rate at depth,
    % y is sample concentration; eq. 4 and 6
    G(:,1)=1;
    G(:,2)=x;
    [Up, Lp, Vp] = svd( G);   % singular value decomposition
    mest = lsqnonneg(G, y(:,i));
    Te=mest(2,1);        % effective exposure age (age without decay or erosion)
    C_inh=mest(1,1);        % inherited concentration
    B=density(i)*r(i)/La(1)+decay;    %eq. 7
    t=-log(1-Te*B)/B;
    Cinh_o = C_inh/exp(-decay*t); %correct inheritance for decay
    De = r(i)*t;
    %       Te, inheretance, t1, origional inheritance, erosion rate, total erosion
    Re(i,:)=[Te, C_inh, t, Cinh_o, r(i), De];

end
%==========================Save results====================
xlswrite(filename,Re); %results saved as excel, please change file name.
%first column: Te value (effective exposure age,yr)
%second column: Cinh (inhereted concentration in samples)
%third column: t(exposure age),
%forth column: Cinh_o(original inhereted concentration, corrected for decay),
%fifth column: r(erosion rate)
%sixth column: total erosion



%==========================Display results=================
Te=Re(:,1)/1000; %effective exposure age in kyr
C_inh=Re(:,2);
t=Re(:,3)/1000; %exposure age in kyr
Cinh_o=Re(:,4);
r=Re(:,5)*1000;  %erosion rate in cm/kyr
De=Re(:,6); %eroded thickness

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



%==============plot data===========

figure(1)
% histogram of t1
histogram(t);
xlabel('Exposure age (kyr)','Fontsize',20);
ylabel('Frequency','Fontsize',20);
t1_y=get(gca,'ylim');
hold on
plot([t_ci(1,1) t_ci(1,1)],t1_y,'r','LineWidth',2);
plot([t_ci(1,2) t_ci(1,2)],t1_y,'r','LineWidth',2);
plot([t_mean t_mean],t1_y,'b','LineWidth',2);
plot([t_median t_median],t1_y,'g','LineWidth',2);
set(gca,'FontSize',15)
set(gca,'linewidth',2)
hold off


figure(2)
%histogram of Te
subplot(1,2,1);
histogram(Te);
xlabel('T_e value (kyr)','Fontsize',20);
ylabel('Frequency','Fontsize',20);
Te_y=get(gca,'ylim');
hold on
plot([Te_ci(1,1) Te_ci(1,1)],Te_y,'r','LineWidth',2);
plot([Te_ci(1,2) Te_ci(1,2)],Te_y,'r','LineWidth',2);
plot([Te_mean Te_mean],Te_y,'b','LineWidth',2);
plot([Te_median Te_median],Te_y,'g','LineWidth',2);
set(gca,'FontSize',15)
set(gca,'linewidth',2)
hold off

%histogram of C_inh
subplot(1,2,2);
histogram(C_inh);
hold on
xlabel({'Inherited concentration (atoms/g)'},'Fontsize',20);
ylabel('Frequency','Fontsize',20);
Cinh_y=get(gca,'ylim');
plot([Cinh_ci(1,1) Cinh_ci(1,1)],Cinh_y,'r','LineWidth',2);
plot([Cinh_ci(1,2) Cinh_ci(1,2)],Cinh_y,'r','LineWidth',2);
plot([Cinh_mean Cinh_mean],Cinh_y,'b','LineWidth',2);
plot([Cinh_median Cinh_median],Cinh_y,'g','LineWidth',2);
set(gca,'FontSize',15)
set(gca,'linewidth',2)
hold off

figure(3)
%histogram of r
subplot(1,2,1);
histogram(r);
xlabel('Erosion rate (cm/kyr)','Fontsize',20);
ylabel('Frequency','Fontsize',20);
r_y=get(gca,'ylim');
hold on
plot([r_ci(1,1) r_ci(1,1)],r_y,'r','LineWidth',2);
plot([r_ci(1,2) r_ci(1,2)],r_y,'r','LineWidth',2);
plot([r_mean r_mean],r_y,'b','LineWidth',2);
plot([r_median r_median],r_y,'g','LineWidth',2);
set(gca,'FontSize',15)
set(gca,'linewidth',2)
hold off

%histogram of De
subplot(1,2,2);
histogram(De);
xlabel('Eroded thickness (cm)','Fontsize',20);
ylabel('Frequency','Fontsize',20);
De_y=get(gca,'ylim');
hold on
plot([De_ci(1,1) De_ci(1,1)],De_y,'r','LineWidth',2);
plot([De_ci(1,2) De_ci(1,2)],De_y,'r','LineWidth',2);
plot([De_mean De_mean],De_y,'b','LineWidth',2);
plot([De_median De_median],De_y,'g','LineWidth',2);
set(gca,'FontSize',15)
set(gca,'linewidth',2)
hold off



%---------------Figure 4 inverted results vs data----------------
K=500;   %number of simulation

% best fit lines
x1=zeros(P,1);
x2=P0nR(1)*ones(P,1);
y1=Re(:,2);
TeN=Re(:,1);
y2=TeN.*x2+y1;

% error of production rate at depth
Pz_err=sqrt(P0nR(2)^2+(densityR(1)*z_sd/LanR(1)).^2+(densityR(2)*z_mean/LanR(1)).^2+(LanR(2)*densityR(1)*z_mean/LanR(1)^2).^2).*(P0nR(1)*exp(-densityR(1)*z_mean/LanR(1)));


figure(4)
subplot(1,2,1)
axis([0 P0nR(1) 0 max(y2)]);
xlabel('Production rate at depth (Pz; atome*g^{-1}*yr^{-1})','Fontsize',20)
ylabel('^{10}Be concentration (atoms/g)','Fontsize',20)
hold on
%============line fit================
for j=1:K
    i=round((P-1)*rand(1))+1;
   plot([x1(i),x2(i)],[y1(i),y2(i)],'color',[0.5 0.5 0.5])
end

%data
errorbar(x,y_mean,y_sd,'bo','linewidth',1.5); %vertical errorbar
y_lim=get(gca,'ylim');
y_l=y_lim(1,2)/100;
%this plot the horizontal errorbar. Newer versions of matlab can simplify
%this step
for i=1:6
    plot([x(i)-Pz_err(i),x(i)+Pz_err(i)],[y_mean(i),y_mean(i)],'b','linewidth',1.5);
    plot([x(i)-Pz_err(i),x(i)-Pz_err(i)],[y_mean(i)-y_l,y_mean(i)+y_l],'b','linewidth',1.5);
    plot([x(i)+Pz_err(i),x(i)+Pz_err(i)],[y_mean(i)-y_l,y_mean(i)+y_l],'b','linewidth',1.5);
end

set(gca,'FontSize',15)
set(gca,'linewidth',2)

hold off



% show individual depth profile fits
y_min=y_mean-y_sd;
y_max=y_mean+y_sd;

zmodel=(1:max(z_mean)+10)';
for j=1:K
    i=round((P-1)*rand(1))+1;
    Pzn=P0n(i)*exp(-density(i)*zmodel/Lan(i));
    Pzm1=P0m1(i)*exp(-density(i)*zmodel/Lam1(i));
    Pzm2=P0m2(i)*exp(-density(i)*zmodel/Lam2(i));
    B=density(i)*(r(i)/1000)./[Lan(i), Lam1(i), Lam2(i)]+decay;
    Te=(1-exp(-B*t(i)*1000))./B;
    fmodel(:,j)=Pzn*Te(1)+Pzm1*Te(2)+Pzm2*Te(3)+C_inh(i);
end

subplot(1,2,2)
axis([0 round(max(y_mean)*1.1, -4) -(max(z_mean)+10) 0])
xlabel('^{10}Be concentration (atoms/g)','Fontsize',20)
ylabel('Depth (m)','Fontsize',20)
hold on
for i=1:K
    plot(fmodel(:,i),-zmodel,'Color', [0.5 0.5 0.5]);
end
%this plot the horizontal errorbar. Newer versions of matlab can simplify
%this step
for i=1:N
    plot([y_min(i,1),y_max(i,1)],[-z_mean(i,1),-z_mean(i,1)],'-b','linewidth',1.5);
    plot([y_min(i,1),y_min(i,1)],[-z_mean(i,1)-2,-z_mean(i,1)+2],'-b','linewidth',1.5);
    plot([y_max(i,1),y_max(i,1)],[-z_mean(i,1)-2,-z_mean(i,1)+2],'-b','linewidth',1.5);
end
errorbar(y_mean,-z_mean,z_sd,'bo','linewidth',1.5); %vertical errorbar
set(gca,'FontSize',15)
set(gca,'linewidth',2)

hold off
