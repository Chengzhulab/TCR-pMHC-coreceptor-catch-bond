% %  Matlab package for analysis for fitting to bond lifetime vs force data %%
% %  -------------------------------------------------------------------------
% %  Georgia Tech license was used to use MATLAB
% %  ==============================================
% %  ----------------------------------------------
% %  Introduction of basic tutorial code
% %  ----------------------------------------------
% %  1. This code is based on our previous paper (Nature Communications : doi.org/10.1038/s41467-023-38267-1).
% %  2. For simplifying, this demo is implemented through scanning only discrete mean values of bond lifetime vs force data set
% %  3. One can calculate error by considering sem of bond lifetime of bond with Jacobian and residual matrix
% %  4. This script is encoded in "chronological order for analyzing"
% %  5. All sub-codes were written using Matlab-bulit in function or home-built codes
% %  6. One can select one of three Models
% %
% %  Any modification of code is available upon user's preference
% %  
                    %%%%      Enjoy your analysis     %%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %------------------------------------%
                    %                                    %
                    %              |/ / / /              %
                    %          //  -    -   //           %
                    %             (  @ @  )              %
                    % +---------oOOo- & -oOOo--------+   %
                    % |          Good Luck           |   %
                    % |         Dr. H-K.Choi         |   %
                    % +-------------------Oooo-------+   %
                    %           oooO      (    )         %
                    %          (    )      )  /          %
                    %            )  (     (__/           %
                    %            (__/                    %
                    %------------------------------------%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear all previous workspace metadata
if 1
   clear all
   close all
   clc
end
%% Select Model
Modelselect = 2; % Biomolecular interaction = 1 // Trimolecular interaction = 2 // Two-pathway model = 3.

%% Loading raw data

rawdata = dir('E8 TCR+CD4+MHC_Tot_mean.txt');
S = struct(rawdata);
numdata= size(S,1);
dataraw={};

for i=1:numdata
    dataraw{i} = load(S(i).name)';
end

for i=1:numdata
    BlT_Force_set{i} = dataraw{i}';
end

n_cycle = 1;
%% Data setting
disp(['Start analysis of ', S(n_cycle).name,'!!'])
force_data_platform = BlT_Force_set{n_cycle}(:,1);
forcetranslation = 1;
% To avoid divergence/signularity  of fitting, negative or zero force was set by 1fN
for fij = 1:numel(force_data_platform )
    if force_data_platform(fij) < 0 || force_data_platform(fij) == 0
        force_data_platform_re(fij) = 0.001; % pN
    else
        force_data_platform_re(fij) = force_data_platform(fij); % pN
    end
end
bond_lifetime_platform(:,1) = BlT_Force_set{n_cycle}(:,3); % s in dimension
bond_lifetime_platform(:,2) = BlT_Force_set{n_cycle}(:,4); % s in dimension
%% Compartmentalization in force range
figure(2); hold on; 
sc = 25;
figure(2);
scatter(force_data_platform,bond_lifetime_platform(:,1),sc,'ko');
errorbar(force_data_platform,bond_lifetime_platform(:,1),bond_lifetime_platform(:,2),'r');
e.Marker = 'o';
e.Color = 'red';
e.CapSize = 15;
xlabel('Force (pN)')
ylabel('Bond lifetime (s)')
%% Basic info of biophysical properties of unstructured peptide
Temp = 273.15 + 30; % Absolute temperature, K
kbT = 4.114*Temp./298.15 ; % pN*nm
beta = 1/kbT;
lp = 0.36; % contour length of one amino-acid (nm)
Pp = 0.39; % persistance length of one amino-acid (nm)
g_r = 100; % elastic modulus of 3D folded interaction (pN)
Kaa = 50*10^(6); % elastic modulus of peptide (pN)
bk = 0.55; % Kuhn length of amino-acids (nm)
%% Construct-dependent constants related TCR-pMHC

d_N = 12.3 ;
pi_alpha3 = 23.6 ;
deld1 = 3.35;
d_Lz = 0.001;
pi_Cab = 0; % Degree (not used)
d_MHCa2b2 = 3.79; % nm (not used)
%% Set boundaries of non-linear fitting
if Modelselect == 1
    % order of parameterizations : [d0 theta0 k0 delx0]
    lb = [d_N/2 0.001 0.01 0.01]; % Lower bound can be set by considering actual PDB structure and physical reasonability
    ub = [d_N 89.99 200 2.69]; % upper bound can be set by considering actual PDB structure and physical reasonability
elseif Modelselect  == 2
    % order of parameterizations : [d0 theta0 k0 delx0 k0_s x0_s]
    lb = [d_N/2 0.001 0.01 0.01 0.001 0.001]; % Lower bound can be set by considering actual PDB structure and physical reasonability
    ub = [d_N 89.99 200 2.69 250 d_N]; % upper bound can be set by considering actual PDB structure and physical reasonability
elseif Modelselect == 3
    % order of parameterizations : [k0_c x0_c k0_s x0_s]
    lb = [0.001 -d_N 0.001 0.001]; % Lower bound can be set by considering actual PDB structure and physical reasonability
    ub = [250 -0.001 250 d_N]; % upper bound can be set by considering actual PDB structure and physical reasonability
end
%% Inital guess for parameters:
if Modelselect == 1
    nmax = 7; % Possible maximum unfolded amino-acids
    % order of parameterizations : [d0 theta0 k0 delx0]
    catchbond0 = [8 4 50 1] % Starting values of fitting parameters.
    % General intial value-set is [8.3 4 5 1.5] for class 2 model
    
    cselect = [linspace(0,1,nmax+1)' linspace(1,0,nmax+1)' linspace(1,0,nmax+1)']; % Color-setting for drawing lines
elseif Modelselect  == 2
    nmax = 7; % Possible maximum unfolded amino-acids
    % order of parameterizations : [d0 theta0 k0 delx0 k0_s x0_s]
    catchbond0 = [8 4 50 1 0.5 1] % Starting values of fitting parameters.
    % General intial value-set is [8.3 4 5 1.5] for class 2 model
    
    cselect = [linspace(0,1,nmax+1)' linspace(1,0,nmax+1)' linspace(1,0,nmax+1)']; % Color-setting for drawing lines
elseif Modelselect == 3 
    % order of parameterizations : [k0_c x0_c k0_s x0_s]
    catchbond0 = [5 -2 0.5 0.5]; % Starting values of fitting parameters.
    nmax = 0;
    
    cselect = [linspace(0,1,4)' linspace(1,0,4)' linspace(1,0,4)']; % Color-setting for drawing lines
end

%% Setting fitting

% Simple version using only mean bond lifetime vs force scatters.

options = optimoptions(@lsqcurvefit,'Display','iter-detailed','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',2000);

force_data =  force_data_platform_re';
lifetime_data = bond_lifetime_platform(:,1);

% with extra paramters built in
figure(3);
errorbar(force_data_platform,bond_lifetime_platform(:,1),bond_lifetime_platform(:,2),'k','LineStyle','none');
scatter(force_data_platform,bond_lifetime_platform(:,1),2*sc,'ko','filled');
errorbar(force_data_platform,bond_lifetime_platform(:,1),bond_lifetime_platform(:,2),'k');
hold on;

for in = 1:nmax+1
    n_tran = in-1;
    tvec = @(catchbond,force) modelvec(catchbond,force,n_tran,lp,Pp,beta,Temp,Kaa,g_r,d_N,pi_alpha3,pi_Cab,deld1,d_MHCa2b2,d_Lz,Modelselect);
    [catchbond,resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(tvec,catchbond0,force_data,lifetime_data,lb,ub,options);
    
    if Modelselect == 3
        k0_c = catchbond(1) % 1/s
        x0_c = catchbond(2) % nm
        k0_s = catchbond(3) % 1/s
        x0_s = catchbond(4) % nm
    elseif Modelselect == 1
        deld0(in) = catchbond(1) % nm
        theta0(in) = catchbond(2) % Degree
        k0(in) = catchbond(3) % 1/s
        delx0(in) = catchbond(4) % nm
        delG0(in) = -log(k0(in)/10^6); % -log(k0(in)/10^6) kBT
    elseif Modelselect == 2
        
        deld0(in) = catchbond(1) % nm
        theta0(in) = catchbond(2) % Degree
        k0_c(in) = catchbond(3) % 1/s
        delx0(in) = catchbond(4) % nm
        k0_s(in) = catchbond(5)
        x0_s(in) = catchbond(6)
    end
    % Check fit:
    tplot = linspace(0.01,42.5,426);
    t_fit_res{in} = tvec(catchbond,tplot); 
    figure(3);hold on;
    plot(tplot,t_fit_res{in},'color',cselect(in,:),'LineWidth',2);
    ylim([0 max(lifetime_data)*sqrt(2)]);
    xlabel('Force (pN)')
    ylabel('Bond lifetime (s)')
end

close all;
for in = 1:nmax+1
    figure(1);hold on;
    plot(tplot,t_fit_res{in},'color',cselect(in,:),'LineWidth',2);
    ylim([0 max(lifetime_data)*sqrt(2)]);
    if Modelselect ~= 3
        Legend{in}=strcat('Fitting curve when n* = ', num2str(in-1));
    end
end
figure(1);
scatter(force_data_platform,bond_lifetime_platform(:,1),2*sc,'ko','filled');
xlabel('Force (pN)')
ylabel('Bond lifetime (s)')
if Modelselect ~= 3
    legend(Legend)
else
    legend(['Fitted curve']);
end

%% Drawing results
if Modelselect ~= 3
    for in = 1:nmax+1

        if Modelselect == 1
            phi0(in) = atan(tan(theta0(in)*pi/180)*(d_N-deld1+d_Lz))./(deld1)*180/pi;
            b_n(in) = deld0(in).*(sin(theta0(in)*pi/180)*cot(phi0(in)*pi/180)+cos(theta0(in)*pi/180)) + d_Lz.*((sin(theta0(in)*pi/180)*cot(phi0(in)*pi/180)+cos(theta0(in)*pi/180))-1);
            L_n(in) = ((in-1))*lp*(sin(theta0(in)*pi/180)*cot(phi0(in)*pi/180)+cos(theta0(in)*pi/180));
            delta_expected(in) = ((in-1))*lp;
            n_tot(in) = (abs(((deld0(in)+d_Lz)*sin(theta0(in)*pi/180)./sin(phi0(in)*pi/180)-d_MHCa2b2)./lp+(1+sin(theta0(in)*pi/180)./sin(phi0(in)*pi/180))*in));
            
        elseif Modelselect == 2
            phi0(in) = atan(tan(theta0(in)*pi/180)*(d_N-deld1+d_Lz))./(deld1)*180/pi;
            b_n(in) = deld0(in).*(sin(theta0(in)*pi/180)*cot(phi0(in)*pi/180)+cos(theta0(in)*pi/180)) + d_Lz.*((sin(theta0(in)*pi/180)*cot(phi0(in)*pi/180)+cos(theta0(in)*pi/180))-1);
            L_n(in) = ((in-1))*lp*(sin(theta0(in)*pi/180)*cot(phi0(in)*pi/180)+cos(theta0(in)*pi/180));
            delta_expected(in) = ((in-1))*lp;
            n_tot(in) = (abs(((deld0(in)+d_Lz)*sin(theta0(in)*pi/180)./sin(phi0(in)*pi/180)-d_MHCa2b2)./lp+(1+sin(theta0(in)*pi/180)./sin(phi0(in)*pi/180))*in));
        end
    end
    %% Finding the best-fitting n*
    % Intersections indicates appropriately possible n*
    
    if Modelselect == 1
        figure(4); box on; grid on; hold on;
        plot(0:nmax,L_n,'--g')
        plot(0:nmax,delta_expected,'--r')
        plot(0:nmax,delx0,'--b')
        scatter(0:nmax,delx0,1.5*sc,'kd','filled')
        plot(0:nmax,delx0,'k','LineWidth',2)
        legend('Molecular extension with minimum n*','Molecular extension with maximum n*','Fitted molecular extension at force-free state')
        ylim([-2 12])
        ylabel('Fitted molecular extension (nm)')
        xlabel('Number of unfolded amino acid (n*,#)')
    elseif Modelselect == 2
        figure(4); box on; grid on; hold on;
        plot(0:nmax,L_n,'--g')
        plot(0:nmax,delta_expected,'--r')
        scatter(0:nmax,delx0+x0_s,1.5*sc,'kd','filled')
        plot(0:nmax,delx0+x0_s,'k','LineWidth',2)
        ylim([-2 6])
        legend('Molecular extension with minimum n*','Molecular extension with maximum n*','Fitted molecular extension at force-free state')
        ylabel('Fitted molecular extension (nm)')
        xlabel('Number of unfolded amino acid (n_{MHC}*,#)')
    end
end

%% Select n*
pause;
ans = input('Set number of amino-acid of total unfolding parts :  ','s');
n_tran2 = str2num(ans)+1
if Modelselect == 1
    Fit_result_ = [deld0(n_tran2) 0 theta0(n_tran2) 0 delx0(n_tran2) 0 n_tran2-1 k0(n_tran2) 0]
elseif Modelselect == 2
    Fit_result_ = [deld0(n_tran2) 0 theta0(n_tran2) 0 delx0(n_tran2) 0 n_tran2-1 k0_c(n_tran2) 0 x0_s(n_tran2) 0 k0_s(n_tran2) 0]
end

%% Fitting function %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Low-Level Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tvec = modelvec(catchbond,force,n_tran,lp,Pp,beta,Temp,Kaa,g_r,d_N,pi_alpha3,pi_Cab,deld1,d_MHCa2b2,d_Lz,Modelselect)
    % Vector of y model for a vector of force values:
    tvec = zeros(size(force));
    if Modelselect == 1
        for i = 1:length(force)
            tvec(i) = model1(catchbond,force(i),n_tran,lp,Pp,beta,Temp,Kaa,g_r,d_N,pi_alpha3,pi_Cab,deld1,d_MHCa2b2,d_Lz);
        end
    elseif Modelselect == 2
        for i = 1:length(force)
            tvec(i) = model2(catchbond,force(i),n_tran,lp,Pp,beta,Temp,Kaa,g_r,d_N,pi_alpha3,pi_Cab,deld1,d_MHCa2b2,d_Lz);
        end
    elseif Modelselect == 3
        for i = 1:length(force)
            tvec(i) = model3(catchbond,force(i),n_tran,lp,Pp,beta,Temp,Kaa,g_r,d_N,pi_alpha3,pi_Cab,deld1,d_MHCa2b2,d_Lz);
        end
    end
end

function t = model1(catchbond,F,n_tran,lp,Pp,beta,Temp,Kaa,g_r,d_N,pi_alpha3,pi_Cab,deld1,d_MHCa2b2,d_Lz)
     % Model value, given values
    % Fixed parameters, plus function R1, are transferred in the call
    
    function d_F = rbstret(d_i,F,beta,g_r)
        d_F = d_i.*(coth(beta*F.*(d_i))-1./(beta*F.*(d_i))).*(1+F./g_r);
    end

    function z = eWLC_inv(F,Lo,Lp,T,Ko,corr)
        if Lo == 0
            for j = 1:numel(F)
                z(j) = 0;
            end
        else
            z = zeros(size(F));
            l0 = arrayfun(@(j) WLC_inv(F(j),Lo,Lp,T,false),1:numel(F))/Lo;
            if corr == 0
                % without correction
                for j = 1:numel(F)
                    z(j) = fzero(@(l)(1.38e-2*T/Lp)*(1/4./(1-(l-F(j)/Ko)).^2-1/4+(l-F(j)/Ko))-F(j),l0(j))*Lo;
                end
            elseif corr == 1
                % with correction
                a = [-0.5164228, -2.737418, 16.07497, -38.87607, 39.49944, -14.17718];
                for j = 1:numel(F)
                    z(j) = fzero(@(l)(1.38e-2*T/Lp)*(1/4./(1-(l-F(j)/Ko)).^2-1/4+(l-F(j)/Ko)+sum(a.*l.^(2:7)))-F(j),l0(j))*Lo;
                end
            elseif corr == 2
                % with correction by Ogden et al.
                for j = 1:numel(F)
                    z(j) = fzero(@(l)(1.38e-2*T/Lp)*(1/4./(1-(l-F(j)/Ko)).^2-1/4+(l-F(j)/Ko)-3/4*(l-F(j)/Ko).^2)-F(j),l0(j))*Lo;
                end
                
            end
        end
    end

    if n_tran == 0
        catchbond(1) = 0;
        deld0 = catchbond(1); % nm
        catchbond(2) = 0;
        theta0 = catchbond(2); % Degree
        k0 = catchbond(3); % 1/s
        delx0 = catchbond(4); % nm
        delG0 = -log(k0./10^6); % -log(k0./10^6);% kBT
        t = 1/(k0*sqrt(1-2*delx0*F/(3*delG0/beta))*exp(delG0*(1-(1-2*delx0*F/(3*delG0/beta)).^(3/2))));
    else   
        deld0 = catchbond(1); % nm
        theta0 = catchbond(2); % Degree
        phi0 = atan(tan(theta0*pi/180)*(d_N-deld1+d_Lz))./(deld1)*180/pi;
        k0 = catchbond(3); % 1/s
        delx0 = catchbond(4); % nm
        delG0 = -log(k0./10^6); % -log(k0./10^6);% kBT
        integrand = @(f) eWLC_inv(f,lp*n_tran,Pp,Temp,Kaa,1).*(sin(theta0*pi/180)*cot(phi0*pi/180)+cos(theta0*pi/180))./(lp*n_tran)+rbstret(deld0,f,beta,g_r).*(sin(theta0*pi/180)*cot(phi0*pi/180)+cos(theta0*pi/180))./(lp*n_tran) + rbstret(d_Lz,f,beta,g_r).*((sin(theta0*pi/180)*cot(phi0*pi/180)+cos(theta0*pi/180))-1)./(lp*n_tran) - rbstret(d_N,f,beta,g_r)./(lp*n_tran);
        int_unitx = integral(integrand,0,F);
        t = 1/(k0*sqrt(1-2*delx0*int_unitx/(3*delG0/beta))*exp(delG0*(1-(1-2*delx0*int_unitx/(3*delG0/beta)).^(3/2))));
    end
end

function t = model2(catchbond,F,n_tran,lp,Pp,beta,Temp,Kaa,g_r,d_N,pi_alpha3,pi_Cab,deld1,d_MHCa2b2,d_Lz)
    % Model value, given values
    % Fixed parameters, plus function R1, are transferred in the call
    
    function d_F = rbstret(d_i,F,beta,g_r)
        d_F = d_i.*(coth(beta*F.*(d_i))-1./(beta*F.*(d_i))).*(1+F./g_r);
    end

    function z = eWLC_inv(F,Lo,Lp,T,Ko,corr)
        if Lo == 0
            for j = 1:numel(F)
                z(j) = 0;
            end
        else
            z = zeros(size(F));
            l0 = arrayfun(@(j) WLC_inv(F(j),Lo,Lp,T,false),1:numel(F))/Lo;
            if corr == 0
                % without correction
                for j = 1:numel(F)
                    z(j) = fzero(@(l)(1.38e-2*T/Lp)*(1/4./(1-(l-F(j)/Ko)).^2-1/4+(l-F(j)/Ko))-F(j),l0(j))*Lo;
                end
            elseif corr == 1
                % with correction
                a = [-0.5164228, -2.737418, 16.07497, -38.87607, 39.49944, -14.17718];
                for j = 1:numel(F)
                    z(j) = fzero(@(l)(1.38e-2*T/Lp)*(1/4./(1-(l-F(j)/Ko)).^2-1/4+(l-F(j)/Ko)+sum(a.*l.^(2:7)))-F(j),l0(j))*Lo;
                end
            elseif corr == 2
                % with correction by Ogden et al.
                for j = 1:numel(F)
                    z(j) = fzero(@(l)(1.38e-2*T/Lp)*(1/4./(1-(l-F(j)/Ko)).^2-1/4+(l-F(j)/Ko)-3/4*(l-F(j)/Ko).^2)-F(j),l0(j))*Lo;
                end
                
            end
        end
    end

    if n_tran == 0
        catchbond(1) = 0;
        deld0 = catchbond(1); % nm
        catchbond(2) = 0;
        theta0 = catchbond(2); % Degree
        kc = catchbond(3); % 1/s
        delx0 = catchbond(4); % nm
        ks = catchbond(5); % 1/s
        xs = catchbond(6); % nm
        
        delG0 = -log(kc./10^6); % -log(k0./10^6);% kBT
        
        t = 1/(kc*sqrt(1-2*delx0*F/(3*delG0/beta))*exp(delG0*(1-(1-2*delx0*F/(3*delG0/beta)).^(3/2)))+ks.*exp(-beta*xs.*F));
        
    else   
        deld0 = catchbond(1); % nm
        theta0 = catchbond(2); % Degree
        phi0 = atan(tan(theta0*pi/180)*(d_N-deld1+d_Lz))./(deld1)*180/pi;
        kc = catchbond(3); % 1/s
        delx0 = catchbond(4); % nm
        delG0 = -log(kc./10^6); % -log(k0./10^6);% kBT
        integrand = @(f) eWLC_inv(f,lp*n_tran,Pp,Temp,Kaa,1).*(sin(theta0*pi/180)*cot(phi0*pi/180)+cos(theta0*pi/180))./(lp*n_tran)+rbstret(deld0,f,beta,g_r).*(sin(theta0*pi/180)*cot(phi0*pi/180)+cos(theta0*pi/180))./(lp*n_tran) + rbstret(d_Lz,f,beta,g_r).*((sin(theta0*pi/180)*cot(phi0*pi/180)+cos(theta0*pi/180))-1)./(lp*n_tran) - rbstret(d_N,f,beta,g_r)./(lp*n_tran);
        int_unitx = integral(integrand,0,F);
        t1 = (kc*sqrt(1-2*delx0*int_unitx/(3*delG0/beta))*exp(delG0*(1-(1-2*delx0*int_unitx/(3*delG0/beta)).^(3/2))));
        
        ks = catchbond(5); % 1/s
        xs = catchbond(6); % nm

        t = 1/(t1+ks.*exp(-beta*xs.*F));
    end
end

function t = model3(catchbond,F,n_tran,lp,Pp,beta,Temp,Kaa,g_r,d_N,pi_alpha3,pi_Cab,deld1,d_MHCa2b2,d_Lz)
    % Model value, given values
    % Fixed parameters, plus function R1, are transferred in the call

    % Model value, given values
    kc = catchbond(1); % nm
    xc = catchbond(2); % Degree
    ks = catchbond(3); % 1/s
    xs = catchbond(4); % nm

    t = 1/(kc.*exp(beta*xc.*F)+ks.*exp(beta*xs.*F));
end