%% General Code Description of BPM_IndoorUFP_Main.m and Author Introduction
%======================================================================================================================
%> @ Code Description:
%> @ File        BPM_IndoorUFP_Main.m
%> @ Brief       Advanced module for parameter identification in indoor ultrafine particle dynamics using Bayesian inference.
%> @ Details     BPM_IndoorUFP_Main is the main program for analyzing indoor ultrafine particle (UFP) dynamics. 
%>               This module integrates particle number concentration data from particle counters with indoor particle dynamics models.
%>
%>               The parameter identification process using Bayesian inference includes:
%>               (1) Data Integration: Incorporates particle number concentration data measured by particle counters.
%>               (2) Model Setup for Indoor Particle Dynamics: Configures and aligns the indoor particle dynamics model with real-world data.
%>               (3) Bayesian Inference for Parameter Estimation: Employs Bayesian inference to identify and optimize model parameters, ensuring accurate representation of indoor ultrafine particle behavior.
%>               (4) Iterative Analysis and Visualization: Conducts iterative simulations and analyses, culminating in the visualization of model predictions against actual measured data.
%>
%> @ Author Introduction
%>              - Yesol Hyun (yesol2@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
%>              - Jung-Il Choi (jic@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
%>
%> @ Date        November 2023
%> @ Version     1.0
%> @ Copyright   Copyright (c) 2023 Yesol Hyun and Jung-Il Choi, Yonsei University, All rights reserved.
%> @ License     This project is released under the terms of the MIT License (see LICENSE file).
%======================================================================================================================
clc;clear all;close all;warning off; 
rng(100,'twister')
uqlab
%% Experimental Setup and Data Preparation for Indoor Ultrafine Particle Model
%======================================================================================================================
%> @details     This section of the PI_IndoorUFP_Main.m module deals with the preparation of experimental settings and initial data for Bayesian analysis in indoor ultrafine particle models. The process includes:

%>              (1) File Selection and Verification:
%>                  - Determines the appropriate settings file based on the specified emission source.
%>                  - Checks for the existence of this file and proceeds accordingly.

%>              (2) Initial Settings Execution:
%>                  - If the settings file does not exist, executes the PI_IndoorUFP_InitialSetting function to create it.
%>                  - This function sets up initial conditions and parameters necessary for forward model calculations.

%>              (3) Data Loading and Global Variable Declaration:
%>                  - Once the settings file is confirmed or created, loads the data from the file.
%>                  - Declares global variables for use in subsequent parts of the program.

%>              This approach ensures that all necessary experimental settings are correctly configured before proceeding with the Bayesian analysis and simulation of indoor ultrafine particle dynamics.
%======================================================================================================================
% Select and verify the appropriate settings file based on the emission source
Num_Source  = 2; % Emission source identifier
Num_hamaker = 6; % H identifier
if Num_Source ==1
    SourceName = 'Gas'
    filename = strcat(SourceName,'_setting.mat');
    control_volume = 340;
    Num_hamaker_ref =2;
    optpar_ref = [-0.84 2.49 2.29 4.65 1.55 0];
elseif Num_Source ==2
    SourceName = 'Candle'
    filename = strcat(SourceName,'_setting.mat');
    control_volume = 34;
    Num_hamaker_ref =6;
    optpar_ref = [-1.22 1.53 7.96 3.86 1.38 0];
end

% Check if the settings file exists and load or create it as needed
if exist(filename, 'file')
    disp(['File ', filename, ' exists.']);
    load(filename); % Load the existing file
else
    disp(['File ', filename, ' does not exist. Execute initial setting.']);
    PI_IndoorUFP_InitialSetting(Num_Source); % Execute initial settings for the given source
    load(filename); % Load the newly created settings file
end
% Declare global variables for use throughout the program
global ptcl_basic exp_data
%% Processing the Output Data for Bayesian Analysis and Setting Up the Model
%======================================================================================================================
%> @details     This part of the code handles data normalization and prepares for Bayesian analysis in the indoor ultrafine particle model. It includes:

%>              (1) Normalizing Experimental Data:
%>                  - Processes data from particle counters, normalizing particle number concentrations.
%>                  - Selects and normalizes data within a specific range around the off-source scan number.

%>              (2) Compiling and Structuring Data for Bayesian Analysis:
%>                  - Aggregates the normalized data into a structured format suitable for Bayesian analysis.
%>                  - Prepares global variables and sets up the forward model for the Bayesian framework.

%>              (3) Bayesian Analysis Setup:
%>                  - Defines prior distributions for various model parameters, considering expected ranges.
%>                  - Configures the Bayesian inversion solver, specifying the type, sampler, number of chains, and steps.
%>                  - Sets up the discrepancy model and options for the Bayesian analysis.
%======================================================================================================================
ratio = 0.55;
Source_off = ctot_10(nscan_offsource);
Nscan_Start= find(ctot_10(1:idx_max)<=ratio*SourceConcPeak_10,1,'last');
Nscan_End= idx_max+find(ctot_10(idx_max:Nscanexp)<=ratio*SourceConcPeak_10,1,'first')-1;
indices =Nscan_Start:Nscan_End;
idx_start =Nscan_Start;

filter_tot = indices(end);
disp([filter_tot,ctot_10(filter_tot)/ExpInfo(4)])

noise = 100;

valid_indices = [];
for i = indices
    c_new = cnum_exp(idx_10:NBIN, i);
    if any(c_new >= noise)  
        valid_indices = [valid_indices, i]; 
    end
end

for i = valid_indices
    c_new = cnum_exp(idx_10:NBIN, i);
    idx(i)=idx_10+find(c_new<=noise,1,'first')-1;
end

max_bin = max(idx);

for i = indices
    c_new= cnum_exp(idx_10:max_bin,i);
    [max_c(i),idx_maxc] = max(c_new);
    c_data = c_new/max_c(i);
    
    if i==idx_start
        c_tmp = c_data';
    else
        c_tmp = [c_tmp,c_data'];
    end
end

ctot_iuq = ctot_all(idx_max_all:length(ctot_all));
leng=length(indices);
c_tmp =[c_tmp,ctot_10(2:filter_tot)'/ExpInfo(4)];
myData.y = c_tmp;
myData.Name = 'Exp data';
beta = BETA(:,:,Num_hamaker);

% Setting up the forward model for Bayesian analysis
ModelOpts.mHandle = @(par) uq_UFP(par,ptcl_basic.radius,ptcl_basic.Dradius,ptcl_basic.volume,Vfrac,ExpInfo,cnumini_10,filter_tot,NBIN,beta ,control_volume,idx_10,indices,max_c,max_bin,ExpInfo(4));
ModelOpts.isVectorized = true;
myForwardModel1 = uq_createModel(ModelOpts);
BayesOpts1.ForwardModel = myForwardModel1;

% Defining prior distributions for Bayesian parameters
PriorOpts.Marginals(1).Name = 'a';
PriorOpts.Marginals(1).Type = 'Uniform';
PriorOpts.Marginals(1).Parameters = [-1.25 -0.7];

PriorOpts.Marginals(2).Name = 'b';
PriorOpts.Marginals(2).Type = 'Uniform';
PriorOpts.Marginals(2).Parameters = [1 4];

PriorOpts.Marginals(3).Name = 'f';
PriorOpts.Marginals(3).Type = 'Uniform';
PriorOpts.Marginals(3).Parameters = [0 1.9];

PriorOpts.Marginals(4).Name = 'GMD';
PriorOpts.Marginals(4).Type = 'Uniform';
PriorOpts.Marginals(4).Parameters = [2 10];

PriorOpts.Marginals(5).Name = 'Gsigma';
PriorOpts.Marginals(5).Type = 'Uniform';
PriorOpts.Marginals(5).Parameters = [1.3 2];

% Creating the prior distribution as an input object
myPriorDist = uq_createInput(PriorOpts);
BayesOpts.Prior = myPriorDist;

% Setting up the discrepancy model for the Bayesian analysis
SigmaOpts.Marginals(1).Name = 'sigma2';
SigmaOpts.Marginals(1).Type = 'Uniform';
SigmaOpts.Marginals(1).Parameters = [0 (0.05)^2];

mySigmaDist = uq_createInput(SigmaOpts);
DiscrepancyOptsUnknownDisc.Type = 'Gaussian';
DiscrepancyOptsUnknownDisc.Prior = mySigmaDist;
BayesOpts.Discrepancy = DiscrepancyOptsUnknownDisc;

% Configuring the Bayesian inversion solver and sampler
Solver.Type = 'MCMC';
Solver.MCMC.Sampler = 'AIES';
Solver.MCMC.NChains = 100;
Solver.MCMC.Steps = 1000;
Solver.MCMC.Visualize.Interval = 10;
Solver.MCMC.Visualize.Parameters = [1 2 3 4 5];

% Finalizing Bayesian inversion setup
BayesOpts.Type = 'Inversion';
BayesOpts.Data = myData;
BayesOpts.Solver = Solver;
%%  Execution of Bayesian Analysis and Forward Modeling Using IUQ Results
%======================================================================================================================
%> @details     This section of the code executes Bayesian analysis for parameter identification in indoor ultrafine particle models and performs forward modeling based on the inferred parameters. The process includes:

%>              (1) Bayesian Analysis Execution:
%>                  - Executes Bayesian analysis with pre-defined settings.
%>                  - Measures execution time for performance analysis.

%>              (2) Post-Processing of Bayesian Results:
%>                  - Filters out MCMC chains with low acceptance rates to improve the quality of analysis.
%>                  - Visualizes and evaluates the acceptance rates of MCMC chains.

%>              (3) Results Reporting and Visualization:
%>                  - Generates detailed reports of Bayesian analysis results.
%>                  - Visualizes and saves figures representing various aspects of the analysis.

%>              (4) Forward Modeling Based on Inferred Parameters:
%>                  - Applies the inferred parameters to run forward models and simulate indoor ultrafine particle dynamics.
%>                  - Compares the simulation results with experimental data for validation.
%======================================================================================================================
% Execution of Bayesian Analysis
tic
myBayesianAnalysis_Step = uq_createAnalysis(BayesOpts);
uq_time(1)=toc;

% Filtering chains with low acceptance rate
uq_display( myBayesianAnalysis_Step, 'acceptance', true)
acceptance =  myBayesianAnalysis_Step.Results.Acceptance;
[~,tolL,tolU,tolC] = isoutlier(acceptance, 'ThresholdFactor', 2);
TF=acceptance<min(max(tolL,0.1),tolC);
badchains = find(TF);
yline(tolL,'b--');
yline(tolU,'b--');
yline(tolC,'r--');
uq_postProcessInversion( myBayesianAnalysis_Step,'badChains',badchains);
hold on
scatter(badchains, acceptance(badchains), 'red', 'filled')
legend('All chains', 'lb', 'ub', 'median')

% Print out and display the results of Bayesian analysis
uq_print(myBayesianAnalysis_Step)
uq_display(myBayesianAnalysis_Step)

% Save figures from Bayesian analysis
Nfig = 9;
for j =1:Nfig
    fig = figure(j);
    filename_IUQResults = sprintf('IUQ Results of %s,figure%d',SourceName,j);
    saveas(fig,filename_IUQResults,'jpg');
end

optpar_tmp=myBayesianAnalysis_Step.Results.PostProc.PointEstimate.X{1,1};
optpar_tmp= optpar_tmp(1:5) ;

results = myBayesianAnalysis_Step1.Results.PostProc;
filename_data = strcat('cpc_result_',source,'_',date,'.mat');
save (filename_data,'results')
save(filename_data,'SourceConcPeak_all','-append')
save(filename_data,'SourceConcPeak_10','-append')
fclose('all');

a = optpar_tmp(1);
b = optpar_tmp(2);
f = 10.^optpar_tmp(3);
GMD = optpar_tmp(4);
Gsigma = optpar_tmp(5);
delay = 0;
optpar_map = [a,b,f,GMD,Gsigma];

% Save the optimal parameters to a CSV file
filename_csv = sprintf('IUQ_Results_OptimalParameter_%s.csv',SourceName);
header = {'a','b','f','GMD','Gsgima'};
par_save = [a b f GMD Gsigma];
fid1 = fopen(filename_csv, 'w');
fprintf(fid1, '%s,', header{1:end-1});
fprintf(fid1, '%s\n', header{end});
fclose(fid1);
dlmwrite(filename_csv, par_save, 'delimiter', ',', '-append');

% Forward modeling using inferred parameters
tic
[ctot_coag,Qcnum_coag,~] =  Aerosol_Model(control_volume,ptcl_basic.radius,ptcl_basic.Dradius,ptcl_basic.volume,Vfrac,ExpInfo,cnumini_10,NBIN,beta,a,b,f,GMD,Gsigma,delay);
toc
ctot_10_est = Cal_TotNum_new(Qcnum_coag,idx_10,NBIN,1,NumCurve); % True value

fig1 = figure;
plot(t,ctot_all*1e-6,'ro','LineWidth',1.1)
hold on
plot(t,ctot_10'*1e-6,'bo','LineWidth',1.1)
hold on
plot(t,ctot_coag*1e-6,'k-','LineWidth',1.1)
hold on
plot(t,ctot_10_est*1e-6,'k--','LineWidth',1.1)
xlabel('Time [min]','fontsize',13)
ylabel('Total number [x10^6/cm^3]','fontsize',13)
title('Total number concentration','fontsize',13)
h= legend('Exp','','IUQ','Ref','');
set(h,'fontsize',11)

filename = sprintf('Compare_total_number_concentration_of_%s.jpg', SourceName);
saveas(fig1, filename);

ioffsource= ExpInfo(3)+1;
fig2 = figure;
d= ptcl_basic.radius*2*1e+7;
plot(d,cnum_exp(:,ioffsource)*1e-6,'o','LineWidth',1.1)
hold on
plot(d,Qcnum_coag(:,ioffsource)*1e-6,'LineWidth',1.1)
xlabel('Dp[nm]','fontsize',13)
ylabel('Number concentration [x10^6/cm^3]','fontsize',13)
title('Size-resolved particle distribution at the end of source emission','fontsize',13)
h= legend('Exp','IUQ','Ref');
set(h,'fontsize',11)

filename = sprintf('Compare_the_feature_of_source_emission_of_%s.jpg', SourceName);
saveas(fig2, filename);
%% Module for IUQ Procedure in Ultrafine Particle Modeling
%======================================================================================================================
%> @details     This module, uq_UFP, performs Bayesian analysis and forward modeling for ultrafine particle dynamics in indoor environments. The process includes:

%>              (1) Parameter Transformation and Forward Modeling Execution:
%>                  - Transforms Bayesian parameters and runs the Aerosol_Model function for each set of parameters.
%>                  - Calculates total number concentrations and normalizes the estimated concentration data.

%>              (2) Data Stacking and Preparation for IUQ Analysis:
%>                  - Stacks normalized concentration data for Bayesian analysis.
%>                  - Prepares the data in a structured format suitable for subsequent IUQ steps.

%>              (3) IUQ Data Compilation:
%>                  - Compiles all processed data from the forward modeling into a format required for IUQ analysis.
%>                  - Ensures the data is correctly aligned and structured for effective Bayesian inference.
%======================================================================================================================
function t_IUQ= uq_UFP(X,radius,Dradius,volume,Vfrac,expdata,cnumini,nscan_filter,nbin,beta,control_volume,idx_10,indices,max_c,max_bin,SourceConcPeak)
Npar = size(X,1);
sourcepeak = expdata(4);
dt = expdata(1);
dtmin = dt/60;
nscan_final = nscan_filter;
expdata(2) = nscan_final;
% Iterating over each set of parameters
for j =1:Npar
    par_tmp= X(j,1:5);
    a=par_tmp(1);
    b=par_tmp(2);
    f=10.^par_tmp(3);
    GMD=par_tmp(4);
    Gsigma=par_tmp(5);
    delay = 0;
    
     % Run forward modeling using Aerosol_Model function   
    [~,cnum_est,~] = Aerosol_Model(control_volume,radius,Dradius,volume,Vfrac,expdata,cnumini,nbin,beta,a,b,f,GMD,Gsigma,delay); %disp(size(cnum_est))
    ctot_10 = Cal_TotNum_new(cnum_est,idx_10,nbin,1,nscan_final);
    
    dc_dt = diff(ctot_10)/dt;

    if abs((log(2)-log(GMD))/log(Gsigma)) < 1.65 || abs((log(2)-log(GMD))/log(Gsigma)) > 2.58
        for i =indices
            if i==indices(1)
                tmp_stack = 100*ones(1,max_bin-idx_10+1);  %disp(size(tmp_stack))
            else
                c_tmp = 100*ones(1,max_bin-idx_10+1); 
                tmp_stack = [tmp_stack,c_tmp];
            end
        end
        tmp_tot= 100*ones(1,nscan_final);
        tmp_stack = [tmp_stack,tmp_tot];
    else
        for i =indices
            c_tmp=cnum_est(idx_10:max_bin,i)/max_c(i); %disp(size(c_tmp))
            if i==indices(1)
                tmp_stack = c_tmp';
            else
                tmp_stack = [tmp_stack,c_tmp'];
            end
        end
        tmp_tot = ctot_10'/SourceConcPeak;
        tmp_stack = [tmp_stack,tmp_tot];
    end
    t_IUQ(j,:) = tmp_stack; 
end
end
%% Function for Ultrafine Particle Dynamics Modeling
%======================================================================================================================
%> @details     This function, Aerosol_Model, simulates aerosol dynamics in an indoor environment. It includes various physical processes like deposition, air change, and particle emission. The process involves:

%>              (1) Initial Setup and Data Reading:
%>                  - Reads experimental data and initializes variables for the simulation.
%>                  - Calculates deposition efficiency and effective particle loss due to air changes.

%>              (2) Source Emission Rate Calculation:
%>                  - Computes the emission rate of particles from the source, incorporating factors like air change rate and particle size distribution.

%>              (3) Particle Concentration Simulation:
%>                  - Simulates the time-varying concentration of particles using a time marching approach.
%>                  - Incorporates processes such as coagulation, deposition, and air changes.

%>              (4) Result Compilation:
%>                  - Stores the particle concentration at each time step.
%>                  - Calculates the total number concentration for the entire simulation period.
%======================================================================================================================
%%
function [ctot_coag,Qcnum_coag,scnum] = Aerosol_Model(control_volume,radius,Dradius,volume,Vfrac,exp_data,cnum_BG,nbin,beta,a,b,f,GMD,Gsigma,delay)
% *** Read experiment data sets
dt = exp_data(1); %[sec]
control_volume= control_volume*1e+6;
nscan_ini = 1;
nscan_final = exp_data(2);
nscan_offsource = exp_data(3);
SourceConcPeak_10 = exp_data(4);
ACh = exp_data(5); % [h-1]

% *** Background particle distribution
% *** Calculation of deposition efficiency [h-1]
DepoEff = exp(a * log(2 * radius * 1e+7) + b); %disp(size(DepoEff ))
DepoMean = mean(DepoEff);

% *** Particle loss term due to air change rate and deposition
EffDepo = 1 - exp(-DepoEff/3600*dt);
EffACh =  1 - exp(-ACh/3600*dt);

% *** Source emission rate(*1e+12/min)
emission_time = (nscan_offsource-nscan_ini+1)*dt/60 ; %disp(emission_time) %[min]
SourceConcPeak = SourceConcPeak_10*f;
prior_ss = (ACh+DepoMean)/60*SourceConcPeak/(1-exp(-(ACh+DepoMean)/60*emission_time)); % total num of particles [min-1]

SourceStrength = prior_ss;
Source_nt = SourceStrength/60; % Total number of ptcl released per unit time [sec-1]
Source_cmd= GMD*1e-7; %[nm]
Source_lnGMD= log(Source_cmd);
Source_lnSig= log(Gsigma);

% Vectorized calculations
tmp1 = Source_nt * (2 * Dradius) ./ (2 * radius) / sqrt(2 * pi) / Source_lnSig; %disp(size(tmp1))
tmp2 = (log(2 * radius)-Source_lnGMD).^2 / (2 * Source_lnSig^2); %disp(size(tmp2))
scnum = tmp1 .* exp(-tmp2) * dt;

% *** Time varying total number concentration
IDepoACh =1;
% Activation function test
% *** Allocate arrays
Pcnum= zeros(nbin,nscan_final);
% *** Initialize
Pcnum(:,1) = cnum_BG; %initial cnum
cnum = cnum_BG;
cnum0 = cnum;

[H,~]= Act_func_Heavi(delay,nscan_offsource,nscan_final);

% *** Time marching of particle distributions
% Assuming that you have already initialized a parallel pool using parpool.
% Time marching of particle distributions using parallel for loop
for n_iter = nscan_ini:nscan_final-1 % 2~ nscan_final繹먮슣?
    cnum0 = cnum0 + H(n_iter)*scnum;
    % *** Coagualtion
    cnum_coag = particle_balance(Vfrac,volume,cnum0,beta,nbin,dt);
    % *** Decoupled loss terms due to deposition and ACh --> 0: off, 1: on
    if(IDepoACh == 1)
        [cnum,~,~] = loss_term(EffDepo,EffACh,cnum_coag,cnum0);
    end
    cnum0 = cnum;
    % *** Store Concentration at each scan
    Pcnum(:,n_iter+1) = cnum;
end % end of time marching
Qcnum_coag = Pcnum;
ctot_coag = Cal_TotNum_new(Qcnum_coag,1,nbin,1,nscan_final);
end
%% Function for Initial Setup of Indoor Ultrafine Particle Model
%======================================================================================================================
%> @details     The PI_IndoorUFP_InitialSetting function configures initial settings for ultrafine particle modeling based on the source of emission. It involves:

%>              (1) Setting Up Environment and Particle Properties:
%>                  - Defines ambient conditions like air temperature and pressure.
%>                  - Sets up particle properties including density, molecular weight, and diffusion coefficients.

%>              (2) Emission Source Configuration:
%>                  - Based on the emission source, selects appropriate experimental data.
%>                  - Configures control volume, reference parameters, and data files for the chosen source.

%>              (3) Reading and Processing Experimental Data:
%>                  - Reads data from CSV files and extracts necessary parameters.
%>                  - Calculates particle size distributions and initializes particle property arrays.

%>              (4) Kernel Evaluation and Background Concentration Setup:
%>                  - Evaluates coagulation kernels for different particle sizes.
%>                  - Sets up background particle concentrations and total number concentrations.

%>              (5) Finalizing and Saving Initial Settings:
%>                  - Compiles all settings and experimental data into structured formats.
%======================================================================================================================
function BPM_IndoorUFP_InitialSetting (num_source)
rng(100,'twister')
disp('start')

% Ambient and particle properties setup
global amb_cond ptcl_prop ptcl_basic exp_data
amb_cond.Tair= 300; % Air temperature [K]
amb_cond.Pair= 1; % Air pressure [atm]

ptcl_prop.avognum = 6.02252e+23;
ptcl_prop.boltg   = 1.38054e-16; % [g cm^2 s^(-2) K^(-1)]
ptcl_prop.Wair    = 28.966;
ptcl_prop.Rgas    = 0.08206;
ptcl_prop.CONSMU  = 1.8325e-04 * (296.16 + 120);
ptcl_prop.WTMOLEC = 18.02;
ptcl_prop.rhoptcl = 1;   % density of water (fig 15.8 of FAM)
ptcl_prop.muair    = (ptcl_prop.CONSMU / (amb_cond.Tair + 120)) * (amb_cond.Tair / 296.16)^(1.5);
ptcl_prop.rhoair   = amb_cond.Pair * ptcl_prop.Wair * 0.001 / (ptcl_prop.Rgas * amb_cond.Tair);
ptcl_prop.nuair    = ptcl_prop.muair / ptcl_prop.rhoair;
ptcl_prop.MTV_air   = sqrt(8 * ptcl_prop.boltg * amb_cond.Tair * ptcl_prop.avognum / (pi * ptcl_prop.Wair));
ptcl_prop.MFP_air = 2 * ptcl_prop.nuair / ptcl_prop.MTV_air;

% Selecting emission source and corresponding settings
if num_source ==1
    filename_read_csv = 'Gasstove(GAS).csv';
    filename_setting_mat = 'Gas_setting.mat';
    control_volume = 340;
    Num_hamaker =2;
    optpar_ref = [0.84 2.49 2.29 4.65 1.55 0]; % Refer to Rim et al.(2016) 
elseif  num_source==2
    filename_read_csv = 'Candle(CAND).csv';
    filename_setting_mat = 'Candle_setting.mat';
    control_volume = 34;
    Num_hamaker =6;
    optpar_ref = [1.22 1.53 7.96 3.86 1.38 0]; % Refer to Rim et al.(2016) 
end

% Reading experimental data from CSV files
exp_inf = csvread(filename_read_csv,1,0,[1 0 5 0]);
ptcl_inf = csvread(filename_read_csv,7,0);
exp_data.nbinexp = exp_inf(1,1);
exp_data.nscanexp = exp_inf(2,1);
Nbinexp = exp_data.nbinexp;
Nscanexp = exp_data.nscanexp;
exp_data.ioffsource = exp_inf(3,1)-1;
exp_data.ACh = exp_inf(4,1);
exp_data.dtmin = exp_inf(5,1);
D_exp = ptcl_inf(:,1);
rmin = D_exp(1)/2 *1e-7;
exp_data.dt = exp_data.dtmin * 60; %[sec]
exp_data.xcnum = ptcl_inf(:,2:exp_data.nscanexp+1);

% Setting up particle diameter and related calculations
vratavg = 0;
for k=1:Nbinexp-1
    vratloc = (D_exp(k+1)/D_exp(k))^3;
    vratavg = vratavg+ vratloc;
end
vratavg = vratavg/(Nbinexp-1);
vrat = vratavg;

% **** Particle Basics on each sizebin : setting particle size distribution ****
NBIN =  length(D_exp);
ptcl_basic.nbin = NBIN;
ptcl_basic.radius = zeros(NBIN,1);
ptcl_basic.volume = zeros(NBIN,1);
ptcl_basic.vol_Smol=zeros(NBIN,1);
ptcl_basic.volumewidth=zeros(NBIN,1);
ptcl_basic.dlogdp = zeros(NBIN,1);
ptcl_basic.Dradius = zeros(NBIN,1);
ptcl_basic.Dvolume = zeros(NBIN,1);
ptcl_basic.volumelow = zeros(NBIN,1);
ptcl_basic.volumehigh = zeros(NBIN,1);
ptcl_basic.ptclMFPnm = zeros(NBIN,1);
ptcl_basic.diamum =  zeros(NBIN,1);
ptcl_basic.Cdiff = zeros(NBIN,1);
ptcl_basic.ptclAvgTvel = zeros(NBIN,1);

adv = 2 * (vrat- 1) / (vrat + 1);
vrlow = (2 / (1+ vrat))^(1/3);
vrhigh = vrlow * vrat^(1/3);

for i = 1:ptcl_basic.nbin
    ptcl_basic.radius(i) = rmin*vrat^((i-1)/3);
    ptcl_basic.volume(i) = (4*pi/3) * ptcl_basic.radius(i)^3;
    ptcl_basic.volumelow(i) = 2*ptcl_basic.volume(i)/(1+vrat);
    ptcl_basic.volumehigh(i) = ptcl_basic.volumelow(i)*vrat;
    ptcl_basic.volumewidth(i) = ptcl_basic.volumehigh(i) - ptcl_basic.volumelow(i);
    radiuslow = ptcl_basic.radius(i) * vrlow;
    radiushigh = ptcl_basic.radius(i) * vrhigh;
    ptcl_basic.Dradius(i) = radiushigh - radiuslow;
    ptcl_basic.Dvolume(i) = ptcl_basic.volume(i) * adv;
    ptcl_basic.diamum(i) =  ptcl_basic.radius(i) * (2e+4);
    ptcl_basic.dlogdp(i) =  log10(2*radiushigh) - log10(2*radiuslow);
    ptcl_basic.Cdiff(i)  = diffusion(ptcl_basic.radius(i));
    ptcl_basic.ptclAvgTvel(i) = sqrt(8 * ptcl_prop.boltg * amb_cond.Tair / (pi * ptcl_prop.rhoptcl * ptcl_basic.volume(i)));% [cm/s]
    ptcl_basic.ptclMFPnm(i) = 8 * ptcl_basic.Cdiff(i) / (pi * ptcl_basic.ptclAvgTvel(i)) * 1e+7; % [nm]
end

idx_10= find(ptcl_basic.radius<5*1e-7,1,'last'); idx_10 = max(idx_10)+1;
exp_10= find(D_exp>=10,1,'first');

% Kernel evaluation for particle interactions
hamaker = 20:20:200;
BETA = Kernel_evaluation(NBIN,hamaker);
Vfrac = volume_fraction(ptcl_basic.volume,NBIN);

% Additional calculations and data processing for forward modeling  
cnum_exp = zeros(NBIN,Nscanexp);
cnum_exp(1:Nbinexp,:) = exp_data.xcnum;
r_data = D_exp(exp_10:length(D_exp))/2*1e-7;
cnumini= cnum_exp(:,1);
NumCurve = length(cnum_exp(1,:));

ctot_all = Cal_TotNum_new(cnum_exp,1,Nbinexp,1,NumCurve); 
ctot_10 = Cal_TotNum_new(cnum_exp,idx_10,Nbinexp,1,NumCurve); 

[SourceConcPeak_all,idx_max_all] = max(ctot_all);
[SourceConcPeak_10,idx_max] = max(ctot_10);
ratio = SourceConcPeak_all/SourceConcPeak_10;
ExpInfo = zeros(5,1);
ExpInfo(1) = exp_data.dt; %[sec]
ExpInfo(2) = exp_data.nscanexp;
ExpInfo(3) = exp_data.ioffsource;
ExpInfo(4)= SourceConcPeak_10;
ExpInfo(5) = exp_data.ACh ;
nscan_offsource = ExpInfo(3)+1;


if nscan_offsource ~= idx_max_all
    disp('diff')
    expdata_10(4)= ctot_10(nscan_offsource);
else
    expdata_10(4)= SourceConcPeak_10;
end

cnum_exp = zeros(NBIN,Nscanexp);
cnum_exp(1:Nbinexp,:) = exp_data.xcnum;
r_data = D_exp(exp_10:length(D_exp))/2*1e-7;
cnumini= cnum_exp(:,1);
cnumini_10 = zeros(NBIN,1);
cnumini_10(idx_10:NBIN,1)= cnum_exp(idx_10:NBIN,1);
NumCurve = length(cnum_exp(1,:));
nscan_offsource = ExpInfo(3)+1;
ctot_all = Cal_TotNum_new(cnum_exp,1,Nbinexp,1,NumCurve); % True value
ctot_10 = Cal_TotNum_new(cnum_exp,idx_10,Nbinexp,1,NumCurve); % True value
time = 0:2.5:2.5*(ExpInfo(2)-1);

% Finalizing and saving the initial settings
disp('Initial setting finished')
save(filename_setting_mat );
end 
%%
function ctot =Cal_TotNum_new(pcnum,bin_10,nbin,nscan_ini,nscan_fin)
ctot = zeros(nscan_fin,1);
for ns= nscan_ini:nscan_fin
    for k=bin_10:nbin
        ctot(ns) = ctot(ns) + pcnum(k,ns);
    end
end
end
%%
function [f, t] = Act_func_Heavi(a, ioffsource, nscan_fin)
t = 0:nscan_fin;  
f = zeros(1, nscan_fin);  

b = a + ioffsource; 


active_intervals = (t(1:end-1) < b) & (t(2:end) > a);


f(active_intervals) = 1;

if mod(a, 1) ~= 0
    a_index = find(t >= a, 1) - 1;
    f(a_index) = (t(a_index+1) - a); 
end

if mod(b, 1) ~= 0
    b_index = find(t < b, 1, 'last');
    f(b_index) = (b - t(b_index)); 
end
end
%%
function [cnum_new, DDepo, DACh] = loss_term(EffDepo, EffACh, cnum, cnum0)
% Vectorized calculations
total_loss = (EffDepo + EffACh).* cnum0;
cnum_new = cnum - total_loss;

% Ensure cnum_new doesn't go negative
cnum_new = max(cnum_new, 0);

% Calculate losses
DDepo = EffDepo .* cnum0;
DACh = EffACh .* cnum0;

% Adjust for bins where cnum_new is zero
zero_bins = cnum_new == 0;
ratio = EffDepo(zero_bins) ./ (EffDepo(zero_bins) + EffACh);
DDepo(zero_bins) = cnum(zero_bins) .* ratio;
DACh(zero_bins) = cnum(zero_bins) .* (1 - ratio);
end

%%
function [cnum]= particle_balance(Vfrac,volume,cnum0,beta,nbin,dt)
% Particle balance - eq 15.14
cnum = cnum0;
cvol0= cnum0.*volume;
cvol = cvol0;
for k= 1: nbin
    term1 = 0;
    for j = 1: k
        for i = 1: k-1
            term1 = term1 + Vfrac(i,j,k) * beta(i,j) * cvol(i) * cnum0(j);
        end
    end
    term2 = 0;
    for j = 1: nbin
        term2 = term2 + (1 - Vfrac(k,j,k)) * beta(k,j) * cnum0(j);
    end
    cvol(k) = (cvol0(k) + dt * term1) / (1 + dt * term2);
    cnum(k) = cvol(k) / volume(k);
end
end
%% Volume fraction
function Vfrac= volume_fraction(Volume,nbin)
volij = zeros(nbin,nbin);
Vfrac = zeros(nbin,nbin,nbin);

for i = 1 : nbin
    for j = i : nbin
        volij(i,j) = Volume(i) + Volume(j);
        volij(j,i) = volij(i,j);
    end
end

for k = 1: nbin
    for j = 1: nbin
        for i = 1: nbin
            kp = min(k+1,nbin);
            km = max(k-1,1);
            if(k< nbin && volij(i,j) >= Volume(k)&& volij(i,j) < Volume(kp))
                Vfrac(i,j,k) = (Volume(k+1) - volij(i,j))/(Volume(k+1) - Volume(k))*Volume(k)/volij(i,j);
            elseif (k > 1 && volij(i,j) > Volume(km) && volij(i,j) < Volume(k))
                Vfrac(i,j,k) = 1 - Vfrac(i,j,k-1);
            elseif(k==nbin && volij(i,j) >= Volume(k))
                Vfrac(i,j,k) = 1;
            else
                Vfrac(i,j,k) = 0;
            end
        end
    end
end
end
%%
function Cdiff = diffusion(radius)
global ptcl_prop amb_cond

nbin= length(radius);
Cdiff = zeros(nbin,1);
for i =1:nbin
    airKnud = ptcl_prop .MFP_air/ radius(i);
    BPM = 1 + airKnud * (1.249 + 0.42 *exp(-0.87/airKnud));
    Cdiff(i)  = ptcl_prop .boltg * amb_cond.Tair * BPM  / (6* pi * radius(i) * ptcl_prop .muair);
end
end
%% Kernel evaluation
function beta = Kernel_evaluation(nbin,hamaker)
beta= zeros(nbin,nbin,2+2*length(hamaker));

BrKernel = Brownian(nbin);  beta_Br = BrKernel; % Brownian (sphere)
FrKernel = Fractal(nbin);   beta_Fr = FrKernel;  % Fractals (fractal)
for i=1:length(hamaker)
    [vdWaal,FrvdWaal] = vdW_Visc(nbin,hamaker(i));
    beta_vdW(:,:,i) = BrKernel.*vdWaal; % Brownian + van der Waals (sphere)
    beta_vdWFr(:,:,i) = FrKernel.*FrvdWaal; % Brownian + van der Waals (sphere)
end

beta(:,:,1)= beta_Br;
beta(:,:,2:length(hamaker)+1) = beta_vdW;
beta(:,:,length(hamaker)+2) = beta_Fr;
beta(:,:,length(hamaker)+3:2*length(hamaker)+2) = beta_vdWFr;
end
%% Kernel evaluation: Brownian
function BrKernel = Brownian(nbin)
global ptcl_basic
BrKernel = zeros(nbin,nbin);
% Calculate delta, diffusion coefficient
delta= deltafn(ptcl_basic.radius,ptcl_basic.ptclMFPnm);
Cdiff= diffusion(ptcl_basic.radius);
for i = 1: nbin
    for j = i: nbin
        rsum = ptcl_basic.radius(i) + ptcl_basic.radius(j); % [cm]
        Dsum = Cdiff(i) + Cdiff(j); %[cm^2/s]
        delsum2 = delta(i)^2 + delta(j)^2; %[cm^2]
        rdelta = rsum / (rsum + sqrt(delsum2)); % x
        avgTvel = sqrt(ptcl_basic.ptclAvgTvel(i)^2 + ptcl_basic.ptclAvgTvel(j)^2);%[cm/s]
        Dfactor = 4 * Dsum / (rsum * avgTvel); % x
        % Transition : Fuchs interpolation formula
        BrKernel(i,j) = 4 * pi * rsum * Dsum / (rdelta + Dfactor); %[cm^3/(no*s)]: eq.(15.33)
        %BrKernel(i,j) = 4 * pi * rsum * Dsum; % continuum
        BrKernel(j,i) = BrKernel(i,j);
    end
end
end
%% Kernel evaluation: Fractal (eq 15.52)
function FrKernel = Fractal(nbin)

% Calculate Brownian collision kernel modified for fractal geometry
FrKernel = zeros(nbin);

[radi_Fr,Cdiff_modi,ptclAvgTvel_modi,delta_modi] = Fractalbasics(nbin);

for i= 1: nbin
    for j = i: nbin
        rsum = radi_Fr(i) + radi_Fr(j); % [cm]
        Dsum = Cdiff_modi(i) + Cdiff_modi(j); %[cm^2/s]
        delsum2 = delta_modi(i)^2 + delta_modi(j)^2; %[cm^2]
        rdelta = rsum / (rsum + sqrt(delsum2)); % x
        avgTvel = sqrt(ptclAvgTvel_modi(i)^2 + ptclAvgTvel_modi(j)^2);%[cm/s]
        Dfactor = 4.0 * Dsum / (rsum * avgTvel); % x
        FrKernel(i,j) = 4 * pi * rsum * Dsum / (rdelta + Dfactor); %eq. 15.33
        FrKernel(j,i) = FrKernel(i,j);
    end
end
end
%% Kernel evaluation: van der Waals/ viscous force (ch 15.6.5)
function [vdWaal,FrvdWaal] = vdW_Visc(nbin,hamaker)
global ptcl_basic
FrvdWaal = zeros(nbin);
vdWaal = zeros(nbin) ;
A_H = hamaker;
QCONS1   = A_H/ 6;

[radi_Fr,Cdiff_modi,ptclAvgTvel_modi,~] = Fractalbasics(nbin);
Cdiff= diffusion(ptcl_basic.radius);

for i = 1: nbin
    for j = i: nbin
        rinm = ptcl_basic.radius(i) * 1e+7; 
        rjnm = ptcl_basic.radius(j)  * 1e+7; 
        [vdWfactor_con,vdWfactor_free] = vdWintegral(rinm,rjnm,QCONS1);
        
        % sphere
        rsum = ptcl_basic.radius(i) + ptcl_basic.radius(j); 
        Dsum = Cdiff(i) + Cdiff(j); 
        avgTvel = sqrt(ptcl_basic.ptclAvgTvel(i)^2 + ptcl_basic.ptclAvgTvel(j)^2);
        Dfactor = 4 * Dsum / (rsum * avgTvel); 
        vdWaal(i,j) = vdWfactor_con * (1+Dfactor) / (1+(vdWfactor_con/vdWfactor_free)*Dfactor); 
        vdWaal(j,i) = vdWaal(i,j);
        
        %fractal
        rsum = radi_Fr(i) + radi_Fr(j); 
        Dsum = Cdiff_modi(i) + Cdiff_modi(j);
        avgTvel = sqrt(ptclAvgTvel_modi(i)^2 + ptclAvgTvel_modi(j)^2); 
        Dfactor = 4 * Dsum / (rsum * avgTvel); 
        FrvdWaal(i,j) = vdWfactor_con * (1+Dfactor) / (1 +(vdWfactor_con/vdWfactor_free)*Dfactor); 
        FrvdWaal(j,i) = FrvdWaal(i,j);
    end
end
end
%% Fractal basics : for van der Waals
function [r_Fr,Cdiff_modi,ptclAvgTvel_modi,delta_modi] = Fractalbasics(nbin)
global amb_cond ptcl_prop ptcl_basic
r_Fr = zeros(nbin,1);
Cdiff_modi = zeros(nbin,1);
ptclAvgTvel_modi = zeros(nbin,1);
ptclMFPnm_modi = zeros(nbin,1);
delta_modi = zeros(nbin,1);

for i = 1: nbin
    Dim_Fr = 1.7;
    Dim_sphere = 3;
    D_spherule= 24e-7; % [cm] : typical diameter for diesel soot , 24nm
    
    r_spherule = D_spherule/2; 
    V_spherule = 4*pi/3 * r_spherule^3; 
    V_local = 4*pi/3 * ptcl_basic.radius(i)^3;
    Num_spherule = V_local / V_spherule; 
    
    if ptcl_basic.radius(i) < r_spherule
        Dim_Fr = Dim_sphere;
    end
    
    r_Fr(i) = r_spherule * Num_spherule^(1/Dim_Fr); 
    
    % Mobility radius
    area1eq = Num_spherule^(2/3) ;
    area2eq = 1 + 2/3*(Num_spherule-1);
    area3eq = Dim_Fr * Num_spherule^(2/Dim_Fr)/3;
    max_areaeq = max(area1eq,min(area2eq,area3eq));
    r_areaeq  = r_spherule * sqrt(max_areaeq); 
    
    if ptcl_basic.radius(i)< r_spherule
        r_areaeq = ptcl_basic.radius(i);
    end
    
    r_Mob1 = r_Fr(i)/(log(r_Fr(i)/r_spherule) + 1);
    r_Mob2 = r_Fr(i)*((Dim_Fr - 1)/2)^0.7;
    r_Mob  = max([r_Mob1,r_Mob2,r_areaeq]);
    
    if(ptcl_basic.radius(i)< r_spherule)
        r_Mob = ptcl_basic.radius(i);
    end
    
    Cdiff_modi(i) = diffusion(r_Mob);
    ptclAvgTvel_modi(i) = sqrt(8 * ptcl_prop.boltg * amb_cond.Tair / (pi*ptcl_prop.rhoptcl*V_local));% [cm/s]
    delta_modi(i) = deltafn(r_Mob,ptclMFPnm_modi(i));
end
end
%% Calculateed van der Waals collision kernel 
function [vdWfactor_con,vdWfactor_free] = vdWintegral(rinm,rjnm,QCONS1)
rplus  = rinm + rjnm;
rplus2 = rplus^2 ;
rminus  = rinm - rjnm;
rminus2 = rminus^2 ;
rproduct = rinm * rjnm;

numdx = 10000; % number of integration

rval = rplus;
vdWfactor_con  = 0;
vdWfactor_free = 0;
drval     = rplus /500;
drvald2   = drval/2;

for k = 1: numdx % numerical integration
    rval = rval + drval;
    rcursq   = rval^2;
    rvalp1   = rval + drvald2;
    rcursqp1 = rvalp1^2;
    rvalm1   = rval - drvald2;
    rcursqm1 = rvalm1^2;
    term_vdW  = rproduct / (rplus * (rval - rinm - rjnm));
    DinfdDr = 1 + 2.6*rproduct*sqrt(term_vdW)/rplus2 + term_vdW;  
    %  van der Waals interaction potential E 
    term = 2*rproduct/(rcursq - rplus2) + 2*rproduct/(rcursq - rminus2)+ log((rcursq - rplus2)/(rcursq - rminus2));
    termp1 = 2*rproduct/(rcursqp1 - rplus2) + 2*rproduct/(rcursqp1 - rminus2)+ log((rcursqp1 - rplus2)/(rcursqp1 - rminus2));
    termm1 = 2*rproduct / (rcursqm1 - rplus2) + 2*rproduct / (rcursqm1 - rminus2)+ log((rcursqm1 - rplus2)/(rcursqm1 - rminus2));
    %   first order derivative of E w.r.t 'r'
    D_Ep1   = (termp1 - term)/ drvald2;
    D_Em1   = (term   - termm1)/drvald2;
    D_Eavg = (termp1 - termm1)/drval; % for exponential part
    %  second order derivative of E w.r.t 'r'
    D2_E   = (D_Ep1 - D_Em1)/drvald2;
    %  exponential term in eq(15.44)
    exp_term  = exp(-QCONS1 * term);
    %  van der Waals + viscous force coagulation enhancement factor for continuum regime
    vdWfactor_con = vdWfactor_con + rplus * DinfdDr * exp_term * drval / rcursq; 
    term1    = D_Eavg + rval * D2_E;
    term2    = exp(QCONS1 * (0.5 * rval* D_Eavg + term));
    term3    = rcursq * drval;
    %   van der Waals coagulation enhancement factor for free molecular regime
    vdWfactor_free = vdWfactor_free + QCONS1*term1*term2*term3 / (2*rplus2 );
end % end of integral
vdWfactor_con = 1 / vdWfactor_con; 
end
%%
function [vdWfactor_con,vdWfactor_free] = vdWintegral_new(rinm,rjnm,QCONS1)
rplus  = rinm + rjnm;
rplus2 = rplus^2 ;
rminus  = rinm - rjnm;
rminus2 = rminus^2 ;
rproduct = rinm * rjnm;

% For Gaussian quadrature
[nodes, weights] = lgwt(10000, rplus, 2*rplus); % lgwt is a function that gives nodes and weights for given number of points and limits
vdWfactor_con  = 0;
vdWfactor_free = 0;

for k = 1:length(nodes)
    rval = nodes(k);
    wval = weights(k);
    rcursq = rval^2;
    
    term_vdW  = rproduct / (rplus * (rval - rinm - rjnm));
    DinfdDr = 1 + 2.6*rproduct*sqrt(term_vdW)/rplus2 + term_vdW;
    
    term = 2*rproduct/(rcursq - rplus2) + 2*rproduct/(rcursq - rminus2)+ log((rcursq - rplus2)/(rcursq - rminus2));
    
    % Computing derivatives using finite differences
    epsilon = 1e-6;
    term_perturb = 2*rproduct/((rcursq+epsilon) - rplus2) + 2*rproduct/((rcursq+epsilon) - rminus2) + log(((rcursq+epsilon) - rplus2)/((rcursq+epsilon) - rminus2));
    D_Eavg = (term_perturb - term) / epsilon;
    D2_E = (term_perturb - 2*term + term) / epsilon^2;
    
    exp_term  = exp(-QCONS1 * term);
    
    vdWfactor_con = vdWfactor_con + wval * rplus * DinfdDr * exp_term / rcursq;
    
    term1    = D_Eavg + rval * D2_E;
    term2    = exp(QCONS1 * (0.5 * rval* D_Eavg + term));
    term3    = rcursq;
    
    vdWfactor_free = vdWfactor_free + wval * QCONS1*term1*term2*term3 / (2*rplus2);
end

vdWfactor_con = 1 / vdWfactor_con;
end
%%
function [x, w]=lgwt(N, a, b)
% lgwt computes the Legendre-Gauss nodes x and weights w on the interval [a,b] with N nodes.
% The function returns x and w such that the integral from a to b of a given function f
% can be approximated as: integral(f) = sum( w .* f(x) ).

% Calculation of the zeros of the Legendre polynomial
x = cos(pi*(4*(1:N)-1)./(4*N+2))';
P0 = zeros(N,1);
P1 = ones(N,1);

for k = 2:N
    P2 = ((2*k-1)*x.*P1-(k-1)*P0)/k;
    P0 = P1;
    P1 = P2;
end
x = x - P1./(2*P2);
x = 0.5*((b-a)*x + a+b);

% Calculation of the weights
w = 2./((1-x.^2).*P2.^2)*(b-a)/2;
end
%%
function delta = deltafn(radius,ptclMFPnm)
nbin= length(radius);
delta = zeros(nbin,1);

for i = 1: nbin
    ptclMFPcm = ptclMFPnm(i) * 1e-7;
    term1 = 2 * radius(i) + ptclMFPcm; % [cm]
    term2 = 4 * radius(i)^2 + ptclMFPcm^2 ; %[cm^2]
    term3 = radius(i) * ptclMFPcm; %[cm^2]
    delta(i) = (term1^3 - term2^1.5)/(6.0*term3) - 2* radius(i); % [cm]
end
end