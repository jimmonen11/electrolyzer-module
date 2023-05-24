%% File to Load all relevant parameters for Simulink
clc; clear all; close all

%% SOEC Load Parameters

Ppa = 5e5; %pressure we are at in Pa, 5 bara

%Molar mass of species
MM_H2O = 18;
MM_O2 = 32;
MM_H2 = 2;
MM_N2 = 28;

w_H2_SP = 0.01220; %10 mol% setpoint into the SOECs
rec_percent_cathode_init = 0.1220; %initial rec% of cathode stream


w_O2_SP = 0.2756; %25 mol% setpoint into the SOECs
rec_percent_anode_init = 0.1686; %initial rec% of anode stream

%Loads SS values for time intialization
load('concSS.mat')

T_cell_init = 1073; %Initial temperature of cell, K

n_cell = 8; %discretization of HTSE cell
num_cells = 605711/100; %number of HTSE cells

%Initial guesses for Fsolve of current and cell potential
x0 = zeros((5*n_cell),1);
x0(1:n_cell) = 1.28; %V, electrical potential of cell - made to be constant be eqs
x0(n_cell+1:2*n_cell) = 400; %A/m^2, local cathode current density
x0(2*n_cell+1:3*n_cell) = 400;  %A/m^2, local anode current density
x0(3*n_cell+1:4*n_cell) = 0.2; %V, cathode activation overpotential
x0(4*n_cell+1:end) = 0.2; %V, anode activation overpotential

%Fundamental constansts
F = 96485.33; %C/mol, Faraday's constant
sig = 5.6703e-8; %W/(m^2*K^4), stefan boltzman constant
R = 8.314; %J/molK, universal gas constant

%Geometry of HTSE cell
h_C = 0.001; %m, height of cathode channel
h_A = 0.001; %m, height of anode channel
h_I = 7.5e-4; %m, height of interconnect

%85.5 cm^2 active cell area - design fit to give current denisty of 1500 mA/cm^2
L_cell = 0.1; %m, length of cell
W_cell = 0.0855; %m, width of cell
dx_cell = L_cell/n_cell; %m, length of discretized node of cell

%thicknesses of cell
%cathode supported from Strategic Analysis
t_C = 300e-6; %m, thickness of cathode
t_E = 10e-6; %m, thickness of electrode
t_A = 30e-6; %m  thickness of anode

h_S = t_C + t_E + t_A; %m, height of solid stucture (anode, cathode, and electrolyte)

%electric conductivities of cathode, electrolye and anode
sig_C = 80e3; %1/ohm*m
sig_E = 1.9; %1/ohm*m, this is around correct value, model replaces with temp dependent
sig_A = 31047.9; %1/ohm*m, A detailed kinetic model for the reduction of oxygen on LSCF-GDC

%Diffusion constants
Deff_C = 3.93e-5; %m^2/s, effective diffusion coef in cathode
Deff_A = 3.99e-6; %m^2/s, effective diffusion coef in anode

%Solid Structure and interface Properties
eps_S = 0.8; %solid structure emissivity
eps_I = 0.1; %interconnect emissivity

Cp_S = 500; %J/(kg*K), solid structure heat capacity
Cp_I = 500; %J/(kg*K), interface heat capacity

lam_S = 2; %J/(m*s*K), solid stucture conduction coef
lam_I = 25; %J/(m*s*K), interconnect conduction coef

lam_C = 0.194 ;%J/(m*s*K), cathode stream conduction coef
lam_A = 0.069; %J/(m*s*K), anode stream conduction coef

rho_S = 5900; %kg/m^3, solid structure density
rho_I = 8000; %kg/m^3, interface density

%Nusselt number assumed constant
Nu_C = 3.09;
Nu_A = 3.09;

dh_C = (2*W_cell*h_C)/(W_cell+h_C); %m, hydralic diameter of cathode
dh_A = (2*W_cell*h_A)/(W_cell+h_A); %m, hydralic diameter of anode

%Convective HTF in cathode
k_C = (Nu_C*lam_C)/dh_C; %J/(m^2*s*K)

%Convective HTF in anode
k_A = (Nu_A*lam_A)/dh_A; %J/(m^2*s*K)

% %Activation energy and pre-exponential factor of cathode and anode
A_C = 6.54e11; %1/ohm
A_A = 2.35e11; %1/ohm

E_C = 1.4e5; %J/mol
%LSCF-GDC https://reader.elsevier.com/reader/sd/pii/S0013468620300116?token=3D0AD6BD2EAA0B3C5678BBE0393E3966A5CFAE723A7381724A73BB7FA73A7C85AD05D362E587BCCDCC383C842AA1C495&originRegion=us-east-1&originCreation=20230125164343
E_A = 1.21e5; %J/mol

alph = 0.5; %transfer coefficient
eta_elec = 0.99; %electrical heating efficiency

%% HX Load Parameters
%Cathode High Heat Disc. and Length
nhx_wHR = 15; %crossflow passes
nhx_aHR = 31; %m

%Anode High Heat Disc. and Length
Lhx_wHR = 1.25; %crossflow passes
Lhx_aHR = 4.95; %m

%Heat exchanger dimensions and discretization
ttube = 0.002794; %thickness of sch 40 1/2" plain carbon steel tube
ktube = 60.5; %W/m*K, conduction coeff for plain carbon steel at 300K 
              %tabulated in Incropera pg. 930
              
%Heat exchanger design and tube layout
Dtubei = 0.015748; %m, inner diameter of each tube, 1/2" sch 40 tubes
Dtubeo = Dtubei +2*ttube; %m, outer diameter of each tube, 1/2" sch 40 tubes
SlD =  1.5; 
StD = 2; 

Sl = SlD*Dtubeo; %longitudal pitch of tubes, distance back (in direction of flow)
St = StD*Dtubeo; %traverse pitch of tubes, (distance down/up, normal to direction of flow)

ntubes = 500;

%will have to change below parameters to match geometry
ntubel = sqrt(ntubes);
ntubet = sqrt(ntubes);

%longitdal length tubes encompass using longitudal height for each tube plus extra in front and back
Ltubel = Sl* (ntubel + 2);

%transverse length tubes encompass encompass using transverse height for each tube plus extra on top and bottom
Ltubet = St* (ntubet + 2);

%Anode High Heat Recovery
tubenodewidth_aHR = Lhx_aHR/nhx_aHR; %m, length of each tube subject to crossflow
Acsflow_aHR = Ltubet * tubenodewidth_aHR; %cross sectional area of the flow as it passes tubes for each node, estimated as a square
V_tubes_aHR = pi*(Dtubei/2)^2*tubenodewidth_aHR*ntubes; %m^3, volume air takes up in each node
V_shell_aHR = (Acsflow_aHR*Ltubel) - (pi*(Dtubeo/2)^2*tubenodewidth_aHR*ntubes); %m^3, volume HTF takes up in each node
Anode_aHR = pi*Dtubeo*tubenodewidth_aHR*ntubes; %surface area where heat transfer occurs on each node

%Cathode Heat Recovery
tubenodewidth_wHR = Lhx_wHR/nhx_wHR; %m, length of each tube subject to crossflow
Acsflow_wHR = Ltubet * tubenodewidth_wHR; %cross sectional area of the flow as it passes tubes for each node, estimated as a square
V_tubes_wHR = pi*(Dtubei/2)^2*tubenodewidth_wHR*ntubes; %m^3, volume air takes up in each node
V_shell_wHR = (Acsflow_wHR*Ltubel) - (pi*(Dtubeo/2)^2*tubenodewidth_wHR*ntubes); %m^3, volume HTF takes up in each node
Anode_wHR = pi*Dtubeo*tubenodewidth_wHR*ntubes; %surface area where heat transfer occurs on each node

[C1, m] = tubebankpar(SlD, StD);

%Heat Exchanger Initialization Parameters
T_wSHout = 700 + 273.13; %K
T_aSHout = 700 + 273.13; %K



