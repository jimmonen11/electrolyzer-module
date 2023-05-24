function eqsolver(block)
%This is a special fucntion to avoid using intrepreted Matlab functions when calling fsolve make's code run way faster

setup(block);

%endfunction

%%
function setup(block)
n = block.DialogPrm(1).Data; %spacial discretization is first parameter

% Register number of ports
block.NumInputPorts  = 12;
block.NumOutputPorts = 1;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
for i = 1:block.NumInputPorts-1
    block.InputPort(i).Dimensions        = n;
end

%last input dimensions for javg
block.InputPort(block.NumInputPorts).Dimensions  = 1;

% Override output port properties
block.OutputPort(1).Dimensions       = 5*n;

% Register parameters
block.NumDialogPrms     = 17; %number of parameters to input in box on simulink
block.SampleTimes = [0 0];
block.SimStateCompliance = 'DefaultSimState';

%% -----------------------------------------------------------------
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup
%%
function Outputs(block)

n = block.DialogPrm(1).Data;

x0 = zeros((5*n),1);

x0(1:n) = block.InputPort(7).Data; %V, electrical potential of cell - made to be constant be eqs
x0(n+1:2*n) =block.InputPort(8).Data; %A/m^2, local cathode current density
x0(2*n+1:3*n) = block.InputPort(9).Data;  %A/m^2, local anode current density
x0(3*n+1:4*n) = block.InputPort(10).Data; %V, cathode activation overpotential
x0(4*n+1:end) = block.InputPort(11).Data; %V, anode activation overpotential

%parameters to pass to localcur.m
params = zeros(6*n+18,1);

params(1:n) = block.InputPort(1).Data; %c_H2O;
params(n+1:2*n) = block.InputPort(2).Data; %c_H2;
params(2*n+1:3*n) = block.InputPort(3).Data; %c_O2;
params(3*n+1:4*n) = block.InputPort(4).Data; %c_N2;
params(4*n+1:5*n) = block.InputPort(5).Data; %T_A;
params(5*n+1:6*n) = block.InputPort(6).Data; %T_S;
params(6*n+1) = block.InputPort(12).Data; %javg, current density requested

for k = 2:17
    params(6*n+k) = block.DialogPrm(k).Data; %assign all parameter values to pass
end

params(end) = n; %last parameter is the spatial discretization

%If statement handles when applied current is less than 1 (essentially zero)
%Gives U as Urev with all other values as zero as there should be no overpotential
if params(6*n+1) < 1
    
    G =  -56.019*mean(params(5*n+1:6*n)) + 248198; %J/mol H2O, Gibbs free energy of reaction (see interpolation spreadsheet) 
    Unot = G/(2*params(6*n+2)); %V, standard potential
    p_O2 = (mean(params(2*n+1:3*n))/(mean(params(2*n+1:3*n)) + mean(params(3*n+1:4*n))))*params(6*n+17)*1e-5;%  101325e-5 
    Urev = Unot + (params(6*n+3)*mean(params(5*n+1:6*n)))/(2*params(6*n+2))*log((mean(params(n+1:2*n))*p_O2^0.5)/mean(params(1:n))); %V
    
    block.OutputPort(1).Data(1:n) = Urev; %V, electrical potential of cell - made to be constant be eqs
    block.OutputPort(1).Data(n+1:2*n) =  0; %A/m^2, local cathode current density
    block.OutputPort(1).Data(2*n+1:3*n) = 0;  %A/m^2, local anode current density
    block.OutputPort(1).Data(3*n+1:4*n) = 0; %V, cathode activation overpotential
    block.OutputPort(1).Data(4*n+1:end) = 0; %V, anode activation overpotential
else %for all other cases solve the system of equations with fsolve 
    
    %Can play around with tolerance to make faster sims
    %options = optimset('Display','off', 'TolFun', 5e-2);
    options = optimset('Display','off', 'TolFun', 1e-4);

    %solve system of eqs to find voltage and local currents
    [x] = fsolve(@(x)localcur(x,params),x0, options);
    block.OutputPort(1).Data(1:n) = x(1:n); %V, electrical potential of cell - made to be constant be eqs
    block.OutputPort(1).Data(n+1:2*n) = x(n+1:2*n); %A/m^2, local cathode current density
    block.OutputPort(1).Data(2*n+1:3*n) = x(2*n+1:3*n);  %A/m^2, local anode current density
    block.OutputPort(1).Data(3*n+1:4*n) = x(3*n+1:4*n); %V, cathode activation overpotential
    block.OutputPort(1).Data(4*n+1:end) = x(4*n+1:end); %V, anode activation overpotential
end


%end Outputs

%%
function Terminate(block)

%end Terminate

