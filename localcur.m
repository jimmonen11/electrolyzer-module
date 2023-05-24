function Fsolve = localcur(x, params)
%Formulates system of equations to be solved so can get voltage and local current simulataneously as they are linked
%%

%First seperate all the parameters included in params
n = params(end); %always include at the end for clarity!!!

c_H2O = params(1:n);
c_H2 = params(n+1:2*n);
c_O2 = params(2*n+1:3*n);
c_N2 = params(3*n+1:4*n);
T_A = params(4*n+1:5*n);

T_S = params(5*n+1:6*n);
javg = params(6*n+1);

F = params(6*n+2);
R = params(6*n+3);

t_C = params(6*n+4);
t_E = params(6*n+5);
t_A = params(6*n+6);

sig_C = params(6*n+7);
sig_A = params(6*n+8);
sig_E = params(6*n+9); %use temp varying instead, but here if want to not do that

A_C = params(6*n+10);
A_A = params(6*n+11);

E_C = params(6*n+12);
E_A = params(6*n+13);

Deff_C = params(6*n+14);
Deff_A = params(6*n+15);

alph = params(6*n+16);

P = params(6*n+17); %pressure

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%F is a vector with all 5*n equations
Fsolve = ones(length(x),1);

%split up x for clarity in calculations
U = x(1:n); %V, electrical potential of cell - made to be constant be eqs
j_C = x(n+1:2*n); %A/m^2, local cathode current density
j_A = x(2*n+1:3*n);  %A/m^2, local anode current density
n_Cact = x(3*n+1:4*n); %V, cathode activation overpotential
n_Aact = x(4*n+1:end); %V, anode activation overpotential


%=========================================================================%
%Triple phase boundary (TPB)

c_H2_TPB = zeros(size(c_H2O));
c_H2O_TPB = zeros(size(c_H2O)); 
c_O2_TPB = zeros(size(c_H2O));


%triple phase boundary equations
for i = 1:n
    c_H2_TPB(i) = c_H2(i) + t_C/(2*F*Deff_C)*j_C(i);
    c_H2O_TPB(i) = c_H2O(i) - t_C/(2*F*Deff_C)*j_C(i);
    c_O2_TPB(i) = c_O2(i) + c_N2(i) - c_N2(i)*exp((-t_A*j_C(i))/(4*F*Deff_A*(c_O2(i)+c_N2(i))));    
end

%=========================================================================%
%Overpotentials and cell potential

Unot = zeros(size(c_H2O));
Urev = zeros(size(c_H2O));
R_ohm = zeros(size(c_H2O));
n_ohm = zeros(size(c_H2O));
n_Cconc = zeros(size(c_H2O));
n_Aconc =  zeros(size(c_H2O));

p_O2 =  zeros(size(c_H2O));
G =  zeros(size(c_H2O));

for i = 1:n
    %Interpolation spreadsheet from CoolProp
    G(i) =  -56.019*T_S(i) + 248198; %J/mol H2O, Gibbs free energy of reaction 
    
    Unot(i) = G(i)/(2*F); %V, standard potential
    
    p_O2(i) = c_O2(i)/(c_O2(i) + c_N2(i))*P*1e-5;%  101325e-5  partial pressue
    
    %Reversible potential
    Urev(i) = Unot(i) + (R*T_S(i))/(2*F)*log((c_H2(i)*p_O2(i)^0.5)/c_H2O(i)); %V
    
    %Cell reisitivity
    R_ohm(i) = t_C/sig_C + t_E/(33.4e3*exp(-10.3e3/T_S(i))) + t_A/sig_A; %ohms*m^2

    %Cell overpotential
    n_ohm(i) = j_C(i)*R_ohm(i);% V
    
    %Cathode concentration overpotential
    n_Cconc(i) = (R*T_S(i))/(2*F)*log((c_H2_TPB(i)*c_H2O(i))/(c_H2(i)*c_H2O_TPB(i))); %V
    
    %Anode concentration overpotential
    n_Aconc(i) = (R*T_S(i))/(4*F)*log((c_O2_TPB(i)*T_S(i))/(c_O2(i)*T_A(i))); %V
end

%=========================================================================%
j0_C = zeros(size(c_H2O));
j0_A = zeros(size(c_H2O));
jC_curr = zeros(size(c_H2O));
jA_curr = zeros(size(c_H2O));
n_contact = zeros(size(c_H2O));


%loop to calculate local current density
for i = 1:n
    
    %cathode exchange current density (arrhenius eq like)
    j0_C(i) = (R*T_S(i)*A_C)/(2*F)*exp(-E_C/(R*T_S(i)));
    
    %local cathode current density (modified Butler-Volmer eq.)
    jC_curr(i) = j0_C(i)*((c_H2_TPB(i)/c_H2(i))*exp((2*(1-alph)*F*n_Cact(i))/(R*T_S(i)))-...
    (c_H2O_TPB(i)/c_H2O(i))*exp((-2*alph*F*n_Cact(i))/(R*T_S(i)))); 

    %anode exchange current density (arrhenius eq like)
    j0_A(i) = (R*T_S(i)*A_A)/(2*F)*exp(-E_A/(R*T_S(i)));

    %local andoe current density (modified Butler-Volmer eq.)
    jA_curr(i) = j0_A(i)*(exp((2*(1-alph)*F*n_Aact(i))/(R*T_S(i)))-...
    (c_O2_TPB(i)/c_O2(i))*exp((-2*alph*F*n_Aact(i))/(R*T_S(i)))); 

    %contact resistance
    n_contact(i) = 0.4075e-5*jC_curr(i);
  end


for i = 1:n
    %Potential - n equations where U
    Fsolve(i) = -U(i) + n_ohm(i) + Urev(i) + n_Cconc(i) + n_Aconc(i) + n_Cact(i) + n_Aact(i) + n_contact(i);
    %local currents and guesses
    Fsolve(n+i) = (j_C(i) - jC_curr(i));
    Fsolve(2*n+i) = (j_A(i) - jA_curr(i));
    %cathode current is same as anode current
    Fsolve(3*n+i) = (j_C(i) - j_A(i));
    
    if i == n
        Fsolve(5*n) = (javg - mean(j_A)); 
    else
    %forces U to be the same for all spots on cell
    Fsolve(4*n+i) = U(i) - U(i+1);
    end
end
    