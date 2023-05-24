function [rho, mu, k, Cp] = prop_mix(rho_1, rho_2, mu_1, mu_2, k_1, k_2, C_1, C_2, w_1, w_2)

%Function that gives weighted average of properties

rho = w_1*rho_1 + w_2*rho_2;
mu = w_1*mu_1 + w_2*mu_2;
k = w_1*k_1 + w_2*k_2;
Cp = w_1*C_1 + w_2*C_2;