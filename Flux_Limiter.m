function [Flux] = Flux_Limiter(phi, phi_in, u, dx, method)
%Calculates flux using flux limiter techiniques

% Calculate the discrete RHS for the PDE
%   - phi   : the solution at time level n
%   - method: the flux limiter to use

phi = phi'; %Transpose phi

n = length(phi);

% Set up indexing for stencil operations below
pt    = 2:1:n-1;
left  = pt - 1;  % left stencil coefficient indices
right = pt + 1; % right stencil coefficient indices

left2  = left - 1; %shift for left2

% Here we know that u > 0 so we don't have to look at its sign.
% This simplifies a lot of logic

% Upwind fluxes at left and right faces
f_upw_left  = u*phi(left);
f_upw_right = u*phi(pt);

% Central fluxes at left and right faces
f_c_left  = u/2 * ( phi(pt) + phi(left) );
f_c_right = u/2 * ( phi(pt) + phi(right) );

%----------- Limiting function -----------
% avoid divide-by-zero in calculation of r
r_right = ones(1,length(phi)-2);
itest =  abs(phi(right)-phi(pt));
i = find(itest>=1e-10); %indexes where we have number greater than zero

for i = 1:length(r_right)
    r_right(i) = (phi(pt(i)) - phi(left(i)) ) / ( phi(right(i)) - phi(pt(i)));
end
limit_right = limiter_value(r_right,method);


r_left = ones(1,length(phi)-2);
itest =  abs(phi(pt)-phi(left));
i = find(itest>=1e-10); %indexes where we have number greater than zero
i = i(2:end);
r_left(1) = ( phi(left(1)) - phi_in)  / ( phi(pt(1)) - phi(left(1)) ); %whats going on here

for i = 2:length(r_left)
    r_left(i) = ( phi(left(i)) - phi(left2(i)) ) / ( phi(pt(i)) - phi(left(i)) );
end

limit_left = limiter_value(r_left,method);

% Limited fluxes at left and right faces
F_left  = f_upw_left - limit_left.*(f_upw_left - f_c_left);
F_right = f_upw_right - limit_right.*(f_upw_right - f_c_right);

%%
%extra logic to account for very first node!
%Flux to the right
r_r1 = (phi(1) - phi_in)/(phi(2) - phi(1)); 

F_right1 = u*phi(1)-limiter_value(r_r1, method) * (u*phi(1) - u*(phi(1) + phi(2))/2);
F_right = [F_right1 F_right];

%Flux to the left
r_l1 = 0; %r is zero because we phi_i-1 - phi_i-2 = 0

F_left1 = u*phi_in-limiter_value(r_l1, method) * (u*phi_in - u*(phi_in + phi(1))/2);
F_left = [F_left1 F_left];

%%
%extra logic to account for very last node!
%Flux to the right
r_rn = 100; %will get really small denominator going to infinit

F_rightn = u*phi(n)-limiter_value(r_rn, method) * (u*phi(n) - u*(phi(n) + phi(n))/2); %phi(n+1) = phi(n)
F_right = [F_right F_rightn];

%Flux to the left
r_ln = (phi(n-1) - phi(n-2))/(phi(n) - phi(n-1)); %r is zero because we phi_i-1 - phi_i-2 = 0

F_leftn = u*phi(n-1)-limiter_value(r_ln, method) * (u*phi(n-1) - u*(phi(n-1) + phi(n))/2);
F_left = [F_left F_leftn];

%%

Flux =  (F_right - F_left ) / dx;
