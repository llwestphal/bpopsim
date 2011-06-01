% Jeffrey Barrick <jbarrick@msu.edu>
% Copyright (C) 2008-2009.

function out = establishment_probability(t, T, s, No)

%Set these parameters
%  t is the estimate that is being minimized
%  T number of generations per growth phase
%  No population size after dilution
%  r growth rate (Malthusian fitness) of ancestral population
%  The effective population size (Ne) is taken to be No r T


%t time in binary fission generations within growth cycle
r = log(2);
Nt  = No * exp(r*t);
%t
%Nt
D = exp(-r*T);
%D

%this equation is only valid when the mutant is still rare
%after one growth cycle, incorrect results for s > 0.2

%numerically solve equation for P(0)

%low_s_limit = 2*r*s*T
%PestO = fminbnd(@(x) abs(exp(-x*exp(r*s*T))-(1-x)), 0, 1)

x = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99 0.999, 0.9999, 0.99999, 0.999999 1];
fun_res = exp(-exp(r*s*T)*(1-x))-x;
j=size(x,2)-1;
for i=1:size(x,2)-1
    if ((fun_res(i) > 0) && (fun_res(i+1) < 0))
        j=i;
    end;
end;

%j
%x(j)
%fun_res(j)
%x(j+1)
%fun_res(j+1)
%[x(j) x(j+1)]

PextO = fzero(@(x) exp(-exp(r*s*T)*(1-x))-x, [x(j) x(j+1)], optimset('TolX',1e-6,'Display','notify'));
Pextt = exp(-D*exp(r*(1+s)*(T-t))*(1-PextO));

%and use this starting point to calculate P(t) by offset
Pt = 1 - Pextt;
%Pt
out = Nt .* Pt;