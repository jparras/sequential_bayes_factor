function[dec,n]=ump(x,p0,alpha,nt)

% Function that computes the UMP in a Bernoulli sequence
% H0: p=p0
% H1: p>p0
% alpha: significance level
% nt: length of test (fixed)

% Obtain c threshold
c=find(binocdf(0:nt,nt,p0)>=1-alpha,1); % s>=c is the threshold

% Test vector x

s = sum(x);

if s>=c 
    dec=1; %Decide H1
else %Max length of test achieved
    dec=0;
end

n = length(x);