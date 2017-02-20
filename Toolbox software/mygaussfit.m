function [sigma,mu,A]=mygaussfit(x,y,h,optionStab)

%
% [sigma,mu,A]=mygaussfit(x,y)
% [sigma,mu,A]=mygaussfit(x,y,h)
%
% this function is doing fit to the function
% y=A * exp( -(x-mu)^2 / (2*sigma^2) )
%
% the fitting is been done by a polyfit
% the lan of the data.
%
% h is the threshold which is the fraction
% from the maximum y height that the data
% is been taken from.
% h should be a number between 0-1.
% if h have not been taken it is set to be 0.2
% as default.
%


%% threshold
if nargin==2, h=0.2; optionStab=1; end

%% cutting
ymax=max(y);
xnew=[];
ynew=[];
for n=1:length(x)
    if y(n)>ymax*h;
        xnew=[xnew,x(n)];
        ynew=[ynew,y(n)];
    end
end

%% fitting
ylog=log(ynew);
xlog=xnew;

if optionStab==1
    p=polyfit(xlog,ylog,2);
    A2=p(1);
    A1=p(2);
    A0=p(3);
    sigma=sqrt(-1/(2*A2));
    mu=A1*sigma^2;
    A=exp(A0+mu^2/(2*sigma^2));
end

if optionStab==2
    % option for numerical stability
    [p,s,s2]=polyfit(xlog,ylog,2);
    A2=p(1);
    A1=p(2);
    A0=p(3);
    A1=A1/s2(2)-2*A2*s2(1)/s2(2);
    A2=A2/s2(2)^2;
    sigma=sqrt(-1/(2*A2));
    mu=A1*sigma^2;
    A=exp(A0+mu^2/(2*sigma^2));
end
