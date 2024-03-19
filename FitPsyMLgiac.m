function [alpha beta PSE Lapse Likelyhood Version] = FitPsyML(varargin)
% Maximum Likelihood fit (assuming no Lapse Rate)
% [alpha beta Likelyhood] = FitPsyML(x,A,B)
% [alpha beta Likelyhood] = FitPsyML(x,A,B, 'plot') ; % plot the fit in
% current graph
% [alpha beta Likelyhood] = FitPsyML(x,A,B,'FixBeta',betaC) % If you want to
% fix the shape parameter
% x = independent variable (e.g stimulus conditions, Orientation = [-45, 15, 5,15,45] deg)
% A = number A responses (monkey says left) 
% B = number B responses (monkey says right)
% alpha = param gen gauss (Threshold)
% beta = param gen gauss (Shape)
% PSE = param gen gauss (Bias, point of subjective equality)
% likelihood = likelihood of the fit (i.e. goodness of fit)
% 
% To display the fit: 
%  y = gengaussx([alpha beta PSE],x)
%
% Version = 0.1.1
%------------------------------------------------
% by Giacomo Benvenuti 2018 <giacomox@gmail.com>
%------------------------------------------------
Version = '0.1.1';
% Set minimization parameters
% Limits minimixation
% options=optimset('MaxFunEvals',500000000) ;

% varargin{1}  = [5 10 20 30 45]
% varargin{2}  =[ 1 5 10 15 20]
% varargin{3}  = [20 15 10 5 1]
% Get inputs ----------
x = varargin{1} ;
A = varargin{2} ;
B = varargin{3} ;


p = find(strcmp(varargin,'Param0') == 1,1);
if ~isempty(p)
% use the fix beta value provided by the user
a0 = varargin{p+1} ; 
else
a0=[5 1  0.1 0];% initial param values 
end


n = find(strcmp(varargin,'FixBeta') == 1,1);
if ~isempty(n)
    % use the fix beta value provided by the user
a0(2) = varargin{n+1} ; 
beta_min = a0(2) ;
beta_max = a0(2) ;
else
    % find optiamal beta
beta_min = 0 ;
beta_max = 5 ;
end

% Plot
clear n 
n = find(strcmp(varargin,'plot') == 1,1);
if ~isempty(n)
PlotFit = true ; 
else
PlotFit =false ; 
end


%------------------------------------
LB =[1 beta_min -45   ] ;
UB = [45 beta_max 45 ] ;
% Run minimization
% if exist('bads')==2 % bads is a minimization process
%  param=bads(@(param) Cost_Function(param,x,A,B),a0,LB,UB) ; 
% elseif exist('fminsearchbnd')==2
% [param]=fminsearchbnd(@(param) Cost_Function(param,x,A,B),a0,LB,UB) ; %,options
% else
[param]=fminsearch(@(param) Cost_Function(param,x,A,B),a0) ; %,options
    
%end

alpha = param(1);
beta = param(2);
PSE = param(3);
Lapse = param(4);
Likelyhood  = Cost_Function([alpha beta PSE 0],x,A,B);


if PlotFit
   PlotFitFun([alpha beta PSE 0],x)
end    

end

function  PlotFitFun(param,x)
     hold on
     X = linspace(min(x), max(x), 100) ;
     F = gengaussx([param(1) param(2) param(3) param(4)],X) ;
     plot(X,F);
end


function [Cost]= Cost_Function(param,x,A,B)
% Generan gaussian fun (NO LAPSE)
F = gengaussx([param(1),param(2),param(3),0],x);

% calculate probability for each data point
ct = 0;
for i =1:numel(F)
    if F(i) == 0 |  F(i) == 1
        err
    else
        ct = ct+1;
        tm(ct) =  A(:,i).*log(F(i)) + B(:,i).*log(1-F(i))  ;
    end
end
LL =tm(:) ; % log likelihood
Cost = -2*sum(LL) ;

end

function y = gengaussx(param,x)
% generalized gaussian Psychometric Function
% param = [alpha beta lapse]
% y = (1-2*Lapse) * normcdf(  1/2 * sign(x) .* (abs(x)/alpha).^beta  ) + Lapse ;
alpha = param(1)  ;
beta = param(2) ;
u = param(3); %PSE
Lapse = param(4) ;

if alpha<0
    alpha =0;
end

y = (1-2*Lapse) * normcdf(  [ 1/2 * sign(x-u) .* (abs(x-u)/alpha).^beta  ] ) + Lapse ;
y(find(y<0.0001)) = 0.0001;
y(find(y> 0.9999)) = 0.9999;

end



% Versions
% V 0.0.2 added the PSE to the gengauss