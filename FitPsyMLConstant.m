function [alpha beta bias lapse likelihood version threshold] = FitPsyMLConstant(varargin)
% Maximum Likelihood fit (assuming no lapse Rate)

%% Parameters
% Input
% deltaOrientation = independent variable (e.g stimulus conditions, Orientation = [-45, 15, 5,15,45] deg)
% percentCWreport = number percentCWreport responses (monkey says left) 
% percentCCWreport = number percentCCWreport responses (monkey says right)

% Output
% alpha = param gen gauss (Threshold)
% beta = param gen gauss (Shape)
% bias = param gen gauss (Bias, point of subjective equality)
% likelihood = likelihood of the fit (i.e. goodness of fit)
% version = 0.1.1

%% Calling function:
% [alpha beta likelihood] = FitPsyML(deltaOrientation,percentCWreport,percentCCWreport)
% [alpha beta likelihood] = FitPsyML(deltaOrientation,percentCWreport,percentCCWreport, 'plot') ; % plot the fit in
% current graph
% [alpha beta likelihood] = FitPsyML(deltaOrientation,percentCWreport,percentCCWreport,'FixBeta',betaC) % If you want to
% fix the shape parameter

% To display the fit: 
%  y = gengaussx([alpha beta bias],deltaOrientation)

%------------------------------------------------
% by Giacomo Benvenuti 2018 <giacomox@gmail.com>
%------------------------------------------------
version = '0.1.1';

% Plot params
yLimits=[0,100];

% Set minimization parameters
% Limits minimixation
% options=optimset('MaxFunEvals',500000000) ;

% varargin{1}  = [5 10 20 30 45]
% varargin{2}  =[ 1 5 10 15 20]
% varargin{3}  = [20 15 10 5 1]
% Get inputs ----------
deltaOrientation = varargin{1} ;
percentCWreport = varargin{2} ;
percentCCWreport = varargin{3} ;
nTrials = varargin{4};
sessionStr= varargin{5};
verboseTitle=varargin{9};
xLimits=varargin{11};
if ~isempty(varargin{12})
    cmap=varargin{12};
end
if ~isempty(varargin{13})
    mkr=varargin{13};
else
    mkr='.';
end

p = find(strcmp(varargin,'Param0') == 1,1);
if ~isempty(p)
% use the fix beta value provided by the user
a0 = varargin{p+1} ; 
else
a0=[13 1.2 0 0.1 0];%[5 1  0.1 0];% initial param values 
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
%param=bads(@(param) Cost_Function(param,deltaOrientation,percentCWreport,percentCCWreport),a0,LB,UB) ; 
% Run minimization
% if exist('bads')==2 % bads is a minimization process
%  param=bads(@(param) Cost_Function(param,deltaOrientation,percentCWreport,percentCCWreport),a0,LB,UB) ; 
% elseif exist('fminsearchbnd')==2
% [param]=fminsearchbnd(@(param) Cost_Function(param,deltaOrientation,percentCWreport,percentCCWreport),a0,LB,UB) ; %,options
% else

%main function
%[param]=fminsearch(@(param) Cost_Function(param,deltaOrientation,percentCWreport,percentCCWreport),a0,[30,1,20,0]) ; %,options
% y = (1-2*lapse) * normcdf(  1/2 * sign(deltaOrientation) .* (abs(deltaOrientation)/alpha).^beta  ) + lapse ;
%alpha = param(1)  ;
%beta = param(2) ;
%u = param(3); %bias
%lapse = param(4) ;
%Yconstant=param(5);

paramLB=[0 0 0 -.3 -.3];%[0 1.25 -1 -.3 -.3];
paramUB=[30 3 5 .3 .3];%[30 2.75 1 .3 .3];
[param]=fminsearchbnd(@(param) Cost_Function(param,deltaOrientation,percentCWreport,percentCCWreport),a0,...
    paramLB,...
    paramUB);

%[param]=fminsearchbnd(@(param) Cost_Function(param,deltaOrientation,percentCWreport,percentCCWreport),a0,[13.5 .95 -7 -.1],[16 1.1 7 .1]);%[12 .9 -45 -2],[15 1.1 45 2]) ; %,options
%end

alpha = param(1)
beta = param(2)
bias = param(3)
lapse = param(4)
yConstant=param(5);
likelihood  = Cost_Function([alpha beta bias lapse yConstant],deltaOrientation,percentCWreport,percentCCWreport);%Cost_Function([alpha beta bias 0],deltaOrientation,percentCWreport,percentCCWreport);


if PlotFit

   threshold=PlotFitFun([alpha beta bias lapse yConstant],likelihood,deltaOrientation,percentCWreport,nTrials,...
       xLimits,yLimits,sessionStr,verboseTitle,cmap,mkr);
end    

end

function threshold=PlotFitFun(param,likelihood,deltaOrientation,percentCWreport,nTrials,xLimits,yLimits,sessionStr,...
    verboseTitle,cmap,mkr)
     % --- Psychometric fit ---
     hold on
     rangeDeltaOrientation = linspace(xLimits(1), xLimits(2), 100); %linspace(min(deltaOrientation), max(deltaOrientation), 100) ;
     Yfit = gengaussx([param(1) param(2) param(3) param(4) param(5)],rangeDeltaOrientation) ;
     plot(rangeDeltaOrientation,Yfit*100,'Color',cmap,'LineWidth',2.5,'HandleVisibility','off');
     [~, indexOfMin] = min(abs(Yfit-.70));
     threshold=rangeDeltaOrientation(indexOfMin);

     offset=.1;
     line([-100 100] , [50 50],'Color',.75*[1 1 1], 'LineStyle','--','LineWidth',.25,'HandleVisibility','off')

     % --- scatter ---
     if sum(cmap)==0
         cmap=[1 1 1];
         alphalevel=1; %baseline in black
         mkrSz=100;
     else
         mkrSz=100;
         alphalevel=0.7;
     end
     h=scatter(deltaOrientation,percentCWreport,mkrSz,cmap,mkr,'filled','MarkerEdgeColor','k');
     set(h, 'MarkerFaceAlpha', alphalevel, 'MarkerEdgeAlpha', alphalevel) % set transparency level

     % display threshold
     if (xLimits(2)-xLimits(1))<10
        %addSkippedTicks(xLimits(1),xLimits(2),.5,'x')
     elseif (xLimits(2)-xLimits(1))>10
         %addSkippedTicks(xLimits(1),xLimits(2),10,'x')
     end
     %addSkippedTicks(yLimits(1)-10,yLimits(2),10,'y')
     %addSkippedTicks(0,35,5,'x')
     %addSkippedTicks(0,100,10,'y')
     %yticks([0:10:yLimits(2)])
     xlim([xLimits(1),xLimits(2)])
     ylim([yLimits(1),yLimits(2)])
     switch verboseTitle
         case {0}
              title({sessionStr,...
                 ['Threshold = ' num2str(threshold,'%.1f') ', Bias = ' num2str(param(3),'%.1f')]},...
                 'FontWeight','normal','FontName','Arial')
         case {1}
              title({sessionStr,...
                 ['Threshold = ' num2str(threshold,'%.1f') ', Bias = ' num2str(param(3),'%.1f') ', LL = ' num2str(likelihood,'%g')],...
                 ['a = ' num2str(param(1),'%g') ', B = ' num2str(param(2),'%g')  ', lse = ' num2str(param(4),'%g')]},...
                 'FontWeight','normal','FontName','Arial')
             %{
             title({sessionStr,...
                 ['Sz=',num2str(gaborSize),'\circ, E=',num2str(gaborEccentricity,'%.1f'),'\circ, SF=',num2str(gaborSF) 'cpd'],...
                 ['Threshold = ' num2str(param(1),'%.1f') '\circ' ', Bias = ' num2str(param(3),'%.1f')]},'FontWeight','normal')
             %}
     end
     %pubfig(h1,12,2,0.005);
     pubfig(h,12,2,0.005);
     %hold off
     %legend off
     hLeg=legend;
     hLeg.ItemTokenSize = [7,25];
     set(hLeg,'visible','off')
     %colormap
     %colormap(c(1:maxNTrials,:)); % apply new colormap
     %hColorbar=colorbar();
     %hColorbar.Limits=[0 1];
     %hColorbar.Ticks = linspace(0, 1, 5); %Create 8 ticks from zero to 1
     %hColorbar.TickLabels = num2cell((0:maxNTrials/4:maxNTrials+1)); 
    fprintf('%s: a=%g, B=%g, bias=%g, lapse=%g, LL=%g\n',sessionStr,param(1), param(2), param(3), param(4), likelihood)
end


function [Cost]= Cost_Function(param,deltaOrientation,percentCWreport,percentCCWreport)
% Generan gaussian fun (NO lapse)
F = gengaussx([param(1),param(2),param(3),param(4),param(5)],deltaOrientation);

% calculate probability for each data point
ct = 0;
for i =1:numel(F)
    if F(i) == 0 |  F(i) == 1
        err
    else
        ct = ct+1;
        tm(ct) =  percentCWreport(:,i).*log(F(i)) + percentCCWreport(:,i).*log(1-F(i))  ;
    end
end
LL =tm(:) ; % log likelihood
Cost = -2*sum(LL) ;

end

function y = gengaussx(param,x)
% generalized gaussian Psychometric Function
% param = [alpha beta lapse]
% y = (1-2*lapse) * normcdf(  1/2 * sign(deltaOrientation) .* (abs(deltaOrientation)/alpha).^beta  ) + lapse ;
alpha = param(1)  ;
beta = param(2) ;
u = param(3); %bias
lapse = param(4) ;
yConstant=param(5);

if alpha<0
    alpha =0;
end

y = ((1-2*lapse) * normcdf((0.5*sign(x-u) .* (abs(x-u)/alpha).^beta)) + lapse) + yConstant;
y(find(y<0.0001)) = 0.0001;
y(find(y> 0.9999)) = 0.9999;

end



% Versions
% V 0.0.2 added the bias to the gengauss