function [threshold,beta,mu,sigma,param,x,y] = fitPsySingleCDF(varargin)
%% README
% Fits psychometric function for contrast, 2 asymmetrical curves for x<mu and x>mu. Maximum Likelihood fit. Created
% by PK (2023), modified from Giac's script.

% --- Input ---
%       x = contrasts, signed according to orientation (0=-ve, 90=+ve)
%       y = % vertical reports

% --- Output ---
%       beta = vertical bias, scaling factor that modulates the height of normal cdf
%       mu = midpoint of cdf
%       sigma = spread of cdf
%       logLikelihood = goodness of fit

%% Extract data from varagin
x=varargin{1}; % gabor contrast, signed accroding to orientation (horz = -ve)
y=varargin{2}; % % vertical reports
sessionID=varargin{3}; % date as a str, for display only
markerColor=varargin{4}; % e.g. 'r'
markerType=varargin{5}; % e.g. 'o'
hAx=varargin{6};
%% Find optimal parameter values
% --- Initial values, lower and upper bounds --- 
smoothingFlag=0;
mu = 40/2; %is signed during fitting
sigma = 40/6;
beta = 0.5;
lapseLower=0;
lapseUpper=0;
initialParams=[mu     sigma    beta     lapseLower lapseUpper];
paramLB      =[10     10/6       0          0.0               0.0];
paramUB      =[100     40/6      1          00               0];
%% Run parameter optimization
sOptOptm = ...
    optimset('LargeScale','On', ...
    'TolFun',1e-8,'TolX',1e-8, ...
    'MaxFunEvals',1e4,'MaxIter',1e4);

% Adjusted Sigmoidal: [midpoint, sigma, beta]
%initialParams = [10, 5, 0.5];
%paramLB = [20   1    0];
%paramUB = [40   15   1];
%param = lsqcurvefit('dualNormCDF',initialParams,x,y,paramLB,paramUB,sOptOptm);

% --- Run parameter optimization --- 
opts = optimset('fminsearch');
opts.MaxFunEvals = 100000;
opts.MaxIter = 100000;
[param]=fminsearchbnd(@(param) costFunction(param,x,y,smoothingFlag),...
    initialParams,...
    paramLB,...
    paramUB);

%  --- Extract param values and likelihood of best fit --- 
mu = param(1);
sigma = param(2);
beta = param(3);
lapseLower = param(4);
lapseUpper = param(5);
logLikelihood  = costFunction([mu sigma beta lapseLower lapseUpper],x,y,smoothingFlag);%Cost_Function([alpha beta bias 0],betaOrientation,y,yInverse);

%% Plot fitted curve with scatter plot overlay
threshold=plotFitCurve([mu sigma beta lapseLower lapseUpper],logLikelihood,x,y,...
   sessionID,markerColor,markerType,hAx,smoothingFlag);

% Return biasing (beta)
beta = beta+lapseLower; %total vertical displacement of curve at midpoint

end

%% Cost function: to calculate likelihood of fit vs actual data
function cost = costFunction(param,x,y,smoothingFlag)
%sort y and inverse of y
yInverse=100-y;

% Fit curve to data
yFit = dualNormCDF([param(1),param(2),param(3),param(4),param(5)],x,smoothingFlag);

% calculate probability for each data point
ct = 0;
for i =1:numel(yFit)
    if yFit(i) == 0 | yFit(i) == 1
        err
    else
        ct = ct+1;
        LLcont(ct) =  y(:,i).*log(yFit(i)) + yInverse(:,i).*log(1-yFit(i))  ;
    end
end
LL =LLcont(:) ; % log likelihood
cost = -2*sum(LL) ;
end

% --- Curve fitting for data ---
function yFit = normCDF(param,x,smoothingFlag)
% Version: Feb 2023
% Psychometric Function for contrast, 2 asymmetrical curves for x<mu and x>mu
mu = param(1);
sigma = param(2);
beta = param(3);
lapseLower = param(4);
lapseUpper = param(5);

% Split data into 2 parts
xHorzConds=x(x<=0);%x(x<=0);
xVertConds=x(x>=-0);%x(x>=0);
xSplitConds=[xHorzConds;xVertConds];

for condSet=1:size(xSplitConds,1)
  conds=xSplitConds(condSet);
    
  % Fit 2 parts (<0 and >0), combine
  yFit=((beta) * normcdf(conds,mu,sigma))+lapseLower;
  
  % Combine to get single curve spanning x-range
  zeroPointExist=~isempty(find(xLower==0));
  switch zeroPointExist
      case {0}
          if smoothingFlag==1
              cont=nan(2, numel(unique([xLower, xUpper])));
              cont(1,1:numel(yFitLower))=yFitLower;
              cont(2,end-numel(yFitUpper)+1:end)=yFitUpper;
              yFit=mean(cont,'omitnan');
          else
              yFitLowerFinal=yFitLower(xLower<0);
              yFitUpperFinal=yFitUpper(xUpper>=0);
              yFit=[yFitLowerFinal yFitUpperFinal];
          end
      case {1}
          if smoothingFlag==1
              cont=nan(2, numel(unique([xLower, xUpper])));
              cont(1,1:numel(yFitLower))=yFitLower;
              cont(2,end-numel(yFitUpper)+1:end)=yFitUpper;
              yFit=mean(cont,'omitnan');
          else
              yFitLowerFinal=yFitLower(xLower<0);
              yFitUpperFinal=yFitUpper(xUpper>=0);
              yFit=[yFitLowerFinal yFitUpperFinal];
          end
  end
  cont=nan(2, numel(unique([xLower, xUpper])));
  cont(1,1:numel(yFitLower))=yFitLower;
  cont(2,end-numel(yFitUpper)+1:end)=yFitUpper;
  yFit=mean(cont,'omitnan');

  % clipping, carryover from Giac's script
  yFit(yFit<0.0001) = 0.0001;
  yFit(yFit> 0.9999) = 0.9999;
end

end

%% --- Plotting curve and scatter --- 
function threshold=plotFitCurve(param,logLikelihood,x,y,sessionStr,...
     markerColor,markerType,hAx,smoothingFlag)
     %axes(hAx)

     hold on
     
     % --- Plot curve --- 
     xFit = linspace(min(x),max(x),100);
     yFit = dualNormCDF(param,xFit,smoothingFlag);
     plot(xFit,yFit*100,'Color',markerColor,'LineWidth',2.5);hold on;
     [~, indexOfMin] = min(abs(yFit-.70));
     threshold=xFit(indexOfMin);
     
     % --- Scatter plot ---
     lw=3;%.25;
     baseline=sum(markerColor)==0;
     if baseline % Baseline in black
         markerColor=[1 1 1];
         alphalevel=1;
         mkrSz=300;%200;
     else % Optostim transparent
         mkrSz=250;%80;
         alphalevel=1;
     end
     h=scatter(x,y,mkrSz,markerColor,markerType,'filled','MarkerEdgeColor','k','HandleVisibility','off','LineWidth',lw);hold on;
     set(h, 'MarkerFaceAlpha', alphalevel, 'MarkerEdgeAlpha', alphalevel) % set transparency level
     set(h, 'MarkerFaceAlpha', alphalevel, 'MarkerEdgeAlpha', 1) % set transparency level
     
     % --- Chance guide lines --- 
     line([0 0],[0 100],'Color',.4*[1 1 1],'LineStyle','--','LineWidth',.25,'HandleVisibility','off')
     line([-100 100],[50 50],'Color',.4*[1 1 1], 'LineStyle','--','LineWidth',.25,'HandleVisibility','off')
     
     
     % Add ticks
     if max(x)-min(x)<10
        addSkippedTicks(min(x),max(x),.5,'x')
     else
         addSkippedTicks(-100,100,10,'x')
     end
     addSkippedTicks(0,100,100/8,'y')
     % Add limits
     xlim([min(x),max(x)])
     ylim([0,100])
     % Add title
     h=title(sessionStr,'FontWeight','normal','FontSize',14,'FontName','Arial','interpreter','latex');
     % Add verbose title (edit to add params 1-7)
     %{
     h=title({sessionStr,...
      ['${\mu}$ = ' num2str(param(1),'%.1f'),', ${\sigma}$ = ' num2str(param(2),'%.1f')...
      ', ${\beta}$ = ' num2str(param(3),'%.1f'),', $\mathcal{L}$ = ' num2str(logLikelihood,'%.1f')]},...
      'FontWeight','normal','FontSize',14,'FontName','Arial','interpreter','latex');
     %}
  
     %pubfig(h,12,2,0.005);

     % Add legend 
     hLeg=legend('location','northwest');
     hLeg.ItemTokenSize = [7,25];
     set(hLeg,'visible','off')
     hold on
end