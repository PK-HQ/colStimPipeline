
[11 13 20 26 35]
% colormap
mkrCond={'o','v','^'};
cmap=fireice;
cmap=cmap([22 1+10 43-10],:);
%newcolors = {'black','black','blue','blue','red','red'};
%colororder(newcolors)
for condSet=1
  %{
    conds2extract=[conds0(condSet,:), conds90(condSet,:)]; %[-5 -12 -30] [5 12 30]
    deltaOrientation=gaborCon(conds2extract);
    percentHorReport=100-vertReportPrctData(conds2extract);
    percentVertReport=vertReportPrctData(conds2extract);
    nTrials=completedData(conds2extract);

    % Maximum likelihood fit with a generalized gaussian fun
    %[muFit, betaFit, sigmaFit, Lapse, Likelihood] = FitPsyML(deltaOrientation,percentVertReport,percentHorReport,'plot')
    %}
    
    deltaOrientation=[12 15 18 20 25 30 40];%[0 8 10 12 14 25];
    deltaOrientation=[-deltaOrientation,deltaOrientation];
    percentVertReport=[30 10 20 0 0 0 0 50 60 80 90 90 70 100];%[50 40 20 20 0 0 0 30 60 60 80 90 100 100];%[70 30 20 40 10 0 80 40 80 90 90 100];
    [~,betaFit,muFit,sigmaFit] = fitPsyDualCDFv2(deltaOrientation,percentVertReport,'20230719',cmap(condSet,:),mkrCond{condSet});
    delta(condSet)=betaFit;
    mu(condSet)=muFit;
    sigma(condSet)=sigmaFit;
end