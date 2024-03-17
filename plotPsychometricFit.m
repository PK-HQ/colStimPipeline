function [betaSorted,muSorted,sigmaSorted,contrastsSorted,percentVerticalSorted]=plotPsychometricFit(dataStruct,runData,taskName,plotFlags,hAx)
%% Plots psychometric curves, thresholds, biases across sessions per monkey
%analysis params
gaborEccentricityWanted=1.46;%1.46;%5;
gaborSizeWanted=1;%1;%1
gaborSFWanted=6;%6;%4

fprintf('Analysis: Z=%g, SF=%g, E=%g \n',gaborSizeWanted,gaborSFWanted,gaborEccentricityWanted)

%filenameStruct=generateFilenames(dataDir);
%TSFilename=filenameStruct.TS;

%cd(dataDir)
%% Plot raw psychometric curve per session
switch plotFlags(1)
    case {1}
        plotFlag='plot';
    case {0}
        plotFlag='';
end

% For each session get all runs


%pdf filename
%pdfFilename=[dataDir 'M' monkeyIDs{monkey} 'Z' num2str(gaborSizeWanted*100,'%05g') 'S' num2str(gaborSFWanted*100,'%05g')...
%   'E' num2str(gaborEccentricityWanted*100,'%05g') 'N' num2str(nSessions,'%03g')];

%load single run data from datafolder
%runData=load(TSFilename);

if isfield(runData.TS.Header.Conditions,'GaborSize')%'GaborSize')
    switch taskName
      case {'optostim'}
        % relevant DV
        condData=runData.TS.Header.Conditions;
        condIDs=1:numel(condData.TypeCond);
        stimType=condData.TypeCond; %0 1 3
        gaborOrt=condData.GaborOrt;
        projectedImageConds=runData.TS.Header.Conditions.ProjImg;
        visualControl=find(stimType==1);
        optoControl=contains(projectedImageConds,'L000');
        opto0=contains(projectedImageConds,'O00000');
        opto90=contains(projectedImageConds,'O09000');
        successData=runData.TS.Header.Outcomes.CountCondSuccess;
        completedData=runData.TS.Header.Outcomes.CountCondTotalValid;

        negative0=-1.*(stimType>0 & gaborOrt==0);
        tmp=negative0'.*-1;
        negative0(negative0==0)=1;


        gaborCon=condData.StimCon.*negative0;
        vertReportData=successData.*negative0' + completedData.*tmp;
        vertReportPrctData=(vertReportData./completedData)'.*100;


        % --- Sort them into blank/gabor 0 (noopto, 0opto, 90opto)/gabor 90 (noopto, 0opto, 90opto) ---
        blankCond=find(stimType==0);

        baselineNoOpto=~isempty(visualControl)==1 && isempty(optoControl)==1;

        switch baselineNoOpto
            case {1}
                visCond0control=find(stimType==1 & gaborOrt==0); %[2:4]
                visCond90control=find(stimType==1 & gaborOrt==90); %[5:7]
            case {0}
                visCond0control=find(optoControl==1 & gaborOrt==0);
                visCond90control=find(optoControl==1 & gaborOrt==90);
        end

        vis0_opto0=find(stimType==3 & gaborOrt==0 & opto0);
        vis0_opto90=find(stimType==3 & gaborOrt==0 & opto90);

        vis90_opto0=find(stimType==3 & gaborOrt==90 & opto0);
        vis90_opto90=find(stimType==3 & gaborOrt==90 & opto90);

        % --- compile into mat --- (change to con/incon optostim if need be)
        conds0=[visCond0control;vis0_opto0;vis0_opto90];
        conds90=[visCond90control;vis90_opto0;vis90_opto90];
        condsOptoControl=[visCond0control;visCond90control];
        condsOpto0=[vis0_opto0;vis90_opto0];
        condsOpto90=[vis0_opto90;vis90_opto90];

        vertReportOptoControl=vertReportData(condsOptoControl); %control
        vertReportOpto0=vertReportData(condsOpto0); %horz
        vertReportOpto90=vertReportData(condsOpto90); %vert

        completedDataOptoControl=completedData(condsOptoControl); %horz
        completedDataOpto0=completedData(condsOpto0); %horz
        completedDataOpto90=completedData(condsOpto90); %vert

        prct90opto=sum(vertReportOpto90(:))*100/sum(completedDataOpto90(:));
        prct0opto=sum(vertReportOpto0(:))*100/sum(completedDataOpto0(:));
        prctBaselineopto=sum(vertReportOptoControl(:))*100/sum(completedDataOptoControl(:));

        fprintf('\nH-0 optostim: %f %% vertical reports\n',...
            prct0opto)
        fprintf('Blank optostim: %f %% vertical reports\n',...
            prctBaselineopto)
        fprintf('V-90 optostim: %f %% vertical reports\n\n',...
            prct90opto)
        %
        
        %% Fit psychometric function here
        nOptoConds=size(conds0,1);
        mkrCond={'o','v','^'};
        % colormap
        cmap=fireice;
        if nOptoConds==3 %when optostim and baseline are in same blocks
            mkrCond={'o-','v-','^-'};
            cmap=cmap([22 1+10 43-10],:); %black, blue red
        elseif nOptoConds==2 %when optostim and baseline are in different blocks
            mkrCond={'v-','^-'};
            cmap=cmap([1+10 43-10],:); %blue & red
        elseif nOptoConds==1 %when optostim and baseline are in different blocks
            mkrCond={'o-'};
            cmap=cmap(22,:); %black
        end

        %init
        nContrasts=numel(unique(gaborCon(stimType==3)));
        beta=nan(1,3);
        mu=nan(1,3);
        sigma=nan(1,3);
        contrasts=nan(nContrasts,3);
        percentVertical=nan(nContrasts,3);
        for half=1:3
            axes(hAx(half))
            for condSet=1:nOptoConds %baseline, opto0, opto90
                conds2extract=[conds0(condSet,:), conds90(condSet,:)]; %[-5 -12 -30] [5 12 30]
                betaOrientation=gaborCon(conds2extract);
                percentHorReport=100-vertReportPrctData(conds2extract);
                percentVertReport=vertReportPrctData(conds2extract);
                nTrials=completedData(conds2extract);

                % Maximum likelihood fit with a generalized gaussian fun
                %[muFit, betaFit, sigmaFit, Lapse, Likelihood] = FitPsyML(betaOrientation,percentVertReport,percentHorReport,'plot')

                %[~,betaFit,muFit,sigmaFit] = fitPsyDualCDFv2(betaOrientation,percentVertReport,dataDir.date,cmap(condSet,:),mkrCond{condSet});
                [betaOrientation,percentVertReport]=mergeZeros(betaOrientation,percentVertReport);
                % define separate axes for each half
                [~,betaFit,muFit,sigmaFit,param,x,y] = fitPsyCDF(betaOrientation,percentVertReport,dataStruct.date,cmap(condSet,:),mkrCond{condSet},hAx,half);
                fprintf('mu=%.2f, sigma=%.2f, beta=%.2f, lapseLower=%.2f, lapseUpper=%.2f\n',param(1),param(2),param(3),param(4),param(5))
                hold on;
                beta(condSet)=betaFit;
                mu(condSet)=muFit;
                sigma(condSet)=sigmaFit;
                contrasts(:,condSet)=betaOrientation;
                percentVertical(:,condSet)=percentVertReport;            
                upFontSize(24,0.01)
            end
            % axis labels
            if half==1
                title('Horizontal stimuli');
                ylabel('Correct (%)')
                xlabel('Gabor contrast (%)')
                upFontSize(24,0.01)


            
                xlims=[-1 1] *100;%.*roundup(abs(unique(xlim)),10);
                xlim([0 round(max(betaOrientation)/50)*50])
                addSkippedTicks(0,100,12.5,'x')
                addSkippedTicks(0,100,12.5,'y')
                % legend
                [h,icons] = legend({'Baseline','Con. optostim (H-opto)','Incon. optostim (V-opto)'},'location','southeast');
                icons = findobj(icons, '-property', 'Marker', '-and', '-not', 'Marker', 'none');
                set(icons,'MarkerSize',15,'linewidth',3);
            elseif half==2
                title('Vertical stimuli');
                ylabel('Correct (%)')
                legend off 
                upFontSize(24,0.01)

            
                xlims=[-1 1] *100;%.*roundup(abs(unique(xlim)),10);
                xlim([0 round(max(betaOrientation)/50)*50])
                addSkippedTicks(0,100,12.5,'x')
                addSkippedTicks(0,100,12.5,'y')
                % legend
                [h,icons] = legend({'Baseline','Incon. optostim (H-opto)','Con. optostim (V-opto)'},'location','southeast');
                icons = findobj(icons, '-property', 'Marker', '-and', '-not', 'Marker', 'none');
                set(icons,'MarkerSize',15,'linewidth',3);
            elseif half==3
                title('Combined stimuli');
                ylabel('Correct (%)')
                legend off 
                upFontSize(24,0.01)
                
            
                xlims=[-1 1] *100;%.*roundup(abs(unique(xlim)),10);
                xlim([0 round(max(betaOrientation)/50)*50])
                addSkippedTicks(0,100,12.5,'x')
                addSkippedTicks(0,100,12.5,'y')
                % legend
                [h,icons] = legend({'Baseline','Incon. optostim','Con. optostim'},'location','southeast');
                icons = findobj(icons, '-property', 'Marker', '-and', '-not', 'Marker', 'none');
                set(icons,'MarkerSize',15,'linewidth',3);
                
            end
        end
        
        betaSorted=[beta(2) beta(1) beta(3)];
        muSorted=[mu(2) mu(1) mu(3)];
        sigmaSorted=[sigma(2) sigma(1) sigma(3)];
        contrastsSorted=[contrasts(:,2) contrasts(:,1) contrasts(:,3)];
        percentVerticalSorted=[percentVertical(:,2) percentVertical(:,1) percentVertical(:,3)];
        
        % Add title
        %{
        h=title({dataStruct.date,...
         ['${\beta}$: ' num2str(betaSorted(1),'%.2f') ', '  num2str(betaSorted(2),'%.2f') ', ' num2str(betaSorted(3),'%.2f')]},...
         'FontWeight','normal','FontSize',14,'FontName','Arial','interpreter','latex');
        %}
        %pubfig(h,12,2,0.005);


      case {'split'}
        % relevant DV
        condData=runData.TS.Header.Conditions;
        condIDs=1:numel(condData.TypeCond);
        stimType=condData.TypeCond; %0 1 3
        gaborOrt=condData.GaborOrt;
        projectedImageConds=runData.TS.Header.Conditions.ProjImg;
        visualControl=find(stimType==1);
        optoControl=contains(projectedImageConds,'L000');
        opto0=contains(projectedImageConds,'O00000');
        opto90=contains(projectedImageConds,'O09000');
        successData=runData.TS.Header.Outcomes.CountCondSuccess;
        completedData=runData.TS.Header.Outcomes.CountCondTotalValid;

        negative0=-1.*(stimType>0 & gaborOrt==0);
        tmp=negative0'.*-1;
        negative0(negative0==0)=1;


        gaborCon=condData.StimCon.*negative0;
        vertReportData=successData.*negative0' + completedData.*tmp;
        vertReportPrctData=(vertReportData./completedData)'.*100;


        % --- Sort them into blank/gabor 0 (noopto, 0opto, 90opto)/gabor 90 (noopto, 0opto, 90opto) ---
        blankCond=find(stimType==0);

        baselineNoOpto=~isempty(visualControl)==1 && isempty(optoControl)==1;

        switch baselineNoOpto
            case {1}
                visCond0control=find(stimType==1 & gaborOrt==0); %[2:4]
                visCond90control=find(stimType==1 & gaborOrt==90); %[5:7]
            case {0}
                visCond0control=find(optoControl==1 & gaborOrt==0);
                visCond90control=find(optoControl==1 & gaborOrt==90);
        end

        vis0_opto0=find(stimType==3 & gaborOrt==0 & opto0);
        vis0_opto90=find(stimType==3 & gaborOrt==0 & opto90);

        vis90_opto0=find(stimType==3 & gaborOrt==90 & opto0);
        vis90_opto90=find(stimType==3 & gaborOrt==90 & opto90);

        % --- compile into mat --- (change to con/incon optostim if need be)
        conds0=[visCond0control;vis0_opto0;vis0_opto90];
        conds90=[visCond90control;vis90_opto0;vis90_opto90];
        condsOptoControl=[visCond0control;visCond90control];
        condsOpto0=[vis0_opto0;vis90_opto0];
        condsOpto90=[vis0_opto90;vis90_opto90];

        vertReportOptoControl=vertReportData(condsOptoControl); %control
        vertReportOpto0=vertReportData(condsOpto0); %horz
        vertReportOpto90=vertReportData(condsOpto90); %vert

        completedDataOptoControl=completedData(condsOptoControl); %horz
        completedDataOpto0=completedData(condsOpto0); %horz
        completedDataOpto90=completedData(condsOpto90); %vert

        prct90opto=sum(vertReportOpto90(:))*100/sum(completedDataOpto90(:));
        prct0opto=sum(vertReportOpto0(:))*100/sum(completedDataOpto0(:));
        prctBaselineopto=sum(vertReportOptoControl(:))*100/sum(completedDataOptoControl(:));

        fprintf('\nH-0 optostim: %f %% vertical reports\n',...
            prct0opto)
        fprintf('Blank optostim: %f %% vertical reports\n',...
            prctBaselineopto)
        fprintf('V-90 optostim: %f %% vertical reports\n\n',...
            prct90opto)
        %
        
        %% Fit psychometric function here
        nOptoConds=size(conds0,1);
        mkrCond={'o','v','^'};
        % colormap
        cmap=fireice;
        if nOptoConds==3 %when optostim and baseline are in same blocks
            mkrCond={'o','v','^'};
            cmap=cmap([22 1+10 43-10],:); %black, blue red
        elseif nOptoConds==2 %when optostim and baseline are in different blocks
            mkrCond={'v','^'};
            cmap=cmap([1+10 43-10],:); %blue & red
        elseif nOptoConds==1 %when optostim and baseline are in different blocks
            mkrCond={'o'};
            cmap=cmap(22,:); %black
        end

        %init
        nContrasts=numel(unique(gaborCon(stimType==3)));
        beta=nan(1,3);
        mu=nan(1,3);
        sigma=nan(1,3);
        contrasts=nan(nContrasts,3);
        percentVertical=nan(nContrasts,3);
        for condSet=1:nOptoConds %baseline, opto0, opto90
            conds2extract=[conds0(condSet,:), conds90(condSet,:)]; %[-5 -12 -30] [5 12 30]
            betaOrientation=gaborCon(conds2extract);
            percentHorReport=100-vertReportPrctData(conds2extract);
            percentVertReport=vertReportPrctData(conds2extract);
            nTrials=completedData(conds2extract);

            % Maximum likelihood fit with a generalized gaussian fun
            %[muFit, betaFit, sigmaFit, Lapse, Likelihood] = FitPsyML(betaOrientation,percentVertReport,percentHorReport,'plot')
                       
            %[~,betaFit,muFit,sigmaFit] = fitPsyDualCDFv2(betaOrientation,percentVertReport,dataDir.date,cmap(condSet,:),mkrCond{condSet});
            [betaOrientation,percentVertReport]=mergeZeros(betaOrientation,percentVertReport);
            [~,betaFit,muFit,sigmaFit,param,x,y] = fitPsyCDF(betaOrientation,percentVertReport,dataStruct.date,cmap(condSet,:),mkrCond{condSet},hAx);

            
            fprintf('mu=%.2f, sigma=%.2f, beta=%.2f, lapseLower=%.2f, lapseUpper=%.2f\n',param(1),param(2),param(3),param(4),param(5))
            hold on;
            beta(condSet)=betaFit;
            mu(condSet)=muFit;
            sigma(condSet)=sigmaFit;
            contrasts(:,condSet)=betaOrientation;
            percentVertical(:,condSet)=percentVertReport;            
        end

        betaSorted=[beta(2) beta(1) beta(3)];
        muSorted=[mu(2) mu(1) mu(3)];
        sigmaSorted=[sigma(2) sigma(1) sigma(3)];
        contrastsSorted=[contrasts(:,2) contrasts(:,1) contrasts(:,3)];
        percentVerticalSorted=[percentVertical(:,2) percentVertical(:,1) percentVertical(:,3)];
        
        % Add title
        %{
        h=title({dataStruct.date,...
         ['${\beta}$: ' num2str(betaSorted(1),'%.2f') ', '  num2str(betaSorted(2),'%.2f') ', ' num2str(betaSorted(3),'%.2f')]},...
         'FontWeight','normal','FontSize',14,'FontName','Arial','interpreter','latex');
        %}
        %pubfig(h,12,2,0.005);

        % axis labels
        ylabel('% Vertical reports')
        xlabel('Gabor orientation')

        % legend
        legend({'Baseline','H stim (0\circ)','V-stim (90\circ)'})

        upFontSize(14,0.01)
        % ax ticks
        
        xlims=[-1 1].*roundup(abs(unique(xlim)),10);
        addSkippedTicks(xlims(1),xlims(2),20,'x')
        addSkippedTicks(0,100,12.5,'y')
    end
end
end