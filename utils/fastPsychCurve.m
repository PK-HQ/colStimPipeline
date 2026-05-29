dateStr='20260410';
run='2';
runStr=['\run' run '\'];
load(['Y:\Pepper\Pepper' dateStr runStr '\M32D' dateStr 'R' run 'TS.mat'])
optoFlag=TS.Header.Conditions.TypeCond==3;
visualOrt=(TS.Header.Conditions.GaborOrt(optoFlag) == 90) * 1 + (TS.Header.Conditions.GaborOrt(optoFlag) == 0) * -1;
x=TS.Header.Conditions.StimCon(optoFlag).*visualOrt; %contrasts
y= TS.Header.Outcomes.CountCondSuccess(optoFlag)'./TS.Header.Outcomes.CountCondTotalValid(optoFlag)' *100; % Percent correct by contrast
zeroConds=find(x==0);
if numel(zeroConds)>0
    y(zeroConds(1))=mean(y(zeroConds)); y(zeroConds(2))=[];
    x(zeroConds(1))=mean(x(zeroConds)); x(zeroConds(2))=[];
end

figure
% Split curve
subplot(1,2,1)
scatter(x(x>=0),y(x>=0),200,'ko', 'LineWidth',2, 'markerfacecolor','r'); hold on
scatter(x(x<=0),100-y(x<=0),200,'ko', 'LineWidth',2, 'markerfacecolor','b'); hold on
xlim([-100 100])
ylim([0 100])
xline(0,'--','linewidth',2,'color',[.65 .65 .65],'HandleVisibility','off'); hold on

yline(50,'--','linewidth',2,'color',[.65 .65 .65],'HandleVisibility','off'); hold on
axis square
addSkippedTicks(-100, 100, 12.5,'x')
addSkippedTicks(0, 100, 12.5,'y')
upFontSize
title('Raw', 'FontWeight','Normal')
ylabel('% vertical')
xlabel('Signed contrast (%)')
legend('H','V')

% Merged curve
subplot(1,2,2)
xMerged=mean([x(x>=0); -x(x<=0)]);
yMerged=mean([y(x>=0); y(x<=0)]);
scatter(xMerged,yMerged,200,'ko', 'LineWidth',2, 'markerfacecolor',' magenta'); hold on
xlim([0 100])
ylim([0 100])
yline(50,'--','linewidth',2,'color',[.65 .65 .65],'HandleVisibility','off'); hold on
axis square
addSkippedTicks(0, 100, 12.5,'x')
addSkippedTicks(0, 100, 12.5,'y')
title('Mean', 'FontWeight','Normal')
ylabel('% correct')
xlabel('Contrast (%)')
upFontSize
h = suplabel('Pepper 20260403', 't', [0.1 0.1 0.8 0.78]);
upFontSize