load('Y:\Pepper\Pepper20250403\M32D20250403R3TS.mat')
optoFlag=TS.Header.Conditions.TypeCond==3;
visualOrt=(TS.Header.Conditions.GaborOrt(optoFlag) == 90) * 1 + (TS.Header.Conditions.GaborOrt(optoFlag) == 0) * -1;
x=TS.Header.Conditions.StimCon(optoFlag).*visualOrt; %contrasts
y= TS.Header.Outcomes.CountCondSuccess(optoFlag)'./TS.Header.Outcomes.CountCondTotalValid(optoFlag)' *100; % Percent correct by contrast
zeroConds=find(x==0);
y(zeroConds(1))=mean(y(zeroConds)); y(zeroConds(2))=[];
x(zeroConds(1))=mean(x(zeroConds)); x(zeroConds(2))=[];

figure
scatter(x(x>=0),y(x>=0),200,'ko', 'LineWidth',2, 'markerfacecolor','r'); hold on
scatter(-x(x<=0),y(x<=0)-1,200,'ko', 'LineWidth',2, 'markerfacecolor','b'); hold on
xlim([0 50])
ylim([0 100])

yline(50,'--','linewidth',2,'color',[.65 .65 .65],'HandleVisibility','off'); hold on
axis square
addSkippedTicks(0, 50, 12.5,'x')
addSkippedTicks(0, 100, 12.5,'y')
upFontSize
title('Pepper 20250403R3', 'FontWeight','Normal')
ylabel('% Vertical')
xlabel('% Correct')