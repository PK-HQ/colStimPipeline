%{
monkeyName='Chip';
monkeyID='28';
sessionDate='20230113';
runNo='0';
load(['Y:\' monkeyName '\' monkeyName sessionDate '\run' runNo '\M' monkeyID 'D' sessionDate 'R' runNo 'TS.mat'])
%}

if isfield(TS.Header.Conditions, 'ProjImg')
    disp(TS.Header.Comments)
    disp(unique(TS.Header.Conditions.ProjImg'))
else
    disp('Not biasing block')
    disp(TS.Header.Comments)
end