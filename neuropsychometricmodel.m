%% Naka-rushton function, i.e. Rmax*[contrast^n]/[contrast^n + C50^n]
respX=0:0.01:1;
rmax=max(respX);
c50=0.4;
nPowers=[.5 1 2 3 4];
figure('name','Neuro-Psychometric modeling')
for n=1:numel(nPowers)
    nPower=nPowers(n);
    parameters=[rmax, c50, nPower];
    probVertical=ComputeNakaRushton(parameters, respX);
    plot(respX, probVertical, '-','linewidth',3,...
        'DisplayName', [sprintf('R_{max}=%.0g, ',rmax),...
        sprintf('C_{50}=%.0g, ',c50),...
        sprintf('n=%.0g',nPower)]); hold on;
    if n==1
        title('Columnar response function (Naka-Rushton)','FontWeight','Normal');
        ylabel('Columnar response')
        xlabel('Visual- & opto-stimulus value')
    end
    axis square;upFontSize(24,.005)
end
[h,icons] = legend('location','northwest');
icons = findobj(icons, 'Type', 'Line');
nColors=plasma(numel(icons)/2 +1); %colororder(nColors);

%% Neurometric response function parameters
% --- Naka Rushton params ---
c50=[.5];
nPower=20; %2-4
parameters=[c50, nPower];

% --- Parameters: Visual-driven responses ---
% visual stimulus inputs (0 to 100%)
contrasts=[-100:1:100];
% column baseline activity (normalized a.u., 0 to 1) *estimate from dF/F
contrasthBL=0;
contrastvBL=0;
% visual stimulus weights (normalized a.u., 0 to 1) *estimate from dF/F
alphaHh=[.1]; %H = for horz columns, h = for horizontal stimuli
alphaHv=0;
alphaVh=0; %V = for vert columns, v = for horizontal stimuli
alphaVv=alphaHh;

% --- Parameters: Optostim-driven responses ---
% opto stimulus inputs (normalized a.u., 0 to 1) *estimate from dF/F
optoH=0;
optoV=0;
optoBL=[.1];
% opto stimulus weights
betaHh=[.1]; %H = for horz columns, h = for horizontal stimuli
betaHv=0;
betaVh=0; %V = for vert columns, v = for horizontal stimuli
betaVv=betaHh;

%% make struct
%{
% Define unique field names as a cell array
fieldNames = {'contrasts','contrastHBL', 'contrastVBL', 'alphaHh', 'alphaHv', 'alphaVh', 'alphaVv', 'optoH', 'optoV', 'optoBL', 'betaHh', 'betaHv', 'betaVh', 'betaVv'};

% Create a sample 10x13 array (you can replace this with your data)
arrayData = rand(10, numel(fieldNames));

% Initialize an empty struct with unique field names
myStruct = struct();

% Populate the struct fields with data from the array
for i = 1:numel(fieldNames)
    myStruct.(fieldNames{i}) = arrayData(:, i);
end

% Access the struct fields
disp(myStruct.contrastHBL);  % Access the first field
disp(myStruct.betaVv);       % Access the last field
%}

%% Predicted vertical reports
for nContrast=1:numel(contrasts)
    contrast=contrasts(nContrast);
    if contrast<0
        contrastH=contrast*-1;
        contrastV=0;
    elseif contrast>0
        contrastH=0;
        contrastV=contrast;
    else
        contrastH=0;
        contrastV=0;
    end

    % --- Responses: horizontal and vertical columns ---
    % Columnar responses (H and V)
    respH=contrasthBL + alphaHh*contrastH + alphaHv*contrastV + ...
      betaHh*(optoH + optoBL) + betaHv*(optoV + optoBL);
    respV=contrastvBL + alphaVh*contrastH + alphaVv*contrastV + ...
        betaVh*(optoH + optoBL) + betaVv*(optoV + optoBL);

    % feed through naka rushton
    RespHmean=ComputeNakaRushton(parameters, respH);
    RespVmean=ComputeNakaRushton(parameters, respV);
    RespHstd=sqrt(RespHmean);
    RespVstd=sqrt(RespVmean);

    deltaResp=RespVmean-RespHmean;
    scaling=sqrt(RespHstd^2 + RespVstd^2);
    respTrial=2*deltaResp/scaling;
    prctHorizontal(nContrast)=normcdf(respTrial);
end
figure('Name','Model')
scatter(contrasts, prctHorizontal,50,'ko','linewidth',2)
ylabel('Vertical reports (%)')
xlabel('Contrast')
upFontSize(24,0.005)
ylim([0 1])
xlim([-100 100])