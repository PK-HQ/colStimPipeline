%
% Visual + Opto-stim Psychometric Functions
%   random contrast
%
clearvars;
figure
testParams=[0.1069    0.1363    0.0    0.0   60.7054    3.548   99.9961    0.388    20];
testFlag=~isempty(testParams);

switch testFlag
    case 0
        fprintf('Default\n')
        % excitatory parameters
        a = 0.066; % alpha, visual excitatory cross-talk parameter
        b = 0.024; % beta, opto-stim excitatory cross-talk parameter
        %
        % normalization parameters
        l = 0.0; % lambda, visual normalization cross-talk parameter
        w = 0.0; % omega, opto-stim normalization cross-talk parameter
        g0 = 30; % normalization constant
        %
        % other parameters
        n = 2;   % spiking exponent
        rmx = 3; % max response
        e = 0.5; % eta, relative weight of optostim on normalization vs. excitation
        o = 40;     % opto-stim effective contrast
    case 1
        fprintf('Test\n')
        % excitatory parameters
        a = testParams(1); % alpha, visual excitatory cross-talk parameter
        b = testParams(2); % beta, opto-stim excitatory cross-talk parameter
        %
        % normalization parameters
        l = testParams(3); % lambda, visual normalization cross-talk parameter
        w = testParams(4); % omega, opto-stim normalization cross-talk parameter
        g0 = testParams(5); % normalization constant
        %
        % other parameters
        n = testParams(6);   % spiking exponent
        rmx = testParams(7); % max response
        e = testParams(8); % eta, relative weight of optostim on normalization vs. excitation
        o = testParams(9);     % opto-stim effective contrast
end
%
% input: (0,0)=(Vs,Vo), (1,1)=(Hs,Ho), (0,1)=(Vs,Ho), (1,0)=(Hs,Vo)
inp = zeros(1,2);

%% stimulation condition
ntrl = 10000; % trials per level of stimulus contrast
cmax = 100; % max stimulus contrast
cmin = 0;   % min stimulus contrast
cstp = 5;   % contrast step size
nlev = (cmax-cmin)/cstp + 1;

c = (cmin:cstp:cmax);
% numerator mean storage
uhhv = zeros(1,nlev);
uvhv = zeros(1,nlev);
uhhh = zeros(1,nlev);
uvhh = zeros(1,nlev);
% denominator mean storage
uhvv = zeros(1,nlev);
uvvv = zeros(1,nlev);
uhvh = zeros(1,nlev);
uvvh = zeros(1,nlev);
for i = 1:nlev
% numerator means  
  ehhv = c(i) + b*(1-e)*o;
  ghhv = c(i) + w*e*o;
  uhhv(i) = rmx * (ehhv^n)/(ghhv^n + g0^n);
  %
  evhv = a*c(i) + (1-e)*o;
  gvhv = l*c(i) + e*o;
  uvhv(i) = rmx * (evhv^n)/(gvhv^n + g0^n);
  %
  ehhh = c(i) + (1-e)*o;
  ghhh = c(i) + e*o;
  uhhh(i) = rmx * (ehhh^n)/(ghhh^n + g0^n);
  %
  evhh = a*c(i) + b*(1-e)*o;
  gvhh = l*c(i) + w*e*o;
  uvhh(i) = rmx * (evhh^n)/(gvhh^n + g0^n);
% denominator means
  ehvv = a*c(i) + b*(1-e)*o;
  ghvv = l*c(i) + w*e*o;
  uhvv(i) = rmx * (ehvv^n)/(ghvv^n + g0^n);
  %
  evvv = c(i) + (1-e)*o;
  gvvv = c(i) + e*o;
  uvvv(i) = rmx * (evvv^n)/(gvvv^n + g0^n);
  %
  ehvh = a*c(i) + (1-e)*o;
  ghvh = l*c(i) + e*o;
  uhvh(i) = rmx * (ehvh^n)/(ghvh^n + g0^n);
  %
  evvh = c(i) + b*(1-e)*o;
  gvvh = c(i) + w*e*o;
  uvvh(i) = rmx * (evvh^n)/(gvvh^n + g0^n);
end
% run trials
trls = zeros(ntrl,4); % inp(1),inp(2),ic,hvrsp
for j = 1:ntrl
    inp = randi([0 1],1,2,"logical");
    ic = randi([1 nlev]);   % index to contrast
    trls(j,1:2) = inp;
    trls(j,3) = ic;
    if inp(1) == 0 && inp(2) == 0     % vert stm vert optstm
      rh = randn() + uhvv(ic);
      rv = randn() + uvvv(ic);
    elseif inp(1) == 0 && inp(2) == 1 % vert stm horz optstm
      rh = randn() + uhvh(ic);
      rv = randn() + uvvh(ic);
    elseif inp(1) == 1 && inp(2) == 0 % horz stm vert optstm
      rh = randn() + uhhv(ic);
      rv = randn() + uvhv(ic);      
    else                              % horz stm horz optstm
      rh = randn() + uhhh(ic);
      rv = randn() + uvhh(ic);  
    end
    num = sum(exp( -0.5*( (rh-uhhv).^2 + (rv-uvhv).^2 ) )...
        + exp( -0.5*( (rh-uhhh).^2 + (rv-uvhh).^2 ) ));
    den = sum(exp( -0.5*( (rh-uhvv).^2 + (rv-uvvv).^2 ) )...
        + exp( -0.5*( (rh-uhvh).^2 + (rv-uvvh).^2 ) ));
    lr = num/den;
    if lr > 1.0     % horz behav response
      trls(j,4) = 1;    
    elseif lr == 1  % random behav response
      if rand() > 0.5
         trls(j,4) = 1;
      end
    end             % else leave vert behav response
end
% process trials and save results
% trls = inp(1),inp(2),ic,hvrsp
% output = c,o,nhh,nhv,nvh,nvv,chh,chv,cvh,cvv
output = zeros(nlev,10);
for j = 1:ntrl
  ic = trls(j,3);
  output(ic,1) = c(ic);
  output(ic,2) = o;
  if trls(j,1) == 1 && trls(j,2) == 1
    output(ic,3) = output(ic,3) + 1;    % hh trial count
    if trls(j,1) == trls(j,4)
      output(ic,7) = output(ic,7) + 1;  % hh correct response count
    end
  elseif trls(j,1) == 1 && trls(j,2) == 0
    output(ic,4) = output(ic,4) + 1;    % hv trial count
    if trls(j,1) == trls(j,4)
      output(ic,8) = output(ic,8) + 1;  % hv correct response count
    end
  elseif trls(j,1) == 0 && trls(j,2) == 1
    output(ic,5) = output(ic,5) + 1;    % vh trial count
    if trls(j,1) == trls(j,4)
      output(ic,9) = output(ic,9) + 1;  % vh correct response count
    end
  elseif trls(j,1) == 0 && trls(j,2) == 0
    output(ic,6) = output(ic,6) + 1;    % vv trial count
    if trls(j,1) == trls(j,4)
      output(ic,10) = output(ic,10) + 1;% vv correct response count
    end
  end
end
pc = (output(:,7) + output(:,10))...
    ./(output(:,3) + output(:,6));
pic = (output(:,8) + output(:,9))...
    ./(output(:,4) + output(:,5));
%
plot(output(:,1),pc*100,"r-o",'LineWidth',2); hold on;
plot(output(:,1),pic*100,"b-o",'LineWidth',2);
% plot(output(:,1),output(:,7)./output(:,3),"b-o"); hold on;
% plot(output(:,1),output(:,8)./output(:,4),"r-o");
% plot(output(:,1),output(:,9)./output(:,5),"r-o");
% plot(output(:,1),output(:,10)./output(:,6),"b-o");
%
%% control condition
o = 0;
% numerator mean storage
uhhv = zeros(1,nlev);
uvhv = zeros(1,nlev);
uhhh = zeros(1,nlev);
uvhh = zeros(1,nlev);
% denominator mean storage
uhvv = zeros(1,nlev);
uvvv = zeros(1,nlev);
uhvh = zeros(1,nlev);
uvvh = zeros(1,nlev);
for i = 1:nlev
% numerator means  
  ehhv = c(i) + b*o;
  ghhv = c(i) + w*o;
  uhhv(i) = rmx * (ehhv^n)/(ghhv^n + g0^n);
  %
  evhv = a*c(i) + o;
  gvhv = l*c(i) + o;
  uvhv(i) = rmx * (evhv^n)/(gvhv^n + g0^n);
  %
  ehhh = c(i) + o;
  ghhh = c(i) + o;
  uhhh(i) = rmx * (ehhh^n)/(ghhh^n + g0^n);
  %
  evhh = a*c(i) + b*o;
  gvhh = l*c(i) + w*o;
  uvhh(i) = rmx * (evhh^n)/(gvhh^n + g0^n);
% denominator means
  ehvv = a*c(i) + b*o;
  ghvv = l*c(i) + w*o;
  uhvv(i) = rmx * (ehvv^n)/(ghvv^n + g0^n);
  %
  evvv = c(i) + o;
  gvvv = c(i) + o;
  uvvv(i) = rmx * (evvv^n)/(gvvv^n + g0^n);
  %
  ehvh = a*c(i) + o;
  ghvh = l*c(i) + o;
  uhvh(i) = rmx * (ehvh^n)/(ghvh^n + g0^n);
  %
  evvh = c(i) + b*o;
  gvvh = c(i) + w*o;
  uvvh(i) = rmx * (evvh^n)/(gvvh^n + g0^n);
end
% run trials
trls = zeros(ntrl,4);
for j = 1:ntrl
    inp = randi([0 1],1,2,"logical");
    ic = randi([1 nlev]);
    trls(j,1:2) = inp;
    trls(j,3) = ic;
    if inp(1) == 0 && inp(2) == 0
      rh = randn() + uhvv(ic);
      rv = randn() + uvvv(ic);
    elseif inp(1) == 0 && inp(2) == 1
      rh = randn() + uhvh(ic);
      rv = randn() + uvvh(ic);
    elseif inp(1) == 1 && inp(2) == 0
      rh = randn() + uhhv(ic);
      rv = randn() + uvhv(ic);      
    else
      rh = randn() + uhhh(ic);
      rv = randn() + uvhh(ic);  
    end
    num = sum(exp( -0.5*( (rh-uhhv).^2 + (rv-uvhv).^2 ) )...
        + exp( -0.5*( (rh-uhhh).^2 + (rv-uvhh).^2 ) ));
    den = sum(exp( -0.5*( (rh-uhvv).^2 + (rv-uvvv).^2 ) )...
        + exp( -0.5*( (rh-uhvh).^2 + (rv-uvvh).^2 ) ));
    lr = num/den;
    if lr > 1.0
      trls(j,4) = 1;
    elseif lr == 1
      if rand() > 0.5
         trls(j,4) = 1;
      end
    end
end
% process trials and save results
% output = c,o,nhh,nhv,nvh,nvv,chh,chv,cvh,cvv
output = zeros(nlev,10);
for j = 1:ntrl
  ic = trls(j,3);
  output(ic,1) = c(ic);
  output(ic,2) = o;
  if trls(j,1) == 1 && trls(j,2) == 1
    output(ic,3) = output(ic,3) + 1;
    if trls(j,1) == trls(j,4)
      output(ic,7) = output(ic,7) + 1;
    end
  elseif trls(j,1) == 1 && trls(j,2) == 0
    output(ic,4) = output(ic,4) + 1;
    if trls(j,1) == trls(j,4)
      output(ic,8) = output(ic,8) + 1;
    end
  elseif trls(j,1) == 0 && trls(j,2) == 1
    output(ic,5) = output(ic,5) + 1;
    if trls(j,1) == trls(j,4)
      output(ic,9) = output(ic,9) + 1;
    end
  elseif trls(j,1) == 0 && trls(j,2) == 0
    output(ic,6) = output(ic,6) + 1;
    if trls(j,1) == trls(j,4)
      output(ic,10) = output(ic,10) + 1;
    end
  end
end
%
pcntrl = (output(:,7)+output(:,8)+output(:,9)+output(:,10))...
    ./(output(:,3)+output(:,4)+output(:,5)+output(:,6));
plot(output(:,1),pcntrl*100,"k-o",'LineWidth',2);
ylim([0,100]); axis square
legend('congruent','inconguent','control');
na = num2str(a); nb = num2str(b); nl = num2str(l); nw = num2str(w);
ne = num2str(e);
ttl = append('(a,b,l,w,e) (',na,',',nb,',',nl,',',nw,',',ne,')');
title(ttl);
% num2str(a);
% plot(output(:,1),output(:,7)./output(:,3),"k-o");
% plot(output(:,1),output(:,8)./output(:,4),"k-o");
% plot(output(:,1),output(:,9)./output(:,5),"k-o");
% plot(output(:,1),output(:,10)./output(:,6),"k-o");
% plot(cp,pc,"^"); hold on;
% plot(cp,pic,"v");
% ylim([0.4,1.0]);
% plot(cp,pc,'--k');

% Cosmetic
xlabel('Gabor contrast')
ylabel('Correct')
upFontSize(24,.005)