function yPredicted = weibullMdl(x, params)
% yPred = weibullMdl(x, p)
% ---------------------------------------------
%  Three +ve-side cumulative-Weibull psychometric curves:
%     • Baseline     (black)
%     • Opto-congruent  (red, shifted up at x = 0)
%     • Opto-incongruent (blue, shifted down at x = 0)
%  Parameter vector p (length = 12) follows your earlier convention:
%     1  A_bl      (lapse; 1 – upper asymptote of baseline)
%     2  B_bl      (fixed at 0.5 in this model – keep as placeholder)
%     3  α_bl      (C50 of baseline)
%     4  β_bl      (slope of baseline)
%     5  ΔA_con
%     6  ΔB_con            <-- used symmetrically ±
%     7  Δα_con
%     8  Δβ_con
%     9  ΔA_inc
%    10  ΔB_inc   (ignored; symmetry ties B_inc = 0.5 – ΔB_con)
%    11  Δα_inc
%    12  Δβ_inc
%
%  Returns a struct with fields
%     yPred.pc      – congruent (%)
%     yPred.pic    – incongruent (%)
%     yPred.pcntrl – baseline (%)
%
%  All outputs are in 0-to-100 % units.

% ---------- baseline parameters ----------
A_bl   = params(1);
B_bl   = 0.5;                 % chance, fixed
alpha_bl = params(3);
beta_bl  = params(4);

% ---------- congruent (Δ from baseline) ----------
A_con   = A_bl   + params(5);
deltaB  = params(6);               % single ΔB controls both sides
B_con   = 0.5 + deltaB;
alpha_con = alpha_bl + params(7);
beta_con  = beta_bl  + params(8);

% ---------- incongruent (Δ from baseline, symmetric B) ----------
A_inc   = A_bl   + params(9);
B_inc   = 0.5 - deltaB;       % symmetric shift
alpha_inc = alpha_bl + params(11);
beta_inc  = beta_bl  + params(12);

% Weibull helper: Y = B + (1 – exp(-(x/α)^β)) * ( (1–A) – B )
weib = @(x, A, B, a, b) ...
       100 * ( B + (1 - exp(-(x ./ a).^b)) .* ( (1 - A) - B ) );

% ---------- generate predictions ----------
yPredicted.pc      = weib(x, A_con,  B_con,  alpha_con,  beta_con);
yPredicted.pic    = weib(x, A_inc,  B_inc,  alpha_inc,  beta_inc);
yPredicted.pcntrl = weib(x, A_bl,  B_bl,   alpha_bl,   beta_bl);

end
