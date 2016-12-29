function [gammaA, gammaB] = calculateABGammas(sS, sP, OmegaQ)
    % Q1 -> A1, Q1 -> B1, S -> A1, B1 -> L2
    % Contribution of first Quinone to RC relaxation rate
    xQ = sS.quinone1Position;
    Qs = sS.systemStates.Quinone1;
    x0 = sP.q1MM.x0;
    NLevRC = sP.numberOfLevels.RC;
    
    a1 = sP.CAOperators.Q1a1;
    a2 = sP.CAOperators.Q1a2;
    aa12 = abs(a1).^2 + abs(a2).^2;
    
    eA = sS.energyShifts.eA1;
    eB = sS.energyShifts.eB1;
    LamAQ = sP.lambdas.LamA1Q1;
    LamBQ = sP.lambdas.LamB1Q1;
    
    % Sums in formula Feb/3/2011/A-11 for "0 -> 1" and "1 -> 0"
    % transitions to A1 and L2 sites
    PhA = zeros(NLevRC, NLevRC);
    PhB = zeros(NLevRC, NLevRC);
    PhA(1,2) = sum((aa12' .* exp(-(OmegaQ - eA + LamAQ).^2./(4 * LamAQ * sP.TT))) * Qs);
    PhA(2,1) = sum((aa12  .* exp(-(OmegaQ + eA + LamAQ).^2./(4 * LamAQ * sP.TT))) * Qs);

    PhB(1,2) = sum((aa12' .* exp(-(OmegaQ - eB + LamBQ).^2./(4 * LamBQ * sP.TT))) * Qs);
    PhB(2,1) = sum((aa12  .* exp(-(OmegaQ + eB + LamBQ).^2./(4 * LamBQ * sP.TT))) * Qs);

    % dependence of tunneling rates Q1-RC (electrons) on change of
    % position
    % ?????????????????????? check if this deltas are right
    % to both: A1 and H1 sites
    delta = sP.transferDistances.q1ToAB1;
    WNel = exp(-2 * abs(xQ + x0) / delta);
    % to both: B1 and L1 sites
    WPel = exp(-2 * abs(xQ - x0) / delta);

    gamAQ = WNel .* sP.delam.DeLamA1Q1 .* PhA;
    gamBQ = WPel .* sP.delam.DeLamB1Q1 .* PhB;

    % Total evolution of RC (gammas for evolution of A1 and L2 sites)
    gammaA = (gamAQ + sS.gammas.GamSA1)  .* sP.meVtoTime;
    gammaB = (gamBQ + sS.gammas.GamB1L2) .* sP.meVtoTime;
end