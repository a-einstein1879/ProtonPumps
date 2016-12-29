function [gammaA, gammaB] = calculateABGammas(sS, sP, OmegaQ)
    % Q -> A, Q -> D, S -> A, D -> B
    % Contribution of Quinone to RC relaxation rate
    xQ = sS.quinonePosition;
    Qs = sS.systemStates.Quinone;
    x0 = sP.qMM.x0;
    NLevRC = sP.numberOfLevels.RC;
    
    a1 = sP.CAOperators.a1;
    a2 = sP.CAOperators.a2;
    aa12 = abs(a1).^2 + abs(a2).^2;
    
    eA = sS.energyShifts.eA;
    eD = sS.energyShifts.eD;
    LamAQ = sP.lambdas.LamAQ;
    LamDQ = sP.lambdas.LamDQ;
    
    % Sums in formula Feb/3/2011/A-11 for "0 -> 1" and "1 -> 0"
    % transitions to A and D sites
    PhA = zeros(NLevRC, NLevRC);
    PhD = zeros(NLevRC, NLevRC);
    PhA(1,2) = sum((aa12' .* exp(-(OmegaQ - eA + LamAQ).^2./(4 * LamAQ * sP.TT))) * Qs);
    PhA(2,1) = sum((aa12  .* exp(-(OmegaQ + eA + LamAQ).^2./(4 * LamAQ * sP.TT))) * Qs);

    PhD(1,2) = sum((aa12' .* exp(-(OmegaQ - eD + LamDQ).^2./(4 * LamDQ * sP.TT))) * Qs);
    PhD(2,1) = sum((aa12  .* exp(-(OmegaQ + eD + LamDQ).^2./(4 * LamDQ * sP.TT))) * Qs);

    % dependence of tunneling rates Q-RC (electrons) on change of
    % position
    % ?????????????????????? check if this deltas are right
    % to both: A and H sites
    delta = sP.transferDistances.qToAB;
    WNel = exp(-2 * abs(xQ + x0) / delta);
    % to both: B and L sites
    WPel = exp(-2 * abs(xQ - x0) / delta);

    gamAQ = WNel .* sP.delam.DeLamAQ .* PhA;
    gamDQ = WPel .* sP.delam.DeLamDQ .* PhD;

    % Total evolution of RC (gammas for evolution of A and D sites)
    gammaA = (gamAQ + sS.gammas.GamSA) .* sP.meVtoTime;
    gammaB = (gamDQ + sS.gammas.GamDB) .* sP.meVtoTime;
end