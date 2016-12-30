function [gammaQ, WNpr, WPpr] = calculateQuinoneGamma(sP, sS, OmegaQ)
    %% Calculate total gamma of quinone
    % GamQA + GamQD + GamQLH + GamQNP
    % contribution of RC to Quinone evolution, Feb/3/2011/A-6
    % contribution of LH system to Quinone evolution
    NLevQ = sP.numberOfLevels.Q;

    a1 = sP.CAOperators.a1;
    a2 = sP.CAOperators.a2;
    b1 = sP.CAOperators.b1;
    b2 = sP.CAOperators.b2;
    aL = sP.CAOperators.aL;
    aH = sP.CAOperators.aH;
    aaL = abs(aL).^2;
    aaH = abs(aH).^2;
    aa12 = abs(a1).^2 + abs(a2).^2;
    bb12 = abs(b1).^2 + abs(b2).^2;
    
    LHs = sS.systemStates.LHSystem;
    DeLamAQ = sP.delam.DeLamAQ;
    DeLamBQ = sP.delam.DeLamBQ;
    DeLamLQ = sP.delam.DeLamLQ;
    DeLamHQ = sP.delam.DeLamHQ;
    TT = sP.TT;
    muN = sS.chemicalPotentials.N;
    muP = sS.chemicalPotentials.P;
    GamN = sP.gammas.gamN;
    GamP = sP.gammas.gamP;
    omegaLH = sS.hamiltonians.omegaLH;
    LamLQ = sP.lambdas.LamLQ;
    LamHQ = sP.lambdas.LamHQ;
    LamAQ = sP.lambdas.LamAQ;
    LamBQ = sP.lambdas.LamBQ;
    eA = sS.energyShifts.eA;
    eB = sS.energyShifts.eB;

    QhL = zeros(NLevQ,NLevQ);
    QhH = zeros(NLevQ,NLevQ);

    for M=1:NLevQ
        for N=1:NLevQ
            OmRunQ = OmegaQ(M,N); 
            % LH-to-Quinone
            ePart = exp(-(OmRunQ + omegaLH + LamLQ).^2./(4 * LamLQ * TT));
            QhL(M,N) = sum(((aaL.*aa12(N,M) + aaL'.*aa12(M,N)).*ePart) * LHs);
            ePart = exp(-(OmRunQ + omegaLH + LamHQ).^2./(4 * LamHQ * TT));
            QhH(M,N) = sum(((aaH.*aa12(N,M) + aaH'.*aa12(M,N)).*ePart) * LHs);
        end
    end
    
    xQ = sS.quinonePosition;
    x0 = sP.qMM.x0;
    delta = sP.transferDistances.qToAB;
    WPel = exp(-2 * abs(xQ - x0) / delta);
    WNel = exp(-2 * abs(xQ + x0) / delta);
    
    % contributions of A and B to Quinone evolution
    ePart1 = exp(-(OmegaQ + eA + LamAQ).^2./(4 * LamAQ * TT));
    ePart2 = exp(-(OmegaQ - eA + LamAQ).^2./(4 * LamAQ * TT));
    n1 = sS.systemStates.ASite;
    GamQA = WNel .* DeLamAQ .* (aa12 .* ePart1 .* (1 - n1) + aa12' .* ePart2 .* n1);

    ePart1 = exp(-(OmegaQ + eB + LamBQ).^2./(4*LamBQ*TT));
    ePart2 = exp(-(OmegaQ - eB + LamBQ).^2./(4*LamAQ*TT));
    n2 = sS.systemStates.BSite;
    GamQB = WPel .* DeLamBQ .* (aa12 .* ePart1 .* (1 - n2) + aa12' .* ePart2 .* n2);

    % LH-to-Q relaxation matrix, NLevQ X NLevQ
    % GamQLH = WPel.*(DeLamL1*QhL1 + DeLamL2*QhL2) + WNel.*(DeLamH1*QhH1 + DeLamH2*QhH2);
    GamQLH = WPel .* DeLamLQ * QhL + WNel .* DeLamHQ * QhH;

    % contribution of proton transitions betwee N/P sides of membrane and
    % Quinone to Quinone evolution
    OmQT = exp(OmegaQ ./ TT); % NLevQ X NLevQ matrix
    ExxN = exp(-muN / TT);
    ExxP = exp(-muP / TT); % numbers

    % Fermi distributions of protons in N and P reservoirs, matrices
    FN = 1./(OmQT.*ExxN + 1);
    FP = 1./(OmQT.*ExxP + 1);

    GamFN = GamN.*bb12.*FN';
    GamFP = GamP.*bb12.*FP';
    GamQN = GamN.*bb12 - GamFN + GamFN';
    GamQP = GamP.*bb12 - GamFP + GamFP';

    % dependence of tunneling rates Q-NP (protons) on change of
    % position
    % to both: N and P sites
    delta = sP.transferDistances.qToN;
    WNpr = (exp((x0 + xQ) / delta) + 1) .^ (-2); % Near N-side (source of protons)
    WPpr = (exp((x0 - xQ) / delta) + 1) .^ (-2); % Near P-side (drain of protons)

    % Quinone relaxation matrix due to coupling to leads
    GamQNP = WNpr .* GamQN + WPpr .* GamQP;
    % Total quinone relaxation matrix
    gammaQ = GamQA + GamQB + GamQLH + GamQNP;
    gammaQ = gammaQ .* sP.meVtoTime;
end