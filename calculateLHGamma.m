function gammaLH = calculateLHGamma(sP, sS, OmegaQ)
    % Total gamma of LH system, is summ of
    % LH tunneling, LQ tunneling and HQ tunneling
    % Contribution of Quinone to LH evolution
    NLevLH = sP.numberOfLevels.LH;
    LamLQ = sP.lambdas.LamLQ;
    LamHQ = sP.lambdas.LamHQ;
    LamLH = sP.lambdas.LamLH;
    dLH = sP.tunnelingAmplitudes.dLH;
    
    a1 = sP.CAOperators.a1;
    a2 = sP.CAOperators.a2;
    aL = sP.CAOperators.aL;
    aH = sP.CAOperators.aH;
    aaL = abs(aL).^2;
    aaH = abs(aH).^2;
    aaLH = abs(aL' * aH).^2;
    aa12 = abs(a1).^2 + abs(a2).^2;
    
    TT = sP.TT;
    omegaLH = sS.hamiltonians.omegaLH;
    DeLamLQ = sP.delam.DeLamLQ;
    DeLamHQ = sP.delam.DeLamHQ;
    Qs = sS.systemStates.Quinone;

    PhL = zeros(NLevLH,NLevLH);
    PhH = zeros(NLevLH,NLevLH);

    for M = 1:NLevLH
        for N = 1:NLevLH
            % see Feb/3/2011/A-4
            omLun = sS.hamiltonians.omegaLH(M, N);
            ePart = exp(-(OmegaQ + omLun + LamLQ).^2./(4 * LamLQ * TT));
            PhL(M,N) = sum(((aaL(N,M).*aa12 + aaL(M,N).*aa12').*ePart) * Qs);

            ePart = exp(-(OmegaQ + omLun + LamHQ).^2./(4 * LamHQ * TT));
            PhH(M,N) = sum(((aaH(N,M).*aa12 + aaH(M,N).*aa12').*ePart) * Qs);
        end
    end

    % see Feb/3/2011/A-5, matrix NLevRC X NLevRC
    % gamLHQ = WPel.*(DeLamL1*PhL1 + DeLamL2*PhL2) + WNel.*(DeLamH1*PhH1 + DeLamH2*PhH2)
    xQ = sS.quinonePosition;
    x0 = sP.qMM.x0;
    delta = sP.transferDistances.qToAB;
    WPel = exp(-2 * abs(xQ - x0) / delta);
    WNel = exp(-2 * abs(xQ + x0) / delta);
    gamLHQ = WPel .* DeLamLQ * PhL + WNel .* DeLamHQ * PhH;

    % Contribution of LH tunneling to LH evolution
    % Feb/3/2011/A-3
    DeLamLHsym = abs(dLH).^2 * sqrt(pi / (LamLH * TT)) * (aaLH + aaLH');
    gamLHtun = DeLamLHsym .* exp(-(omegaLH + LamLH).^2./(4 * LamLH * TT));

    % Total evolution of LH system
    gammaLH = gamLHtun + gamLHQ;
    gammaLH = gammaLH .* sP.meVtoTime;
end