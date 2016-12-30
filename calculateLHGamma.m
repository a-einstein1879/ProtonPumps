function gammaLH = calculateLHGamma(sP, sS, OmegaQ)
    % Total gamma of LH system, is summ of
    % LH tunneling, LQ1 tunneling and HQ1 tunneling
    % Contribution of first Quinone to LH evolution
    NLevLH = sP.numberOfLevels.LH;
    
    a1 = sP.CAOperators.Q1a1;
    a2 = sP.CAOperators.Q1a2;
    aL = sP.CAOperators.aL1;
    aH = sP.CAOperators.aH1;
    aaL = abs(aL).^2;
    aaH = abs(aH).^2;
    aaLH = abs(aL' * aH).^2;
    aa12 = abs(a1).^2 + abs(a2).^2;
    
    TT = sP.TT;
    omegaLH = sS.hamiltonians.omegaLH;
    DeLamLQ = sP.delam.DeLamL1Q1;
    DeLamHQ = sP.delam.DeLamH1Q1;
    LamLQ = sP.lambdas.LamL1Q1;
    LamHQ = sP.lambdas.LamH1Q1;
    LamLH = sP.lambdas.LamL1H1;
    dLH = sP.tunnelingAmplitudes.dL1H1;
    
    Qs = sS.systemStates.Quinone1;

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
    xQ = sS.quinone1Position;
    x0 = sP.q1MM.x0;
    delta = sP.transferDistances.q1ToAB1;
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