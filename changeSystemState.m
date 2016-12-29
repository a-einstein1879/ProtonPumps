function sS = changeSystemState(sS, sP, gammaA, gammaB, gammaLH, gammaQ)
    NLevQ = sP.numberOfLevels.Q;
    NLevLH = sP.numberOfLevels.LH;
    
    % Change of state of A and D sites
    nA = sS.systemStates.A1Site;
    dAState = gammaA(2, 1) .* (1 - nA) - gammaA(1, 2) .* nA;
    nB = sS.systemStates.B1Site;
    dBState = gammaB(2, 1) .* (1 - nB) - gammaB(1, 2) .* nB;

    % Change of LH-system state
    gamMLH = sum(gammaLH);
    ReLH = gammaLH * sS.systemStates.LHSystem;
    dLHState = zeros(NLevLH, 1);
    for j = 1:NLevLH
        dLHState(j)= -(gamMLH(j)) * sS.systemStates.LHSystem(j) + ReLH(j);
    end

    % Change of Quinone state
    GamQM = sum(gammaQ);
    ReLQ = gammaQ * sS.systemStates.Quinone1;
    dQState = zeros(NLevQ, 1);
    for j = 1:NLevQ
        dQState(j) = - GamQM(j)*sS.systemStates.Quinone1(j) + ReLQ(j); 
    end

    deltaX = quinoneMechanicalMotion(sS, sP);

    sS.systemStates.A1Site    = sS.systemStates.A1Site   + dAState  * sP.dt;
    sS.systemStates.B1Site    = sS.systemStates.B1Site   + dBState  * sP.dt;
    sS.systemStates.Quinone1  = sS.systemStates.Quinone1 + dQState  * sP.dt;
    sS.systemStates.LHSystem  = sS.systemStates.LHSystem + dLHState * sP.dt;
    sS.quinone1Position       = sS.quinone1Position      + deltaX   * sP.dt;
end