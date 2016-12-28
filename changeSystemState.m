function sS = changeSystemState(sS, sP, gammaA, gammaB, gammaLH, gammaQ)
    NLevQ = sP.numberOfLevels.Q;
    NLevLH = sP.numberOfLevels.LH;
    
    % Change of state of A and D sites
    nA = sS.systemStates.ASite;
    dAState = gammaA(2, 1) .* (1 - nA) - gammaA(1, 2) .* nA;
    nD = sS.systemStates.DSite;
    dDState = gammaB(2, 1) .* (1 - nD) - gammaB(1, 2) .* nD;

    % Change of LH-system state
    gamMLH = sum(gammaLH);
    ReLH = gammaLH * sS.systemStates.LHSystem;
    dLHState = zeros(NLevLH, 1);
    for j = 1:NLevLH
        dLHState(j)= -(gamMLH(j)) * sS.systemStates.LHSystem(j) + ReLH(j);
    end

    % Change of Quinone state
    GamQM = sum(gammaQ);
    ReLQ = gammaQ * sS.systemStates.Quinone;
    dQState = zeros(NLevQ, 1);
    for j = 1:NLevQ
        dQState(j) = - GamQM(j)*sS.systemStates.Quinone(j) + ReLQ(j); 
    end

    deltaX = quinoneMechanicalMotion(sS, sP);

    sS.systemStates.ASite    = sS.systemStates.ASite    + dAState  * sP.dt;
    sS.systemStates.DSite    = sS.systemStates.DSite    + dDState  * sP.dt;
    sS.systemStates.Quinone  = sS.systemStates.Quinone  + dQState  * sP.dt;
    sS.systemStates.LHSystem = sS.systemStates.LHSystem + dLHState * sP.dt;
    sS.quinonePosition       = sS.quinonePosition       + deltaX   * sP.dt;
end