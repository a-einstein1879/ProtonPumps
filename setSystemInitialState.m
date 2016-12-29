function sS0 = setSystemInitialState(Voltage, SDChain, sP)
    %% System state
    field1 = 'systemStates';
    % A1-site state
    % nA1 = 1, occupied A1-site
    f1 = 'A1Site'; v1 = 1;
    % B1-site state
    % nB1 = 0, empty B1-site
    f2 = 'B1Site'; v2 = 0.5;
    % S-site state
    % nS = 0, empty S-site
    f3 = 'Source'; v3 = 0;
    % L2-site state
    % nL2 = 0, empty L2-site
    f4 = 'L2Site'; v4 = 0;
    % First quinone state
    % First quinone initial state (from 1 to 16) (Q1-basis on Feb/2/2011/A-6)
    MQ1 = 1;
    f5 = 'Quinone1'; v5 = zeros(sP.numberOfLevels.Q1, 1);
    v5(MQ1) = 1;
    % LH-system state
    % LH system at t=0 (from 1 to 4)(LH-basis on Feb/2/2011/A-5)
    MLH = 3;
    f6 = 'LHSystem'; v6 = zeros(sP.numberOfLevels.LH, 1);
    v6(MLH) = 1;
    value1 = struct(f1, v1, f2, v2, f3, v3, f4, v4, f5, v5, f6, v6);
    
    %% First quinone position
    % initial position of first Quinone on the N-side on membrane
    field2 = 'quinone1Position'; value2 = -sP.q1MM.x0;
    
    %% Chemical potentials
    % Source-L2 (electrons) 
    muSDPlus = 1620;
    field3 = 'chemicalPotentials';
    f1 = 'S';  v1 = 0.5 * (muSDPlus + SDChain);
    f2 = 'L2'; v2 = 0.5 * (muSDPlus - SDChain);
    f3 = 'N';  v3 = -sP.dM / 2;
    f4 = 'P';  v4 = sP.dM / 2;
    value3 = struct(f1, v1, f2, v2, f3, v3, f4, v4);
    
    %% Membrane potentials
    field4 = 'membranePotentials';
    f1 = 'dV'; v1 = 20;
    f2 = 'VN'; v2 = 0.5 * (Voltage - v1);
    f3 = 'VP'; v3 = 0.5 * (Voltage + v1);
    value4 = struct(f1, v1, f2, v2, f3, v3);
    
    %% Energy shifts of A,B,L,H energies due to NP field gradient
    field5 = 'energyShifts';
    % eA1 = eA1 - VN
    f1 = 'eA1'; v1 = sP.energyLevels.eA1 - value4.VN;
    % eD1 = eB1 + VP
    f2 = 'eB1'; v2 = sP.energyLevels.eB1 + value4.VP;
    % eH1 = eH1 - VN
    f3 = 'eH1'; v3 = sP.energyLevels.eH1 - value4.VN;
    % eL1 = eL1 + VP
    f4 = 'eL1'; v4 = sP.energyLevels.eL1 + value4.VP;
    % EQ1 = 1.5 * Uep - 0.5 * Up + 0.5 * (VP - VN) + 
    % + 0.5 * (muN + muP) - 25
    f5 = 'EQ1'; v5 = 1.5 * sP.Uep - 0.5 * sP.columbQuinone1.Up + ...
        0.5 * (value4.VP - value4.VN) + 0.5 * (value3.N + value3.P) - 25;
    value5 = struct(f1, v1, f2, v2, f3, v3, f4, v4, f5, v5);
    
    %% Gammas
    field6 = 'gammas';
    % ??????????????
    % Fermi distribution
    % fSeA1 = 1./( exp((eA-muS)/TT) + 1 )
    f1 = 'fSeA1';  v1 = 1 ./ (exp((value5.eA1 - value3.S) / sP.TT) + 1);
    % fB1eL2 = 1./( exp((eD-muD)/TT) + 1 )
    f2 = 'fB1eL2'; v2 = 1 ./ (exp((value5.eB1 - value3.L2) / sP.TT) + 1);

    % normed gammas
    GamSA1 = zeros(2,2); GamB1L2 = zeros(2,2);
    GamSA1(1,2)  = sP.gammas.gamS  .* (1 - v1);
    GamSA1(2,1)  = sP.gammas.gamS  .* v1;
    GamB1L2(1,2) = sP.gammas.gamL2 .* (1 - v2);
    GamB1L2(2,1) = sP.gammas.gamL2 .* v2;
    f3 = 'GamSA1';  v3 = GamSA1;
    f4 = 'GamB1L2'; v4 = GamB1L2;
    
    value6 = struct(f1, v1, f2, v2, f3, v3, f4, v4);
    
    %% Hamiltonians
    field7 = 'hamiltonians';
    % Feb/1/2011/A-2; no dependence on first Quinone position xQ
    % HQ0 = eQ * (n1 + n2) + EQ * (N1 + N2) + ue * n1 * n2 + ...
    % + Up * N1 * N2 - u11 * n1 * N1 - u12 * n1 * N2 - ...
    % - u21 * n2 * N1 - u22 * n2 * N2
    f1 = 'HQ10';
    n1 = sP.populationOperators.Q1n1;
    n2 = sP.populationOperators.Q1n2;
    N1 = sP.populationOperators.Q1N1;
    N2 = sP.populationOperators.Q1N2;
    nL1 = sP.populationOperators.nL1;
    nH1 = sP.populationOperators.nH1;
    
    eQ  = sP.energyLevels.eQ1;
    ue  = sP.columbQuinone1.ue;
    Up  = sP.columbQuinone1.Up;
    u11 = sP.columbQuinone1.u11;
    u12 = sP.columbQuinone1.u12;
    u21 = sP.columbQuinone1.u21;
    u22 = sP.columbQuinone1.u22;
    uLH = sP.columbQuinone1.uLH;
    v1 = ( n1 + n2 ) * eQ + ( N1 + N2 ) * value5.EQ1 + ...
           n1 * n2 * ue  + N1 * N2 * Up - ...
           n1 * N1 * u11 - n1 * N2 * u12 - ...
           n2 * N1 * u21 - n2 * N2 * u22;
    % Hamiltonian HLH
    % HLH0 = eL1*nL1 + eH1*nH1 + uLH*nL1*nH1
    f2 = 'HLH0';
    v2 = nL1 * value5.eL1 + nH1 * value5.eH1 + nL1 * nH1 * uLH;
    % energy spectrum of LH system
    f3 = 'EnLH'; v3 = diag(v2);
    EnLHd = repmat(v3, 1, sP.numberOfLevels.LH);
    % frequencies of LH system
    f4 = 'omegaLH'; v4 = EnLHd - EnLHd';
    value7 = struct(f1, v1, f2, v2, f3, v3, f4, v4);
    
    %% General state
    % First quinone charge squared
    % qqQ1 = (n1 + n2 - N1 - N2) * (n1 + n2 - N1 - N2)
    field8 = 'qqQ1'; value8 = (n1 + n2 - N1 - N2) * (n1 + n2 - N1 - N2);

    
    %% sS for system state
    sS0 = struct(field1, value1, field2, value2, field3, value3, ...
        field4, value4, field5, value5, field6, value6, field7, value7, ...
        field8, value8);
end