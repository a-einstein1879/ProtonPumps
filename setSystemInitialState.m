function sS0 = setSystemInitialState(sP)
    %% System state
    field1 = 'systemStates';
    % A-site state
    % nA = 1, occupied A-site
    f1 = 'ASite'; v1 = 1;
    % B-site state
    % nB = 0, empty B-site
    f2 = 'BSite'; v2 = 0.5;
    % S-site state
    % nS = 0, empty S-site
    f3 = 'Source'; v3 = 0;
    % D-site state
    % nD = 0, empty D-site
    f4 = 'DSite'; v4 = 0;    
    % Quinone state
    % Quinone initial state (from 1 to 16) (Q-basis on Feb/2/2011/A-6)
    MQ = 1;
    f5 = 'Quinone'; v5 = zeros(sP.numberOfLevels.Q, 1);
    v5(MQ) = 1;
    % LH-system state
    % LH system at t=0 (from 1 to 4)(LH-basis on Feb/2/2011/A-5)
    MLH = 3;
    f6 = 'LHSystem'; v6 = zeros(sP.numberOfLevels.LH, 1);
    v6(MLH) = 1;
    value1 = struct(f1, v1, f2, v2, f3, v3, f4, v4, f5, v5, f6, v6);
    
    %% Quinone position
    % initial position of Quinone on the N-side on membrane
    field2 = 'quinonePosition'; value2 = -sP.qMM.x0;
    
    %% Chemical potentials
    % Source-Drain (electrons) 
    muSDPlus = 1620;
    field3 = 'chemicalPotentials';
    f1 = 'S';  v1 = 0.5 * (muSDPlus + sP.SDChain);
    f2 = 'D';  v2 = 0.5 * (muSDPlus - sP.SDChain);
    f3 = 'N';  v3 = -sP.dM / 2;
    f4 = 'P';  v4 = sP.dM / 2;
    value3 = struct(f1, v1, f2, v2, f3, v3, f4, v4);
    
    %% Membrane potentials
    field4 = 'membranePotentials';
    f1 = 'dV'; v1 = 20;
    f2 = 'VN'; v2 = 0.5 * (sP.VSet - v1);
    f3 = 'VP'; v3 = 0.5 * (sP.VSet + v1);
    value4 = struct(f1, v1, f2, v2, f3, v3);
    
    %% Energy shifts of A,B,L,H energies due to NP field gradient
    field5 = 'energyShifts';
    % eA = eA - VN
    f1 = 'eA'; v1 = sP.energyLevels.eA - value4.VN;
    % eB = eB + VP
    f2 = 'eB'; v2 = sP.energyLevels.eB + value4.VP;
    % eH = eH - VN
    f3 = 'eH'; v3 = sP.energyLevels.eH - value4.VN;
    % eD = eD + VP
    f4 = 'eL'; v4 = sP.energyLevels.eL + value4.VP;
    % EQ = 1.5 * Uep - 0.5 * Up + 0.5 * (VP - VN) + 
    % + 0.5 * (muN + muP) - 25
    f5 = 'EQ'; v5 = 1.5 * sP.potentials.Uep - 0.5 * sP.potentials.Up + ...
        0.5 * (value4.VP - value4.VN) + 0.5 * (value3.N + value3.P) - 25;
    value5 = struct(f1, v1, f2, v2, f3, v3, f4, v4, f5, v5);
    
    %% Gammas
    field6 = 'gammas';
    % ??????????????
    % Fermi distribution
    % fSeA = 1./( exp((eA-muS)/TT) + 1 )
    f1 = 'fSeA'; v1 = 1 ./ (exp((value5.eA - value3.S) / sP.TT) + 1);
    % fBeD = 1./( exp((eD-muD)/TT) + 1 )
    f2 = 'fBeD'; v2 = 1 ./ (exp((value5.eB - value3.D) / sP.TT) + 1);

    % normed gammas
    GamSA = zeros(2,2); GamBD = zeros(2,2);
    GamSA(1,2) = sP.gammas.gamS .* (1 - v1);
    GamSA(2,1) = sP.gammas.gamS .* v1;
    GamBD(1,2) = sP.gammas.gamD .* (1 - v2);
    GamBD(2,1) = sP.gammas.gamD .* v2;
    f3 = 'GamSA'; v3 = GamSA;
    f4 = 'GamBD'; v4 = GamBD;
    
    value6 = struct(f1, v1, f2, v2, f3, v3, f4, v4);
    
    %% Hamiltonians
    field7 = 'hamiltonians';
    % Feb/1/2011/A-2; no dependence on first Quinone position xQ
    % HQ0 = eQ * (n1 + n2) + EQ * (N1 + N2) + ue * n1 * n2 + ...
    % + Up * N1 * N2 - u11 * n1 * N1 - u12 * n1 * N2 - ...
    % - u21 * n2 * N1 - u22 * n2 * N2
    f1 = 'HQ0';
    n1 = sP.populationOperators.n1;
    n2 = sP.populationOperators.n2;
    N1 = sP.populationOperators.N1;
    N2 = sP.populationOperators.N2;
    nL = sP.populationOperators.nL;
    nH = sP.populationOperators.nH;
    v1 = ( n1 + n2 ) * sP.energyLevels.eQ + ...
         ( N1 + N2 ) * value5.EQ + ...
           n1 * n2   * sP.potentials.ue + ...
           N1 * N2   * sP.potentials.Up - ...
           n1 * N1   * sP.potentials.u11 - ...
           n1 * N2   * sP.potentials.u12 - ...
           n2 * N1   * sP.potentials.u21 - ...
           n2 * N2   * sP.potentials.u22;
    % Hamiltonian HLH
    % HLH0 = eL*nL + eH*nH + uLH*nL*nH
    f2 = 'HLH0';
    v2 = nL * value5.eL + nH * value5.eH + nL * nH * sP.potentials.uLH;
    % energy spectrum of LH system
    f3 = 'EnLH'; v3 = diag(v2);
    EnLHd = repmat(v3, 1, sP.numberOfLevels.LH);
    % frequencies of LH system
    f4 = 'omegaLH'; v4 = EnLHd - EnLHd';
    value7 = struct(f1, v1, f2, v2, f3, v3, f4, v4);
    
    %% General state
    % First quinone charge squared
    % qqQ = (n1 + n2 - N1 - N2) * (n1 + n2 - N1 - N2)
    field8 = 'qqQ'; value8 = (n1 + n2 - N1 - N2) * (n1 + n2 - N1 - N2);

    
    %% sS for system state
    sS0 = struct(field1, value1, field2, value2, field3, value3, ...
        field4, value4, field5, value5, field6, value6, field7, value7, ...
        field8, value8);
end