function OmegaQ = calculateQuinoneFrequencies(sP, sS)
    % Change of energy levels of first quinone electrons and protons
    % that happen due to effects of surface potential (transmembrane voltage)
    % e - electron, E - proton. Q - first quinone
    xQ = sS.quinone1Position;
    x0 = sP.q1MM.x0;
    VN = sS.membranePotentials.VN;
    VP = sS.membranePotentials.VP;
    
    a1 = sP.CAOperators.Q1a1;
    a2 = sP.CAOperators.Q1a2;
    b1 = sP.CAOperators.Q1b1;
    b2 = sP.CAOperators.Q1b2;
    n1 = a1' * a1;
    n2 = a2' * a2;
    N1 = b1' * b1;
    N2 = b2' * b2;
    
    VseQ =  VN * (xQ - x0) / (2 * x0) + VP * (xQ + x0) / (2 * x0);
    
    if xQ < -x0
        VseQ = -VN;
    elseif xQ > x0
        VseQ = VP;
    end
        
    VsEQ = - VseQ;
        
    % First Quinone Hamiltonian
    % change of energy levels due to 
    % dV-xQ-dependent part HQx
    HQx = VseQ * (n1 + n2) + VsEQ * (N1 + N2);
    
    % penalty for charged shuttle inside the membrane
    % the slope of the walls is aproximated with exponent
    % ???????????????
    ExQk = exp(xQ / sP.q1MM.membraneWallSlope);
    ExQCk = exp(xQ / sP.q1MM.penaltyPotentialSlope);
    Ucharge = sP.q1MM.penaltyPotentialHeight * ...
       (1 / ((ExQCk / sP.q1MM.exCh) + 1) - 1 / ((ExQCk * sP.q1MM.exCh) + 1));
    
    HQcharge = Ucharge * sS.qqQ1;
    % Total Hamiltonian of first Quinone (without tunnelings)
    HQA = sS.hamiltonians.HQ10 + HQcharge + HQx;
        
    % Energy levels of basis states of first Quinone (with xQ-dependence)
    EnQ = diag(HQA);

    % Duplication of column EnQ (result: NLevQ stobzov of EnQ)
    EnQd = repmat(EnQ, 1, sP.numberOfLevels.Q1);
    % frequencies of first Quinone
    OmegaQ = EnQd - EnQd';
end