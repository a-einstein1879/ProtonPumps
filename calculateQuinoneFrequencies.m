function OmegaQ = calculateQuinoneFrequencies(sP, sS)
    xQ = sS.quinonePosition;
    x0 = sP.qMM.x0;
    VN = sS.membranePotentials.VN;
    VP = sS.membranePotentials.VP;
    
    a1 = sP.CAOperators.a1;
    a2 = sP.CAOperators.a2;
    b1 = sP.CAOperators.b1;
    b2 = sP.CAOperators.b2;
    n1 = a1' * a1;
    n2 = a2' * a2;
    N1 = b1' * b1;
    N2 = b2' * b2;

    % Change of energy levels of quinone electrons and protons
    % that happen due to effects of surface potential (transmembrane voltage)
    % e - electron, E - proton. Q - quinone
    VseQ =  VN * (xQ - x0) / (2 * x0) + VP * (xQ + x0) / (2 * x0);
    
    if xQ < -x0
        VseQ = -VN;
    elseif xQ > x0
        VseQ = VP;
    end
        
    VsEQ = - VseQ;
        
    % Quinone Hamiltonian
    % change of energy levels due to 
    % dV-xQ-dependent part HQx
    HQx = VseQ * (n1 + n2) + VsEQ * (N1 + N2);
    
    % penalty for charged shuttle inside the membrane
    % the slope of the walls is aproximated with exponent
    % ???????????????
    ExQk = exp(xQ / sP.qMM.membraneWallSlope);
    ExQCk = exp(xQ / sP.qMM.penaltyPotentialSlope);
    Ucharge = sP.qMM.penaltyPotentialHeight * ...
       (1 / ((ExQCk / sP.qMM.exCh) + 1) - 1 / ((ExQCk * sP.qMM.exCh) + 1));
    
    HQcharge = Ucharge * sS.qqQ;
    % Total Hamiltonian of Quinone (without tunnelings)
    H = sS.hamiltonians.HQ0 + HQcharge + HQx;
        
    % Energy levels of basis states of Quinone (with xQ-dependence)
    EnQ = diag(H);

    % Duplication of column EnQ (result: NLevQ stobzov of EnQ)
    EnQd = repmat(EnQ, 1, sP.numberOfLevels.Q);
    % frequencies of Quinone
    OmegaQ = EnQd - EnQd';
end