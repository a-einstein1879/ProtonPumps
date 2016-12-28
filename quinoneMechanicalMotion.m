function deltaX = quinoneMechanicalMotion(sS, sP)
    % Quinone mechanical motion with potentials and stochastic force
    % average value of quinone charge squared
    qqQ = diag(sS.qqQ)' * sS.systemStates.Quinone;
    % ???????????????? these parameters should be renamed
    xQ = sS.quinonePosition;
    CVw = sP.qMM.CVw;
    ExW = sP.qMM.ExW;
    CVch = sP.qMM.CVch;
    exCh = sP.qMM.exCh;
    ExQk = exp(xQ / sP.qMM.membraneWallSlope);
    ExQCk = exp(xQ / sP.qMM.penaltyPotentialSlope);
    factorDif = sP.qMM.factorDif;
    
    deltaX =  - CVw*( (ExQk/ExW)/( (ExQk/ExW)+1)^2 - (ExQk*ExW)/( (ExQk*ExW)+1)^2 ) + ... % Walls
        qqQ*CVch*( (ExQCk/exCh)/( (ExQCk/exCh)+1)^2 - (ExQCk*exCh)/( (ExQCk*exCh)+1)^2 ) ... % barrier for charged quinone
        + factorDif*randn; % mechanical stochastic force
end