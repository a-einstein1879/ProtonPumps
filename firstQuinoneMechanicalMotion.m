function deltaX = firstQuinoneMechanicalMotion(sS, sP)
    % First quinone mechanical motion with potentials and stochastic force
    % average value of quinone charge squared
    qqQ = diag(sS.qqQ1)' * sS.systemStates.Quinone1;
    % ???????????????? these parameters should be renamed
    xQ = sS.quinone1Position;
    CVw = sP.q1MM.CVw;
    ExW = sP.q1MM.ExW;
    CVch = sP.q1MM.CVch;
    exCh = sP.q1MM.exCh;
    ExQk = exp(xQ / sP.q1MM.membraneWallSlope);
    ExQCk = exp(xQ / sP.q1MM.penaltyPotentialSlope);
    factorDif = sP.q1MM.factorDif;
    
    deltaX =  - CVw*( (ExQk/ExW)/( (ExQk/ExW)+1)^2 - (ExQk*ExW)/( (ExQk*ExW)+1)^2 ) + ... % Walls
        qqQ*CVch*( (ExQCk/exCh)/( (ExQCk/exCh)+1)^2 - (ExQCk*exCh)/( (ExQCk*exCh)+1)^2 ) ... % barrier for charged quinone
        + factorDif*randn; % mechanical stochastic force
end