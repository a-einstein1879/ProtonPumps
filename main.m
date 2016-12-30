function [] = main()
    tic;
    sP = setSystemParameters();
    sS0 = setSystemInitialState(sP);
    [dynamicsOutput, cumulativeOutput] = runSimulation(sP, sS0);
    x = 0:sP.nOS;
    hold on
    plot(x, dynamicsOutput.nA);
    plot(x, dynamicsOutput.nB);
    plot(x, dynamicsOutput.nL);
    plot(x, dynamicsOutput.nH);
    hold off
    toc;
end

function [dynamicsOutput, cumulativeOutput] = runSimulation(sP, sS)
    % Simulation cycle
    % definition of output array
    f1 = 'quinonePosition'; v1 = zeros(1, sP.nOS);
    f2 = 'nA';              v2 = zeros(1, sP.nOS);
    f3 = 'nB';              v3 = zeros(1, sP.nOS);
    f4 = 'nL';              v4 = zeros(1, sP.nOS);
    f5 = 'nH';              v5 = zeros(1, sP.nOS);
    f6 = 'n1Q';             v6 = zeros(1, sP.nOS);
    f7 = 'n2Q';             v7 = zeros(1, sP.nOS);
    f8 = 'N1Q';             v8 = zeros(1, sP.nOS);
    f9 = 'N2Q';             v9 = zeros(1, sP.nOS);
    f10 = 'nD';             v10 = zeros(1, sP.nOS);
    for time = 2:sP.nOS
        OmegaQ = calculateQuinoneFrequencies(sP, sS);
        
        [gammaA, gammaB] = calculateABGammas(sS, sP, OmegaQ);
        
        gammaLH = calculateLHGamma(sP, sS, OmegaQ);
        
        gammaQ = calculateQuinoneGamma(sP, sS, OmegaQ);
        
%        OmegaQ
%        gammaA
%        gammaB
%        gammaLH
%        gammaQ
%        WNpr
%        WPpr
%        pause(10)
        sS = changeSystemState(sS, sP, gammaA, gammaB, gammaLH, gammaQ);
        
        %% Fill in the output array
        v1(time + 1) = sS.quinonePosition;
        v2(time + 1) = sS.systemStates.ASite;
        v3(time + 1) = sS.systemStates.BSite;
        v4(time + 1) = sum(sP.populationOperators.nL * sS.systemStates.LHSystem);
        v5(time + 1) = sum(sP.populationOperators.nH * sS.systemStates.LHSystem);
        v6(time + 1) = sum(sP.populationOperators.n1 * sS.systemStates.Quinone);
        v7(time + 1) = sum(sP.populationOperators.n2 * sS.systemStates.Quinone);
        v8(time + 1) = sum(sP.populationOperators.N1 * sS.systemStates.Quinone);
        v9(time + 1) = sum(sP.populationOperators.N2 * sS.systemStates.Quinone);
        v10(time + 1) = sS.systemStates.DSite;
    end
    
    fileID = fopen('testout.txt','w');
    fprintf(fileID, 'nA\tnB\tnL\tnH\tnD\r\n');
    for time = 2:sP.nOS
        fprintf(fileID, '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\r\n', v2(time), v3(time), v4(time), v5(time), v10(time));
    end
    fclose(fileID);
    
    dynamicsOutput = struct(f1, v1, f2, v2, f3, v3, f4, v4, ...
        f5, v5, f6, v6, f7, v7, f8, v8, f9, v9, f10, v10);
    
    cumulativeOutput = 0;
end