function [] = main()
    tic;
    sP = setSystemParameters();
    sS0 = setSystemInitialState(sP);
    [dynamicsOutput, cumulativeOutput] = runSimulation(sP, sS0);
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
    f10 = 'nL2';            v10 = zeros(1, sP.nOS);
    for time = 2:sP.nOS
        OmegaQ = calculateFirstQuinoneFrequencies(sP, sS);
        
        [gammaA, gammaB] = calculateFirstABGammas(sS, sP, OmegaQ);
        
        gammaLH = calculateLHSystemGamma(sP, sS, OmegaQ);
        
        [gammaQ, WNpr, WPpr] = calculateFirstQuinoneGamma(sP, sS, OmegaQ);
        
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
        v1(time + 1) = sS.quinone1Position;
        v2(time + 1) = sS.systemStates.A1Site;
        v3(time + 1) = sS.systemStates.B1Site;
        v4(time + 1) = sum(sP.populationOperators.nL1 * sS.systemStates.LHSystem);
        v5(time + 1) = sum(sP.populationOperators.nH1 * sS.systemStates.LHSystem);
        v6(time + 1) = sum(sP.populationOperators.Q1n1 * sS.systemStates.Quinone1);
        v7(time + 1) = sum(sP.populationOperators.Q1n2 * sS.systemStates.Quinone1);
        v8(time + 1) = sum(sP.populationOperators.Q1N1 * sS.systemStates.Quinone1);
        v9(time + 1) = sum(sP.populationOperators.Q1N2 * sS.systemStates.Quinone1);
        v10(time + 1) = sS.systemStates.L2Site;
    end
    
    fileID = fopen('testout.txt','w');
    fprintf(fileID, 'nA\tnB\tnL\tnH\tnL2\r\n');
    for time = 2:sP.nOS
        fprintf(fileID, '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\r\n', v2(time), v3(time), v4(time), v5(time), v10(time));
    end
    fclose(fileID);
    
    dynamicsOutput = struct(f1, v1, f2, v2, f3, v3, f4, v4, ...
        f5, v5, f6, v6, f7, v7, f8, v8, f9, v9, f10, v10);
    
    cumulativeOutput = 0;
end