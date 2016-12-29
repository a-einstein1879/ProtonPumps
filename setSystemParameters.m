function sP = setSystemParameters()
    %% Number of each system levels
    NLevQ1 = 16;
    NLevLH = 4;
    NLevRC = 2;
    
    field1 = 'numberOfLevels';
    % number of basis electron/proton states in a first Quinone
    % 4 particles can be or not be there, so 2^4 = 16
    f1 = 'Q1'; v1 = NLevQ1;
    
    % L-H (elements of Q-Cycle): two electron sites, coupled to each other and
    % to e-sites 1,2 of the first Quinone; 
    % Together: 4 basis states for LH system (Feb/2/2011/A-5)
    f2 = 'LH'; v2 = NLevLH;
    
    % number of basis states for A1 and for D1 separately
    % 1 electron can "be" or "not be" in both A1 and D1
    f3 = 'RC'; v3 = NLevRC;
    
    value1 = struct(f1, v1, f2, v2, f3, v3);
    
    %% State matrices
    % CA stands for creation/annihilation
    field2 = 'CAOperators';
    
    % First Quinone matrices
    % creation and annihilation operators for 2 electrons (a)
    % and 2 protons (b) in the first Quinone
    a1 = zeros(NLevQ1, NLevQ1); a2 = zeros(NLevQ1, NLevQ1); 
    b1 = zeros(NLevQ1, NLevQ1); b2 = zeros(NLevQ1, NLevQ1); 

    a1(1,2)=1; a1(3,4)=1; a1(5,8)=1; a1(9,10)=1; a1(6,11)=1;
    a1(12,13)=1; a1(7,14)=1; a1(15,16)=1;

    a2(1,3)=1; a2(2,4)=-1; a2(5,9)=1; a2(8,10)=-1; a2(6,12)=1; 
    a2(11,13)=-1; a2(7,15)=1; a2(14,16)=-1; 

    b1(1,5)=1; b1(6,7)=1; b1(2,8)=1; b1(3,9)=1; b1(4,10)=1; 
    b1(11,14)=1; b1(12,15)=1; b1(13,16)=1; 

    b2(1,6)=1; b2(5,7)=-1; b2(2,11)=1; b2(3,12)=1; b2(4,13)=1; 
    b2(8,14)=-1; b2(9,15)=-1; b2(10,16)=-1;

    f1 = 'Q1a1';  v1 = a1;
    f2 = 'Q1a2';  v2 = a2;
    f3 = 'Q1b1';  v3 = b1;
    f4 = 'Q1b2';  v4 = b2;
    
    % LH-system matrices
    aL = zeros(NLevLH,NLevLH); aH = zeros(NLevLH,NLevLH);
    aL(1,2)=1; aL(3,4)=1;
    aH(1,3)=1; aH(2,4)=-1;
    
    f5 = 'aL1'; v5 = aL;
    f6 = 'aH1'; v6 = aH;
    
    value2 = struct(f1, v1, f2, v2, f3, v3, f4, v4, f5, v5, f6, v6);
    
    %% Population operators
    field3 = 'populationOperators';
    % First quinone electron and proton populations
    % n1 = a1' * a1
    f1 = 'Q1n1'; v1 = a1' * a1;
    % n2 = a2' * a2
    f2 = 'Q1n2'; v2 = a2' * a2;
    % N1 = b1' * b1
    f3 = 'Q1N1'; v3 = b1' * b1;
    % N2 = b2' * b2
    f4 = 'Q1N2'; v4 = b2' * b2;
    
    % LH-system populations
    % nL = aL' * aL
    f5 = 'nL1'; v5 = aL' * aL;
    % nH = aH' * aH
    f6 = 'nH1'; v6 = aH' * aH;
    
    value3 = struct(f1, v1, f2, v2, f3, v3, f4, v4, f5, v5, f6, v6);
    
    %% Transfer distances
    field4 = 'transferDistances';
    % proton transfer length, in nm; Q1-to N or to P-side
    f1 = 'q1ToN'; v1 = 0.25;
    % electron transfer length; Q1-to A1 or B1; and Q1-to L1/H1
    f2 = 'q1ToAB1'; v2 = 0.25;
    value4 = struct(f1, v1, f2, v2);
    
    %% Some general parameters
    % 1 meV = 1520/ns; tau0 = 0.01 microsec
    field10 = 'meVtoTime'; value10 = 1.52 * 10^(4);
    % modeling time in mS
    field11 = 'modelingTime'; value11 = 500;
    numberOfStepsPermS = 20;
    % number of steps
    field12 = 'nOS'; value12 = value11 * numberOfStepsPermS;
    % step size
    field13 = 'dt'; value13 = 1 / numberOfStepsPermS;
    % ????????????????????????
    field14 = 'dM'; value14 = 300;
    % ????????????????????????
    field15 = 'SDChain'; value15 = [1370; 1395; 1420; 1445];
    %value15 = [1220; 1245; 1270; 1295; 1320; 1345; 1370; 1395; 1420; 1445];
    % ????????????????????????
    field16 = 'VSet'; value16 = [160; 260; 360];
    % ????????????????????????
    field17 = 'Uep'; value17 = 800;
    
    %% Temperatures
    % temperature in Kelvin
    field18 = 'TK'; value18 = 298;
    % initial temperature in Kelvin
    field19 = 'TK0'; value19 = value18;
    % temperature in meV
    field20 = 'TT'; value20 = value19 / 11.6;
    
    %% Tunneling amplitudes in mEv
    field21 = 'tunnelingAmplitudes';
    f1 = 'dR';    v1 = 1;
    f2 = 'dA1Q1'; v2 = 0.1  * v1;
    f3 = 'dL1Q1'; v3 = 0.06 * v1;
    f4 = 'dL1H1'; v4 = 0.1  * v1;
    f5 = 'dB1Q1'; v5 = v2;
    f6 = 'dH1Q1'; v6 = v3;
    value21 = struct(f1, v1, f2, v2, f3, v3, f4, v4, f5, v5, f6, v6);
    
    %% Coulomb energies in first Quinone
    field22 = 'columbQuinone1';
    % electron/proton energies in quinone
    f1 = 'Up'; v1 = value17 / 8;
    %(5/8)*Uep; %500;
    f2 = 'ue'; v2 = value17 * 0.5;
    % interaction energies
    f3 = 'u11'; v3 = value17;
    f4 = 'u12'; v4 = value17;
    f5 = 'u21'; v5 = value17;
    f6 = 'u22'; v6 = value17;
    f7 = 'uLH'; v7 = 240;
    value22 = struct(f1, v1, f2, v2, f3, v3, f4, v4, f5, v5, f6, v6, f7, v7);

    %% Energy levels
    field23 = 'energyLevels';
    % Lambda for energy levels only
    LamL1Q1 = 60;
    LamL1H1 = 140; % 260 - 2 * LamL1Q1
    LamA1Q1 = 50;
    LamH1Q1 = LamL1Q1; %100;%200;
    LamB1Q1 = LamA1Q1; %200;%200;
    
    % Parameters: These are energy^(0)
    f1 = 'eB1'; v1 = 0;
    % eQ = eB + 2 * Uep - ue + LamDQ - 50
    f2 = 'eQ1'; v2 = v1 + 2 * value17 - value22.ue + LamB1Q1 - 50;
    % eH = eQ + LamHQ; % - alph*uLH;
    f3 = 'eH1'; v3 = v2 + LamH1Q1;
    % eL = eQ - LamLQ; % - alph*uLH;
    f4 = 'eL1'; v4 = v2 - LamL1Q1;
    % eA = eD + 2*Uep; % + LamDQ;% + LamAQ;
    f5 = 'eA1'; v5 = v1 + 2 * value17;
    value23 = struct(f1, v1, f2, v2, f3, v3, f4, v4, f5, v5);
    
    %% REAL Lambdas
    field24 = 'lambdas';
    f1 = 'LamL1Q1'; v1 = 100;
    f2 = 'LamL1H1'; v2 = 250;
    f3 = 'LamA1Q1'; v3 = 100;
    % LamH1Q1 = LamL1Q1
    f4 = 'LamH1Q1'; v4 = v1;
    % LamDQ = LamAQ
    f5 = 'LamB1Q1'; v5 = v3;
    % ??????????????????????
    f6 = 'LimNum'; v6 = 0.005;
    value24 = struct(f1, v1, f2, v2, f3, v3, f4, v4, f5, v5, f6, v6);
    
    %% Gammas
    field25 = 'gammas';
    % Electron tunneling to Source/L2
    f1 = 'gamS'; v1 = 0.00001*10;
    % gamL2 = gamS
    f2 = 'gamL2'; v2 = v1;
    % Proton rates for Quinone-to-(N/P membrane sides) transitions, in meV
    f3 = 'gamN'; v3 = 0.001*2;
    % GamP = GamN
    f4 = 'gamP'; v4 = v3;
    value25 = struct(f1, v1, f2, v2, f3, v3, f4, v4);
    
    %% First quinone mechanical motion (diffusion)
    field26 = 'q1MM';
    % 2 nm, x = -x0 (N-side); x = + x0 (P-side of membrane)
    f1 = 'x0';                      v1 = 2;
    f2 = 'DifCoef0';                v2 = 0.08;
    % ??????????? do we need to save TK0??
    % delete TK0!
    % FactorDif = sqrt(2*(TK/TK0)*(DifCoef0/dt))
    f3 = 'factorDif';               v3 = sqrt(2 * (value18 / value19) * (v2 / value13));
    % Edif = (TK0/11.6)*(1/DifCoef0)
    f4 = 'Edif';                    v4 = (value19 / 11.6) * (1 / v2);
    
    % membrane potential
    % slope of membrane walls
    f5 = 'membraneWallSlope';       v5 = 0.1;
    % height of membrane walls
    f6 = 'membraneWallHeight';      v6 = 520;
    % Effective width of membrane x00 > x0 = 2 nm, former x00
    f7 = 'membraneEffectiveWidth';  v7 = 2.7;
    
    % penalty for charged shuttle
    % effective width of penalty potential xx0 < x0, former xx0
    f8 = 'penaltyPotentialWidth';   v8 = 1.7;
    % former xxc, slope of penalty potential
    f9 = 'penaltyPotentialSlope';   v9 = 0.05;
    % former Vcc, height of penalty potential
    f10 = 'penaltyPotentialHeight'; v10 = 500;
    
    % ???????????????????
    % ExW = exp(xW/xWk)
    f11 = 'ExW';        v11 = exp(v7 / v5);
    % ???????????????????
    % exCh = exp(xCh/xCk)
    f12 = 'exCh';       v12 = exp(v8 / v9);
    % WALLS, former CFF
    % CVw = (Vw/xWk)*(1/Edif)
    f13 = 'CVw';        v13 = (v6 / v5) * (1 / v4);
    % Penalty for charged quinone, former CFC
    % CVch = (Vch/xCk)*(1/Edif)
    f14 = 'CVch';       v14 = (v10 / v9) * (1 / v4);
    
    value26 = struct(f1, v1, f2, v2, f3, v3, f4, v4, f5, v5, f6, v6, ...
        f7, v7, f8, v8, f9, v9, f10, v10, f11, v11, f12, v12, f13, v13, f14, v14);
    
    % ?????????????????????????
    %randn('state',randNum) %% set the state of randn, keep the same realization
    
    
    %% Delam matrices
    field27 = 'delam';
    % First Quinone-RC transitions, Feb/3/2011/A-3,4
    % root coefficients in A.11 (appendix)
    % DeLamA1Q1 = abs(dA1Q1).^2 * sqrt(pi / (LamA1Q1 * TT))
    f1 = 'DeLamA1Q1'; v1 = abs(value21.dA1Q1).^2 * sqrt(pi / (value24.LamA1Q1 * value20));
    % DeLamB1Q1 = abs(dB1Q1).^2 * sqrt(pi / (LamB1Q1 * TT))
    f2 = 'DeLamB1Q1'; v2 = abs(value21.dB1Q1).^2 * sqrt(pi / (value24.LamB1Q1 * value20));

    % Quinone-LH transitions
    % root coefficients in A.4 (appendix)
    % DeLamL1Q1 = abs(dL1Q1).^2 * sqrt(pi / (LamL1Q1 * TT))
    f3 = 'DeLamL1Q1'; v3 = abs(value21.dL1Q1).^2 * sqrt(pi / (value24.LamL1Q1 * value20));
    % DeLamH1Q1 = abs(dH1Q1).^2 * sqrt(pi / (LamH1Q1 * TT))
    f4 = 'DeLamH1Q1'; v4 = abs(value21.dH1Q1).^2 * sqrt(pi / (value24.LamH1Q1 * value20));
    
    value27 = struct(f1, v1, f2, v2, f3, v3, f4, v4);

    %% sP for system parameters
    % store all parameters in same structure
    sP = struct(field1, value1, field2, value2, field3, value3, field4, value4, ...
        field10, value10, field11, value11, ...
        field12, value12, field13, value13, field14, value14, field15, value15, ...
        field16, value16, field17, value17, field18, value18, field19, value19, ...
        field20, value20, field21, value21, field22, value22, field23, value23, ...
        field24, value24, field25, value25, field26, value26, field27, value27);
end