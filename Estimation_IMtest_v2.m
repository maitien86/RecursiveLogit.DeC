function [] = Estimation_IMtest_v2(nObs, nVariable) 

    globalVar;
    global isLinkSizeInclusive;
    global file_observations;
    global isFixedUturn;
    global Obs;
    tic;
    
    file_linkIncidence = './Input/linkIncidence.txt';
    file_AttEstimatedtime = './Input/ATTRIBUTEestimatedtime.txt';
    file_turnAngles = './Input/ATTRIBUTEturnangles.txt';
    
    global RealOD;
    RealOD = spconvert(load('./Input/observationsForEstimBAI.txt'));
    RealOD = RealOD(:,1:2);
    nObs = str2num(nObs);
    if str2num(nVariable) == 4
        isLinkSizeInclusive = false;
    else
        isLinkSizeInclusive = true;
    end
    file_observations = './Input/observationsForEstimBAI.txt';
    
    loadData;
    
    Op = Op_structure;
    initialize_optimization_structure();
 
    %% Generate Obs 
    
    isFixedUturn = false;
    if isLinkSizeInclusive ==  false
        x0 = [-2.0, -1.0, -1.0, -20.0]';
        Op.n = 4;
    else
        x0 = [-2.0, -1.0, -1.0, -20.0, -1.5 ]';
        Op.n = 5;
    end
    if nObs <= 1832
        nbobsOD = 1;
        ODpairs = RealOD(1:nObs,:);
    else 
        nbobsOD = 10;
        ODpairs = RealOD;
    end
    Obs = generateObs([], x0, ODpairs, nbobsOD);
    [nbobs, maxstates] = size(Obs);
    
    isFixedUturn = true;
    isLinkSizeInclusive = false;
    initialize_optimization_structure();
    Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
    Op.Hessian_approx = OptimizeConstant.BHHH;
    Gradient = zeros(nbobs,Op.n);
    if isLinkSizeInclusive == false
        Op.x = [-1.8, -1.0, -1.0]';
    else
        Op.x = [-1.8, -1.0, -1.0, -1.5]';
    end

    % Optimize with full choice set
    EstimationScript;
    
    TXT = ['RL|',num2str(nObs),'|', num2str(Op.n),'|',nVariable,'|' ];
    Op.x
    %TXT = [TXT,IMtest(Op.x)];
    TXT = [TXT,IM_fulltest(Op.x),'\n'];
    fileID = fopen('IMtestResults.txt','at+');
    fprintf(fileID,TXT);
    fclose(fileID);
    
    fileID = fopen('ComputationalTime.txt','at+');
    fprintf(fileID,['IM|',TXT,'Elapsed time:', num2str(toc),'\n']);
    fclose(fileID);
end
