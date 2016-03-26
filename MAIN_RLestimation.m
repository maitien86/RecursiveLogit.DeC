%  MAIN PROGRAM
Credits;
globalVar;
global resultsTXT; 

%%------------------------------------------------------

file_linkIncidence = './Input/linkIncidence.txt';
file_AttEstimatedtime = './Input/TravelTime.txt';
file_turnAngles = './Input/TurnAngles.txt';
file_observations = './simulatedData/ObservationsAll.txt';


% Initialize the optimizer structure
isLinkSizeInclusive = false;
isFixedUturn = false;
loadData;
Op = Op_structure;
initialize_optimization_structure();
Op.x = [-2.494,-0.933,-0.411,-4.459]';
Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
Op.Hessian_approx = OptimizeConstant.BHHH;
Gradient = zeros(nbobs,Op.n);

%% Starting optimization
tic ;
disp('Start estimating ....')

% print result to string text
header = [sprintf('%s \n',file_observations) Op.Optim_Method];
header = [header sprintf('\nNumber of observations = %d \n', nbobs)];
header = [header sprintf('Hessian approx methods = %s \n', OptimizeConstant.getHessianApprox(Op.Hessian_approx))];
resultsTXT = header;
%------------------------------------------------
options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton','GradObj','on');
[x,fval,exitflag,output,grad] = fminunc(@LL,Op.x,options);
% while (true)    
%   Op.k = Op.k + 1;
%   if strcmp(Op.Optim_Method,OptimizeConstant.LINE_SEARCH_METHOD);
%     ok = line_search_iterate();
%     if ok == true
%         PrintOut(Op);
%     else
%         disp(' Unsuccessful process ...')
%         break;
%     end
%   else
%     ok = btr_interate();
%     PrintOut(Op);
%   end
%   [isStop, Stoppingtype, isSuccess] = CheckStopping(Op);  
%   %----------------------------------------
%   % Check stopping criteria
%   if(isStop == true)
%       isSuccess
%       fprintf('The algorithm stops, due to %s', Stoppingtype);
%       resultsTXT = [resultsTXT sprintf('The algorithm stops, due to %s \n', Stoppingtype)];
%       break;
%   end
% end
%   Compute variance - Covariance matrix
PrintOut(Op);

disp(' Calculating VAR-COV ...');
getCov;

%% Finishing ...
ElapsedTtime = toc
resultsTXT = [resultsTXT sprintf('\n Number of function evaluation %d \n', Op.nFev)];
resultsTXT = [resultsTXT sprintf('\n Estimated time %d \n', ElapsedTtime)];
