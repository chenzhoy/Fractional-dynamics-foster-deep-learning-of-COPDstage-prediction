function [Aout,B,order,u,relErr] = modelEst(varargin)

% input format:
% sensInd: sensors index needs to be selected out of total state dimensions
%          could be set as [1:length(data)] if all states need to be
%          selected for modeling
%
% numInp: unknown input dimensions (must be less than length(sensInd) )
% 
% data: input data of size NxD, where 'N' is number of samples and 'D' is
% state dimension
% 
% silentFlag: for toggling script intermediary status outputs

% output:
% Aout: estimated spatial coupling matrix
% B: input coupling matrix
% order: fractional order vector of size length(sensInd)
% u: estimated unknown inputs delayed by one step
% relErr: relative error vector between data and predicted model with numStep set
% in the script. First error is without taking unknown inputs and second
% error is with taking unknown inputs into consideration


global infit
global numCh
% addpath(genpath(pwd));

nargin = length(varargin);
if nargin > 1
    for i = 1:2:nargin-1
        switch varargin{i}
            case 'sensInd'
                sensInd = varargin{i+1};
            case 'numInp'
                numInp = varargin{i+1};
            case 'data'
                data = varargin{i+1};
            case 'silentFlag'
                silentFlag = varargin{i+1};
        end
    end
%     if ~isempty(patInpEdfStr)
%         [~,data] = edfread(patInpEdfStr);
%     else
% %         patInpStr = 'Data/EEG/S001/mat/S001R03_edfm.mat';
% %         matObj = matfile(patInpStr);
% %         data = matObj.record;
%     end 
    
    calledFromOutside = 1;
else
    sensInd = 1:64;
    numInp = 32;
%     patInpStr = 'Data/EEG/S001/mat/S001R03_edfm.mat';
%     matObj = matfile(patInpStr);
    
    data = [];
    if isempty(data)
       error('default data has to be set here');
    end
    silentFlag = 0;
    calledFromOutside = 0;
end

numCh = length(sensInd);

% K = 400;
% sampleID = [0:K-1]+5000;
% X = data(sensInd,sampleID);
K = size(data,1);
X = data';

order = zeros(numCh,1);
infit = 20;

yCascade = zeros(numCh, K); % [y[0], y[1], ..., y[K-1]]

for i = 1:numCh
    order(i,:) = WT_estimator_v3(X(i,:),1);
    yCascade(i,:) = getFractionalExpan(X(i,1:K),order(i,:),infit);
end

% order = 1.2*ones(1,length(order));

niter = 10;
A = cell(niter,1);

[A{1},mse] = performLeastSq(yCascade, X(:,1:K));

% B = double(A{1} > 0.01); %assumption right now
B = zeros(size(A{1},1),size(A{1},2));
B(abs(A{1})>0.01) = A{1}(abs(A{1})>0.01);
[~,r] = qr(B);
colInd = find(abs(diag(r))>1e-7);
if length(colInd) < numInp
    B = [eye(numInp);zeros(numCh-numInp,numInp)];
else
    colInd = colInd(1:numInp);
    % colInd = [1:4,6:26,28:34];
    B = B(:,colInd);
end
% B = [eye(numInp);zeros(numCh-numInp,numInp)];
if rank(B) < numInp
    error('rank deficient B');
end

u = zeros(size(B,2), K);
if ~silentFlag,fprintf('before iteration, mse = %f\n', mse);end
mseIter = zeros(niter+1,1);
mseIter(1) = mse;
for iterInd = 1:niter
    for kInd = 2:K
        yUse = yCascade(:,kInd) - A{iterInd}*X(:,kInd-1);
        u(:,kInd) = getLassoSoln(B, yUse, 0.5);
    end
    [A{iterInd+1},mseIter(iterInd+1)] = performLeastSq(yCascade - B*u, X(:,1:K));
    
    if ~silentFlag,fprintf('iter ind = %d, mse = %f\n', iterInd, mseIter(iterInd+1));end
    
%     X = X - B*u;
%     for i = 1:numCh
%         order(i) = WT_estimator_v3(X(i,:),1);
%         yCascade(i,:) = getFractionalExpan(X(i,1:K),order(i),infit);
%     end
end
if ~calledFromOutside
    figure;
    plot(1:niter+1, mseIter);grid;
end
Aout = A{end};

% predict values for models, k-step prediction
T = K;
numStep = 1;
chUse = 1;
% without inputs
xPred = predictValues(X, order, T, numStep, A{1}, B, zeros(size(B,2),K));
% mean squared error across all channels 
relErr1 = sqrt(sum(sum((xPred-X).^2))/T/numCh);

% figure;
% subplot(2,1,1);
% plot(1:T, X(chUse,1:K), 'b', 1:T, xPred(chUse,1:T), 'r');grid;
% title(sprintf('model without inputs'));
% legend({'observed', 'predicted'});

% with inputs
xPred = predictValues(X, order, T, numStep, A{end}, B, u);
relErr2 = sqrt(sum(sum((xPred-X).^2))/T/numCh);

% subplot(2,1,2);
% plot(1:T, X(chUse,1:K), 'b', 1:T, xPred(chUse,1:T), 'r');grid;
% title(sprintf('model with inputs'));
% legend({'observed', 'predicted'});


if ~calledFromOutside
    figure;
    plot(sampleID, X(chUse,1:T), 'b', sampleID, xPred(chUse,1:T), 'r');grid;
    legend({'observed', 'predicted'});
end

fprintf('rel Err without inputs = %f\n', relErr1);
fprintf('rel Err with inputs = %f\n', relErr2);
relErr = [relErr1,relErr2];


function [xPred] = predictValues(X, order, T, numStep, A, B, u)

global infit
global numCh

TSteps =  ceil(T/numStep);

xPred = zeros(numCh,T);
xPred(:,1:numStep) = X(:,1:numStep) + B*u(:,1:numStep);

for i = 2:TSteps
    XTemp = zeros(numCh,T);
    XTemp(:,1:(i-1)*numStep) = X(:,1:(i-1)*numStep);
    for stepInd = 1:numStep
        for chInd = 1:numCh
            alpha = order(chInd,:);
            trailLen = min(infit,(i-1)*numStep+stepInd-1);
            j = 1:trailLen;
            if length(alpha)>1
                for aInd = 1:length(alpha)
                    preFact = gamma(-alpha(aInd)+j)./(gamma(-alpha(aInd)).*gamma(j+1));
                    XTemp(chInd,(i-1)*numStep + stepInd) = XTemp(chInd,(i-1)*numStep + stepInd) ...
                        - XTemp(chInd,(i-1)*numStep + stepInd-j)*preFact';
                end
                XTemp(chInd,(i-1)*numStep + stepInd) = XTemp(chInd,(i-1)*numStep + stepInd) / length(alpha);
            else
                preFact = gamma(-alpha+j)./(gamma(-alpha).*gamma(j+1));
                XTemp(chInd,(i-1)*numStep + stepInd) = XTemp(chInd,(i-1)*numStep + stepInd) ...
                    - XTemp(chInd,(i-1)*numStep + stepInd-j)*preFact';
            end
        end
        XTemp(:,(i-1)*numStep + stepInd) = XTemp(:,(i-1)*numStep + stepInd) ...
            + A*XTemp(:,(i-1)*numStep + stepInd-1) ...
            + B*u(:,(i-1)*numStep + stepInd)+ randn(numCh,1);
%         u is already arranged in one step behind order, so no need for
%         '-1' in the index
    end
    xPred(:,(i-1)*numStep+1:i*numStep) = XTemp(:,(i-1)*numStep+1:i*numStep);
end
% relErr = sqrt(mean(abs(xPred(chUse,:) - X(chUse,:)).^2));


function out = getLassoSoln(A, b, lambda)


% method :
% 1: CVX
% 2: lasso MATLAB
% 3: ADMM

method = 3;
switch method
    
    case 1
        cvx_begin quiet
            variable uGG(32)
            minimize(sum_square(A*uGG - b) + lambda*norm(uGG,1))
        cvx_end
        out = uGG;
    case 2
        [out, fInfo] = lasso(A, b, 'Lambda', lambda);

% problem: it doesn't allow to remove the intercept
        
    case 3
%         downloaded from: https://web.stanford.edu/~boyd/papers/admm/lasso/lasso.html 
        % QUIET = 1;
        MAX_ITER = 100;
        ABSTOL   = 1e-4;
        RELTOL   = 1e-2;

        [m, n] = size(A);
        Atb = A'*b;
        % lambda = 1;
        rho = 1/lambda;
        alpha = 1;

        x = zeros(n,1);
        z = zeros(n,1);
        u = zeros(n,1);

        [L, U] = factor(A, rho);

        for k = 1:MAX_ITER

            % x-update
            q = Atb + rho*(z - u);    % temporary value
            if( m >= n )    % if skinny
               x = U \ (L \ q);
            else            % if fat
               x = q/rho - (A'*(U \ ( L \ (A*q) )))/rho^2;
            end

            % z-update with relaxation
            zold = z;
            x_hat = alpha*x + (1 - alpha)*zold;
            z = shrinkage(x_hat + u, lambda/rho);

            % u-update
            u = u + (x_hat - z);

            % diagnostics, reporting, termination checks
            history.objval(k)  = objective(A, b, lambda, x, z);

            history.r_norm(k)  = norm(x - z);
            history.s_norm(k)  = norm(-rho*(z - zold));

            history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
            history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);

        %     if ~QUIET
        %         fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
        %             history.r_norm(k), history.eps_pri(k), ...
        %             history.s_norm(k), history.eps_dual(k), history.objval(k));
        %     end

            if (history.r_norm(k) < history.eps_pri(k) && ...
               history.s_norm(k) < history.eps_dual(k))
                 break;
            end

        end
        out = z;
end


function p = objective(A, b, lambda, x, z)

p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) );


function z = shrinkage(x, kappa)
    
z = max( 0, x - kappa ) - max( 0, -x - kappa );

function [L, U] = factor(A, rho)
    
[m, n] = size(A);
if ( m >= n )    % if skinny
   L = chol( A'*A + rho*speye(n), 'lower' );
else            % if fat
   L = chol( speye(m) + 1/rho*(A*A'), 'lower' );
end

% force matlab to recognize the upper / lower triangular structure
L = sparse(L);
U = sparse(L');


function [A,mse] = performLeastSq(yCascade, X)

numCh = size(X,1);

XUse = [zeros(1,numCh);X(:,1:end-1)'];
% The data may not be zero mean, so take care of intercept as well
% XUse = [ones(size(XUse,1),1),XUse];

A = zeros(numCh,numCh);
mse = zeros(numCh,1);

for i = 1:numCh
    A(i,:) = regress(yCascade(i,:)', XUse);
%     A(i,:) = getLassoSoln(XUse, yCascade(i,:)', 0.5);
    mse(i) = norm(yCascade(i,:)' - XUse*A(i,:)',2)^2 / size(yCascade,2);
end
mse = mean(mse);


function y = getFractionalExpan(x, alpha, infit)

l = length(x);
j = 0:infit;
if length(alpha)>1
    y = zeros(1,l);
    for i = 1:length(alpha)
        preFactVec = gamma(-alpha(i)+j)./(gamma(-alpha(i)).*gamma(j+1));
        yConv = conv(x, preFactVec);
        y = y + yConv(1:l);
    end
    y = y/length(alpha);
else
    preFactVec = gamma(-alpha+j)./(gamma(-alpha).*gamma(j+1));
    y = conv(x, preFactVec);
    y = y(1:l);
end
