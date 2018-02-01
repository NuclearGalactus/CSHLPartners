function model = kalmanRW(X,r,param)
    
    if nargin < 4; param = KTD_defparam; end
    
    % initialization
    [N,D] = size(X);
    w = zeros(D,1);
    
    % parameters
    if nargin < 3 || isempty(param); param = KTD_defparam; end
    C = param.c*eye(D); % prior covariance
    s = param.s;        % noise variance
    q = param.q;        % transition variance
    I = eye(D);
    
    if length(s)==1; s = zeros(N,1)+s; end
    if length(q)==1; q = zeros(N,1)+q; end
    if length(param.lr)==1; param.lr = zeros(N,1)+param.lr; end
    % run Kalman filter
    lastPerr = 0;
    
    for n = 1:N
        h = zeros(1,2);
        h(1:2) = X(n,1:2);
        rhat = h*w;
       
        dt = (r(n)) - rhat;
                %dPerr = dt - lastPerr;
        %h(1,3) = dPerr;
        
                 % prediction error
       
        C = C + q(n)*I;             % a priori covariance
        P = h*C*h'+s(n);            % residual covariance
        K = C*h'/P;                 % Kalman gain
        w0 = w;
        if param.std
            w = w + param.lr(n)*h'*dt;
        else
            w = w + K*dt;               % weight update
        end

       % lastPerr = dt;
        C = C - K*h*C;              % posterior covariance update
        
        % store results
        model(n).w0 = w0;
        model(n).C = C;
        model(n).K = K;
        model(n).dt = dt;
        model(n).rhat = rhat;
        
    end