function [X2,E]=MNSR(A,B)
    opts.tol = 1e-8; 
    opts.max_iter = 500;
    opts.rho = 1.1;
    opts.mu = 1e-4;
    opts.max_mu = 1e10;
    opts.DEBUG = 0;
    lambda = 0.01;
    opts.loss = 'l1'; 
    [X2,E,~,~,~] = SNR(A,B,lambda,opts);
end


function [X,E,obj,err,iter] = SNR(A,B,lambda,opts)
tol = 1e-8; 
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;
DEBUG = 0;
loss = 'l1';
k=10;

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'loss');        loss = opts.loss;            end
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end


[d,na] = size(A);
[~,nb] = size(B);

I1=ones(na,1);
I2=ones(nb,1);
A=[A;I1'];
B=[B;I2'];
[d,na] = size(A);
[~,nb] = size(B);

X = zeros(na,nb);
E = zeros(d,nb);
Z = X;
Y1 = E;
Y2 = X;

AtB = A'*B;
I = eye(na);
invAtAI = (A'*A+I)\I;

iter = 0;
for iter = 1 : max_iter
    Xk = X;
    Ek = E;
    Zk = Z;
    % first super block {X,E}
    X = prox_l1(Z-Y2/mu,lambda/mu);
    X=max(X,0);
    if strcmp(loss,'l1')
        E = prox_l1(B-A*Z-Y1/mu,1/mu);
        E(end,:)=0;
    elseif strcmp(loss,'l2')
        E = mu*(B-A*Z-Y1/mu)/(1+mu);
        E(end,:)=0;
    else
        error('not supported loss function');
    end
    % second  super block {Z}
    Z = invAtAI*(-A'*(Y1/mu+E)+AtB+Y2/mu+X);    
    dY1 = A*Z+E-B;
    dY2 = X-Z;
    chgX = max(max(abs(Xk-X)));
    chgE = max(max(abs(Ek-E)));
    chgZ = max(max(abs(Zk-Z)));
    chg = max([chgX chgE chgZ max(abs(dY1(:))) max(abs(dY2(:)))]);
    if DEBUG        
        if iter == 1 || mod(iter, 10) == 0
            obj = comp_loss(E,loss)+lambda*norm(X(:),1);
            err = sqrt(norm(dY1,'fro')^2+norm(dY2,'fro')^2);
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    if chg < tol
        break;
    end 
    Y1 = Y1 + mu*dY1;
    Y2 = Y2 + mu*dY2;
    mu = min(rho*mu,max_mu);    
end
obj = comp_loss(E,loss)+lambda*norm(X(:),1);
err = sqrt(norm(dY1,'fro')^2+norm(dY2,'fro')^2);
end