function  [min_point,xk] = newton_method(varargin)
%clc

% Default syntax for the function is given as follows;
% newton_method(@your_cost_function,[1 2 3 4 5]) , [1 2 3 4 5] represents
% your start point

delta=10^-6;   % Convergence criteria for the norm of the gradient
epsilon=10^-10;  % Convergence criteria for the consecutive cost values
alpha=0.1;         % Step size coefficient
sigma=10^-5;      % Armijo rule multiplier
grad_func=1;
hessian_func=1;

if (nargin<2)
    disp('Missing parameters')
end
cost_func=varargin{1};
start_point=varargin{2};
if nargin>2
    start_point=varargin{2}(varargin{3});
end
if nargin>3
    delta=varargin{4};
end
if nargin>4
    epsilon=varargin{5};
end
if nargin>5
    alpha=varargin{6};
end
if nargin>6
    sigma=varargin{7};
end
if nargin>7
    grad_func=varargin{8};
end
if nargin>8
    hessian_func=varargin{9};
end

dim=length(start_point);
gradient=ones(1,dim);
hessian=(zeros(dim));
hess_diag=zeros(1,dim);
hess_up=zeros(dim);
numb_of_iter=0;
xk=start_point;
error=realmax;
xkmax=ones(1,dim)*-realmax;
xkmin=ones(1,dim)*realmax;

while( (sum(abs(gradient)) > delta) && (abs(error)>epsilon) && (abs(cost_func(xk))<(10^16)) && numb_of_iter<(dim*10^3))
    recent_cost=cost_func(xk);
    if (grad_func==1)
        for i=1:dim
            diff=zeros(1,dim);
            diff(i)=(10^(log10(eps)*0.5));
            gradient(i)=(cost_func(xk+diff)-cost_func(xk-diff))/(2*diff(i));  % Central difference differentiation formula for gradient
        end
    else
        gradient=grad_func(xk);
    end
    if (hessian_func==1 )
        if (grad_func==1)
            for j=1:dim
                diff1=zeros(1,dim);
                diff1(j)=(10^(log10(eps)*0.25));
                hess_diag(j)=(cost_func(xk+diff1+diff1)-(recent_cost*2)+cost_func(xk-diff1-diff1))/(4*diff1(j)*diff1(j));
                for i=(j):dim
                    diff2=zeros(1,dim);
                    diff2(i)=(10^(log10(eps)*0.25));
                    hess_up(j,i)=(cost_func(xk+diff2+diff1)-cost_func(xk+diff2-diff1)-cost_func(xk+diff1-diff2)+cost_func(xk-diff1-diff2))/(4*diff2(i)*diff1(j));
                    % Central difference formula for partial derivative in Hessian
                end
            end
            hessian=hess_up+hess_up'-diag(hess_diag);
        else
            for i=1:dim
                hessian(i,:)=grad_func(ones(1,dim)*gradient(i));
            end
        end
    elseif (hessian_func~=1 )
        hessian=hessian_func(xk);
    end
    
    [upper,p]=chol(hessian);
    n=0;
    if p>0
        disp('Hessian correction')  % if hessian is negative definite it must be corrected
        hessian       
    end
    c=cond(hessian);
    while (p>0 && c~=Inf && isnan(c)==0)   %accumulate the diagonal elements of hessian with a scalar c^n, where c is the condition number
        c=cond(hessian);        
        n=n+0.1;
        hessian=hessian+(c^n)*eye(dim);
        [upper,p]=chol(hessian);
    end    
    if length(upper)<dim
        upper=zeros(dim);
    end
    lower=upper';   % apply cholesky factorization to find step size directly
    m=-8;           % power term for armijo rule which also provides momentum since being negative   
    step_size= -1*((upper\(lower\gradient')))*(alpha^m);
    while ((recent_cost - cost_func(xk+step_size'))<(-sigma*gradient*step_size) && abs(-sigma*gradient*step_size)>eps)
        m=m+1;
        step_size=step_size*alpha;  % apply armijo rule to make step size adaptive
    end
    % display the parameters of recent iteration
    gradient
    hessian
    if isnan(step_size)
        step_size=zeros(dim,1);
    end
    step_size
    xk=step_size+xk';
    xk=xk'
    cost=cost_func(xk)
    error=cost-recent_cost
    numb_of_iter=numb_of_iter + 1
    disp('------------------------------------------------')
    
end


gradient
xk
min_point=(cost_func(xk))
error
numb_of_iter
start_point
if (sum(abs(gradient))>delta)
    % disp('WARNING: GRADIENT CRITERIA IS NOT SATISFIED')
end
if (error>epsilon)
    % disp('WARNING: FUNCTION CONVERGENCE CRITERIA IS NOT SATISFIED')
end
end
