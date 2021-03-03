function [A,E,B]=TensorEnsemble(F0,W0,lambda)
% F0 is the coherent link matrix
% W0 is the coassociation matrix


tol = 1e-8; 
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;

n=length(F0);

sX = [n, n, 2];

B=zeros(n);

C=zeros(n);

E=B;
Lambda1=B;
Lambda2=B;
Lambda3=B;

A=zeros(n,n,2);
T=A;





for iter=1:max_iter
    % update A
    T(:,:,1)=1*(B-Lambda1/mu);
    T(:,:,2)=0.5*(W0-E- Lambda2/mu + C- Lambda3/mu);
    %%%%%%%%%%%%%%%%%%%%%%%
   
    t = T(:);
    
    [a, ~] = wshrinkObj(t,1/mu,sX,0,3)   ;
    A= reshape(a, sX);
    


    %%%%%%%%%%%%%%%%%%%%%%%
    
    % update E
    Temp=W0-A(:,:,2)-Lambda2/mu;
%     E=prox_l1(Temp,(1*lambda)/(mu));
    E=Temp*mu/(2*lambda+mu);
    
    

    
    % update B
    Temp=A(:,:,1)+Lambda1./mu;
    B=0.5*(Temp+Temp');
    B(F0==1)=1;
%     B(F0==-1)=-1;
    B(B<0)=0;
    B(B>1)=1;
    
    
    % update C
    Temp=A(:,:,2)+Lambda3./mu;
    C=0.5*(Temp+Temp');
%     B(F0==1)=1;
%     B(F0==-1)=-1;
    C(C<0)=0;
    C(C>1)=1;
    

   
    
    d1=A(:,:,1)-B;
    d2=A(:,:,2)+E-W0;
    d3=A(:,:,2)-C;
    disp(['iter: ',num2str(iter)])
    
    chg = max([ max(abs(d1(:))),max(abs(d2(:))),max(abs(d3(:)))]);
    if chg < tol
        break;
    end 
    
    Lambda1=Lambda1+mu*d1;
    Lambda2=Lambda2+mu*d2;
    Lambda3=Lambda3+mu*d3;
    mu = min(rho*mu,max_mu); 
end