function [D,X,Err,A,B,CC]= swbDL(Y,Dp,Xp,nIter,K,spa2,TC,SM)
    A = eye(size(Dp,2),K);
    B = zeros(K,size(Xp,1));
    D = Dp*A;
    X = B*Xp;
%     fprintf('Iteration:     ');    
    U = zeros(size(D));   
    W = U;
    
    beta = 10^-3; %1e-4
    beta_bar = beta * 10^9; %10^15 
    rho = 1.7; %1.7;
    
    alpha = 0.75; %100                   %100  %5000  %10000
    for iter=1:nIter
        A = 0; B = 0; xk=zeros(size(X));
%         fprintf('\b\b\b\b\b%5i',iter);   
        Dold = D;      
               
        for i =1:K
            xk(i,:) = sign(D(:,i)'*Y).*max(0, bsxfun(@minus,abs(D(:,i)'*Y),spa2/2));            
        end
%         tmp1 = inv(D'*D+ 1e-8*speye(size(D,2)))*D'*Y; 
%         xk = sign(tmp1).*max(0, bsxfun(@minus,abs(tmp1),spa2/2));
        B = xk*(Xp'/(Xp*Xp'));
        X = B*Xp;
                
        Do = filter([0 1],1,D);
        for k = 1:1000
            D = (beta*U + Y*X' - W)/(X*X' + beta*eye(size(D,2)));  D = normc(D);
%             U = (beta+ eta)\(W + beta*D+eta*tt); U = normc(U);
            U = (beta*eye(size(D,1))+ alpha*(U*U')-alpha*(Do*Do'))\(W + beta*D); U = normc(U);
            W = W + beta*(D-U);
            beta = min(beta*rho,beta_bar);
            if norm((D - U),'fro') < 0.0001;  break;  end
        end
        
        A = (Dp'*Dp+ 1e-8*speye(size(Dp,2)))\Dp'*D;
        A = (A*diag(1./sqrt(sum((Dp*A).*(Dp*A))))); 
        D = Dp*A;
        
        Err(iter) = (sqrt(trace((D-Dold)'*(D-Dold)))/sqrt(trace(Dold'*Dold)));       
               
        [~,~,ind]=sort_TSandSM_spatial(TC,SM,D,X,K);
        for ii =1:K
            TCcorr(ii) =abs(corr(TC(:,ii),D(:,ind(ii))));
            SMcorr(ii) =abs(corr(SM(ii,:)',X(ind(ii),:)'));
        end
        cTC = sum(TCcorr');
        cSM = sum(SMcorr');
        CC(iter) =cTC+cSM;
    end
end


