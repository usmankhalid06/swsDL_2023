function [D,X,Err,A,B,CC]= swsDL(Y,Dp,Xp,nIter,K,lam,zeta,TC,SM)
    A = eye(size(Dp,2),K);
    B = zeros(K,size(Xp,1));
    D = Dp*A;
    X = B*Xp;
    Err = zeros(1,nIter);
    fprintf('Iteration:     ');    
    for iter=1:nIter
        fprintf('\b\b\b\b\b%5i',iter);   
        Dold = D;      
        D = Dp*A;
        X = B*Xp;
        for j=1:size(D,2)
            X(j,:) = 0; A(:,j) = 0; B(j,:) = 0;
            E = Y-D*X;
            xk = D(:,j)'*E;
            thr = zeta./abs(xk);
            
            xkk = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr/2));
            [~,bb1]= sort(abs(Xp*xkk'),'descend');
            ind1 = bb1(1:lam);
            B(j,ind1)= xkk*Xp(ind1,:)'/(Xp(ind1,:)*Xp(ind1,:)');
            X(j,:) = B(j,:)*Xp;
            
            rInd = find(X(j,:));
            if (length(rInd)<1)
                [~,ind]= max(sum(Y-Dp*A*X.^2));
                D(:,j)= Y(:,ind)/norm(Y(:,ind));
            else
                mag = X(j,rInd)*X(j,rInd)';
                tmp3 = D(:,j)+(1/mag)*E(:,rInd)*X(j,rInd)'; 
            
                [~,bb2]= sort(abs(Dp'*tmp3),'descend');
                ind2 = bb2(1:lam);
                A(ind2,j)= (Dp(:,ind2)'*Dp(:,ind2))\Dp(:,ind2)'*tmp3;         
                
                A(:,j) = A(:,j)./norm(Dp*A(:,j));
                D(:,j) = Dp*A(:,j);
            end
        end
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






