function [D,X,Err,CC]= my_ACSD(Y,Di,spa,nIter,TC,SM)
    D = Di;
    K= size(D,2);
    X = zeros(size(D,2),size(Y,2)); 
%     fprintf('Iteration:     ');
    for iter=1:nIter
%         fprintf('\b\b\b\b\b%5i',iter);
        Dold = D;
        for j =1:size(D,2)
            X(j,:) = 0;
            E = Y-D*X;
            xk = D(:,j)'*E; 
            thr = spa./abs(xk);
            X(j,:) = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr/2));
            rInd = find(X(j,:));
            if (length(rInd)<1)
%                 D(:,j)= randn(size(D(:,j),1), 1);
                [~,ind]= max(sum(Y-D*X.^2)); 
                D(:,j)= Y(:,ind)/norm(Y(:,ind));
            else
                D(:,j) = E(:,rInd)*X(j,rInd)'./norm(E(:,rInd)*X(j,rInd)');
            end                 
        end      
        Err(iter) = sqrt(trace((D-Dold)'*(D-Dold)))/sqrt(trace(Dold'*Dold));

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
