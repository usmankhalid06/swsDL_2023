function my_subplots_horz(nA,S,v,w,TCcorr,SMcorr,rTC,rSM)

    N = size(rTC{1},1);
    axis off
    set(gca,'Units','normalized','Position',[0 0 1 1]);

    text(0.01,0.94, '(A)','Color','b','FontSize',12)
    
    text(0.01,0.80, ['\Sigma |\rho| =' num2str(round(sum(SMcorr(2,:)),2),'%0.2f') ],'Color','m','FontSize',12)
    text(0.01,0.75, '(B)','Color','b','FontSize',12)
    text(0.01,0.70, ['\Sigma |\rho| =' num2str(round(sum(TCcorr(2,:)),2),'%0.2f') ],'Color','m','FontSize',12)
   
    text(0.01,0.59, ['\Sigma |\rho| =' num2str(round(sum(SMcorr(3,:)),2),'%0.2f') ],'Color','m','FontSize',12)
    text(0.01,0.54, '(C)','Color','b','FontSize',12)
    text(0.01,0.49, ['\Sigma |\rho| =' num2str(round(sum(TCcorr(3,:)),2),'%0.2f') ],'Color','m','FontSize',12)

    text(0.01,0.38, ['\Sigma |\rho| =' num2str(round(sum(SMcorr(4,:)),2),'%0.2f') ],'Color','m','FontSize',12)
    text(0.01,0.33, '(D)','Color','b','FontSize',12)
    text(0.01,0.28, ['\Sigma |\rho| =' num2str(round(sum(TCcorr(4,:)),2),'%0.2f') ],'Color','m','FontSize',12)

    text(0.01,0.17, ['\Sigma |\rho| =' num2str(round(sum(SMcorr(5,:)),2),'%0.2f') ],'Color','m','FontSize',12)
    text(0.01,0.13, '(E)','Color','b','FontSize',12)
    text(0.01,0.07, ['\Sigma |\rho| =' num2str(round(sum(TCcorr(5,:)),2),'%0.2f') ],'Color','m','FontSize',12)
     
       
    ivs = 1.54;  %initial_horizontal_shift
    ihs = 1.35;  %initial_vertical_shift
    shz = S-1.55; %subplot_horizontal_size
    svs = 0.08;  %subplot_vertical_size (more the better)
    vs  = 0.6;  %vertical shift of subplots
    nC  = S+1.5; %No. of columns
    shifter = 0.35; %vertical


    for i =1:nA 
        for j=1:S

            if i==1
            %%
                hax=axes();
                imagesc(flipdim(reshape(abs(zscore(rSM{i}(j,:))),v,w),1)); 
                newPos=[(1-1/nC)-(1/nC)*(fix((j-1)/1)+ihs-1),  vs*(mod(j-1,1)+ivs+0.00),     1/shz,   svs];
                set(gca,'outer',newPos), 
                set(gca,'XTickLabel','')
                set(gca,'YTickLabel','')
                
                hax=axes();
                plot(zscore(rTC{i}(:,j)));  axis([0 N -3 3]);
                newPos=[(1-1/nC)-(1/nC)*(fix((j-1)/1)+ihs-1),  vs*(mod(j-1,1)+ivs-0.08),     1/shz,   svs-0.025];
                set(gca,'outer',newPos),
                set(gca,'XTickLabel','')
                set(gca,'YTickLabel','')
                
            else

            %%
                zscore_rxSM = abs(zscore(rSM{i}(j,:)));
                hax=axes();
                imagesc(flipdim(reshape(zscore_rxSM,v,w),1));  
                newPos=[(1-1/nC)-(1/nC)*(fix((j-1)/1)+ihs-1),  vs*(mod(j-1,1)+ivs-shifter*(i-1))+0.04,     1/shz,   svs];
                set(gca,'outer',newPos), 
                set(gca,'XTickLabel','')
                set(gca,'YTickLabel','')
                xlabel(['|\rho|',' = ',num2str(round(SMcorr(i,j),2))],'color','k')
                xh = get(gca,'xlabel'); % handle to the label object
                p = get(xh,'position'); % get the current position property
                p(2) = - 7.5+ p(2);       % double the distance,
                set(xh,'position',p)   % set the new position    
             
                
                

                hax=axes();
                plot(zscore(rTC{i}(:,j))); axis([0 N -3 3]);
                newPos=[(1-1/nC)-(1/nC)*(fix((j-1)/1)+ihs-1),  vs*(mod(j-1,1)+ivs-shifter*(i-1)-0.13)+0.04,     1/shz,   svs-0.025];
                set(gca,'outer',newPos),
                set(gca,'XTickLabel','')
                set(gca,'YTickLabel','')
                xlabel(['|\rho|',' = ',num2str(round(TCcorr(i,j),2))],'color','k')
                xh = get(gca,'xlabel'); % handle to the label object
                q = get(xh,'position'); % get the current position property
                q(2) = 2+ q(2);       % double the distance,
                set(xh,'position',q)   % set the new position   
                
            end
        end
    end
