function [ilower_eig_s,iupper_eig_s,eigenvect_s,steps_s]=homotopy(s,steps_s,op,eigenval_0,eigenvect_0,ilower_eig_0,margin,para)

% Applies the homotopy method to the operator stored in op.
%
% The starting value is given in s (typically s=0), numerical eigenvalues 
% and eigenvvector for this value should be given in eigenval_0 and 
% eigenvect_0, and rigorous lower bounds for the eigenvalues should be in
% ilower_eig_0. 
%
% margin is a parameter (<1) which controls how the breakpoints for the 
% homotopy are selected. 
%
% para.tol_simple_eig is a small parameter which controls when, in the 
% validation procedure using Intlab's verifyeig, close eigenvalues should 
% be validated together as a cluster. 
% 
% para.show_s and para.show_enclosure should be equal to 0 or 1, and 
% control how much information about is displayed during the computations.

%Checking that the intput are valid, and that there is "some room" for a
%first homotopy step.
if not(eigenval_0(end-1)<ilower_eig_0(end))
    error('Something seems wrong with the starting point of the homotopy')
else
    ind_eig=length(eigenval_0);
    if para.show_s
        fprintf('\nind_eig = %d\n',ind_eig)
    end
    eigs_lost=0;
    while ind_eig>2 && abs((eigenval_0(ind_eig-1)-ilower_eig_0(ind_eig))/ilower_eig_0(ind_eig))<1.5*margin
            ind_eig=ind_eig-1;    
            eigs_lost=eigs_lost+1;
    end
    if eigs_lost>1 && para.show_s
        fprintf('\n%d eigenvalues losts\n',eigs_lost)
        fprintf('\nind_eig = %d\n',ind_eig)
    end
end

eigenval_s=eigenval_0(1:ind_eig);
eigenvect_s=eigenvect_0(:,1:ind_eig);
ilower_eig_s=ilower_eig_0(1:ind_eig);
correc=10^-4;%parameter used to slightly adapt the computation of the homotopy breakpoints (not crucial)
pert=10^-5;%finite difference step-size used in delta_s_predictor to approximate the derivate of eigenvalues with respect to s
delta_s_min=10^-6;
while s<1 && ind_eig>1 %%% HOMOTOPY
    eigenval_s=eigenval_s(1:ind_eig-1);
    eigenvect_s=op.Reset_for_BC*eigenvect_s(:,1:ind_eig-1);
    nu_min=min(eigenval_s);
    nu_min=min(1.1*nu_min,nu_min-10^2);
    
    %% Computation of the homotopy "breakpoint" s_{k+1}
    small_nb_eigs=min(ind_eig-1,max(10,floor(ind_eig/100)));
    [delta_s_pred,correc]=delta_s_predictor(s,(1-sign(ilower_eig_s(ind_eig))*margin)*ilower_eig_s(ind_eig),margin,op,small_nb_eigs,pert,delta_s_min,correc,para.show_s);%cheap approximate computation of delta_s = s_{k+1}-s_{k}     
    delta_s=delta_s_pred;
    if para.show_s
        fprintf('\nPredicted delta_s : %f\n',delta_s)
    end
    [eigenval_s,eigenvect_s,s]=homotopy_breakpoint_check(s,delta_s,nu_min,(1-sign(ilower_eig_s(ind_eig))*margin)*ilower_eig_s(ind_eig),margin,ind_eig-1,eigenval_s,eigenvect_s,op,delta_s_min,para.show_s);%Checking that the computed delta_s is appropriate
    
    steps_s=[steps_s s];
    if para.show_s
        fprintf('\nBreakpoint for the homotopy at s=%f\n',s)
    end
    
    eigenvect_s=op.Mat_BC*eigenvect_s;
    eigenvect_s=orthogonalize(eigenvect_s,op.Mat_scal);%Making sure that the numerical eigenvectors are close to orthogonal
    eigenvect_s=clean_data(eigenvect_s);
    
    ieigenvect_s=op.iMat_BC*intval(op.Reset_for_BC*eigenvect_s);
    ieigenvect_s=op.Pad*ieigenvect_s;    

    %% Rayleigh-Ritz method, upper bounds  
    is=intval(s);
    iHs=is*op.iH+(1-is)*op.iH0;
    iv=ieigenvect_s;
    iHv=iHs*iv;
        
    iA0=iv'*op.iMat_scal*iv;
    iA1=iHv'*op.iMat_scal*iv;
    
    [upper_vect_Hs,upper_eig_Hs]=eig(iA1.mid,iA0.mid);%Numerical upper bounds
    [upper_eig_Hs,ind_upper_eig_Hs]=sort(real(diag(upper_eig_Hs)));
    upper_vect_Hs=upper_vect_Hs(:,ind_upper_eig_Hs);
    all_eigs=1;
    iupper_eig_s=validate_eigenvalues(upper_eig_Hs,upper_vect_Hs,iA1,iA0,para.tol_simple_eig,all_eigs);
    iupper_eig_s=sort(sup(real(iupper_eig_s)),'ascend');%Rigorous upper bounds
    
    if iupper_eig_s(ind_eig-1)<ilower_eig_s(ind_eig)
        %% Lehmann Maehly method, lower bounds
        nu=(iupper_eig_s(ind_eig-1)+ilower_eig_s(ind_eig))/2;
        inu=intval(nu);
        
        iHnuv=iHv-inu*iv;   
        
        iB1=iA1-inu*iA0;
        iB2=iHnuv'*op.iMat_scal*iHnuv; 
        [lower_vect_Hs,lower_eig_Hs]=eig(iB2.mid,iB1.mid);%Numerical lower bounds
        [lower_eig_Hs,ind_lower_eig_Hs]=sort(real(diag(lower_eig_Hs)),'descend');
        lower_vect_Hs=lower_vect_Hs(:,ind_lower_eig_Hs);
        
        all_eigs=0;
        if s==1
            all_eigs=1;
        end
        ilower_eig_s=validate_eigenvalues(lower_eig_Hs,lower_vect_Hs,iB2,iB1,para.tol_simple_eig,all_eigs);        
        ilower_eig_s=real(ilower_eig_s);           
        mask=not(isnan(ilower_eig_s));
        if not(ilower_eig_s(mask)<0)
            error('Unable to guarantuee that the eigenvalues are negative. This could be due to some interval coefficients in the operator beeing too wide.')
        end
        ilower_eig_s=inu+flip(ilower_eig_s);
        ilower_eig_s=inf(ilower_eig_s);
        mask=not(isnan(ilower_eig_s));
        ilower_eig_s(mask)=sort(ilower_eig_s(mask),'ascend');%Rigorous lower bounds
        
        if para.show_enclosure %|| s==1
            fprintf('\nEnclosure for the eigenvalues\n')
            [ilower_eig_s eigenval_s iupper_eig_s]
        end
    else
        error('We cannot prove that the last eigenvalue is well separated from the next one. This probably means the rigorous bounds are not tight enough. You may try to increase K, or to increase the size of the domain [a,b].')
    end
    
    if s<1
        ind_eig=ind_eig-1;
        eigs_lost=1;
        while ind_eig>2 && abs((iupper_eig_s(ind_eig-1)-ilower_eig_s(ind_eig))/ilower_eig_s(ind_eig))<1.5*margin && not(isnan(ilower_eig_s(ind_eig-1)))
            ind_eig=ind_eig-1;    
            eigs_lost=eigs_lost+1;
        end
        if ind_eig>2 && abs((iupper_eig_s(ind_eig-1)-ilower_eig_s(ind_eig))/ilower_eig_s(ind_eig))<1.5*margin
            %If we end up here, it means we need to drop more eigenvalues
            %than those for which we had computed rigorous lower bounds, so
            %we compute rigorous lower bounds for all of them (of course we
            %probably only need a few, so this could be somewhat optimized)
            all_eigs=1;
            ilower_eig_s=validate_eigenvalues(lower_eig_Hs,lower_vect_Hs,iB2,iB1,para.tol_simple_eig,all_eigs);        
            ilower_eig_s=real(ilower_eig_s);            
            mask=not(isnan(ilower_eig_s));
            if not(ilower_eig_s(mask)<0)
                error('Unable to guarantuee that the eigenvalues are negative')
            end    
            ilower_eig_s=inu+flip(ilower_eig_s);
            ilower_eig_s=inf(ilower_eig_s);
            mask=not(isnan(ilower_eig_s));
            ilower_eig_s(mask)=sort(ilower_eig_s(mask),'ascend');%Rigorous lower bounds
            if para.show_enclosure
                fprintf('\nEnclosure for the eigenvalues (now with all the lower bounds)\n')
                [ilower_eig_s eigenval_s iupper_eig_s]
            end
        end
        while ind_eig>2 && abs((iupper_eig_s(ind_eig-1)-ilower_eig_s(ind_eig))/ilower_eig_s(ind_eig))<1.5*margin
            ind_eig=ind_eig-1;    
            eigs_lost=eigs_lost+1;
        end
        if para.show_s
            fprintf('\n')
            if eigs_lost>2
                fprintf('\n%d eigenvalues losts\n',eigs_lost)
            end
            fprintf('\nind_eig = %d\n',ind_eig)
        end
    end   
end

if para.show_s
    steps_s
end

if s<1
    steps_s
    s
    error('Could not reach the end of the homotopy')
end
