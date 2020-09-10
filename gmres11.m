function [ur1,ur2,IterationsAtEachToleranceGMRES] = gmres11(A,b,tol,x0,xt,d)

    %initialization
    m = length(b);
    n2b = norm(b);
    iter = 0;
    k=1;
    r0 = b-A*x0;
    beta = norm(r0);
    alpha = norm(xt);
    resvec(k) = beta;
    err(k) = norm(x0-xt);
    rel_res(k) = resvec(k)/n2b;
    rel_err(k) = err(k)/alpha;
    rel_err_est(k) = 0;
    x_curr = x0;
    STORE1=zeros(1,6); STORE2=zeros(1,6); 
    

    
    %initialize q1
    Q = zeros(m,m);
    Q(:,1) = r0/beta;
    U = zeros(m,m);
    R_arr = zeros(m,m);
    v_arr = zeros(1,m);
    H = zeros(m,m);
    I = eye(m,m);
    

    %while( (alpha - norm(x_curr)) > 10^-3 && iter < m)
        %mgs process for hessenberg reduction
 while( norm(r0)/n2b > 10^-6 && iter < m)
            
        v = A*Q(:,k);
        for j=1:k
            H(j,k) = Q(:,j)'*v;
            v = v-H(j,k)*Q(:,j);
        end

        n2v = norm(v);

        % minimize norm(H*y-beta*e1)

        %constructing upper triangular matrix from H
        if k==1
            R_arr(1,1) = H(1,1);
        else
            v_arr(k-1) = H(k,k-1)/R_arr(k-1,k-1);
            R_arr(1:k,k) = H(1:k,k);
            for j=2:k
                R_arr(j,k) = R_arr(j,k) - v_arr(j-1)*R_arr(j-1,k);
            end
        end

        %finding solutions
        z = computeHessenSol(R_arr(1:k,1:k),v_arr(1:k-1),beta*I(1:k,1));
        t1 = solveHessen(transpose(H(1:k,1:k)),I(1:k,k),'lower');
        t = computeHessenSol(R_arr(1:k,1:k),v_arr(1:k-1),t1);

        
        %computing kth iterate and residuals
        h_sq = n2v^2;
        delta = h_sq/(1+h_sq*t(k,1));
        u = delta*t;
        y = z - z(k,1)*u;
        x = x0 + Q(:,1:k)*y;
        res = norm(b-A*x);
        zetta(k) = z(k,1)/beta;
        U(1:k,k) = u;
        s(k) = norm(u)*zetta(k);

        
        if k > d
            Ht = H(k-d+1:k,k-d+1:k);
            a = solveHessen(Ht,I(1:d,1),'upper');
            w = H(1:k-d,k-d+1:k)*a;
            bt = computeHessenSol(R_arr(1:k-d,1:k-d),v_arr(1:k-d-1),w);
            gamma = H(k-d+1,k-d)*zetta(k-d)/(1-H(k-d+1,k-d)*bt(k-d,1));
            meu_est_sq = (gamma^2*(norm(a)^2) + norm(gamma*bt + zetta(k-d)*U(1:k-d,k-d))^2);
            rel_err_est(k-d) = beta*sqrt(abs(meu_est_sq - s(k)^2))/alpha;
        end
        
        
         if (rel_res(k) > 10^-1)
                STORE1(1) = k;
                else if (rel_res(k) > 10^-2)
                STORE1(2) = k;
                    else if (rel_res(k) > 10^-3)
                    STORE1(3) = k;
                        else if (rel_res(k) > 10^-4)
                           STORE1(4) = k;
                            else if (rel_res(k) > 10^-5)
                            STORE1(5) = k;
                                else if (rel_res(k) > 10^-6)
                                STORE1(6) = k;
                                    end
                                end
                            end
                        end
                    end
         end
                
    
 if (k>=d+1)
    if (rel_err_est(k-d)  > 10^-1)
        STORE2(1) = k-d;
    else if (rel_err_est(k-d) > 10^-2)
            STORE2(2) = k-d;
            else if (rel_err_est(k-d) > 10^-3)
                STORE2(3) = k-d;
                 else if (rel_err_est(k-d) > 10^-4)
                    STORE2(4) = k-d;
                     else if (rel_err_est(k-d) > 10^-5)
                         STORE2(5) = k-d;
                        else if (rel_err_est(k-d) > 10^-6)
                            STORE2(6) = k-d;
                            end
                         end
                     end
                end
         end
    end
end
            
          
        %updating values
        k=k+1;
        resvec(k) = res;
        err(k) = norm(x-xt);  
        iter = iter + 1;
        rel_err(k) = err(k)/alpha;
        rel_res(k) = resvec(k)/n2b;
        
        % write the convergence criteria here
        if k > d+1
            if (rel_err(k) < tol)
                break;
            end
        end
                
        %updating values
        H(k,k-1) = n2v;
        Q(:,k) = v/n2v;
    
        x_curr = x;
       
    end
    
    
    
    rel_err = rel_err(1:k-d-1);
    rel_err_est = rel_err_est(1:k-d-1);
    rel_res = rel_res(1:k-d-1);

    %finding ur1 and ur2
    
    Itol = 1;
    while(Itol < k-d && rel_err(Itol) > 1e-3)
        Itol = Itol + 1;
    end
    
    
    Imin = 1;
    while(Imin < k-d && rel_err(Imin) > min(rel_err))
        Imin = Imin + 1;
    end
    
    Iconv = min(Imin,Itol);
    
    rr = rel_res(1:Iconv);
    re = rel_err(1:Iconv);
    rest = rel_err_est(1:Iconv);

    ur1 =  mean(abs((rel_res-rel_err)./min(rel_err,rel_res)));%./(rest-re)));
    ur2 =  mean(abs((rel_err_est-rel_err)./min(rel_err_est,rel_err)));
   % ur2 =  mean(abs((rr.^2-re.^2)./(rest.^2-re.^2)));
    
    
    
    %finding Ir,Ie,Iest
%     y2 = max([min(rel_err) min(rel_res) min(rel_err_est)]);
%     Ir = 1;
%     while( Ir < k-d && rel_res(Ir) > max(y2,1e-6))
%         Ir = Ir + 1;
%     end
%     
%     Ie = 1;
%     while( Ie < k-d && rel_err(Ie) > max(y2,1e-6))
%         Ie = Ie + 1;
%     end
%     
%     Iest = 1;
%     while( Iest < k-d && rel_err_est(Iest) > max(y2,1e-6))
%         Iest = Iest + 1;
%     end
%     
    

% finding iteration saved

%     rel_err = smooth(rel_err,min(5,d));
%     rel_res = smooth(rel_res,min(5,d));
%     rel_err_est = smooth(rel_err_est,min(5,d));
%    
%     [I1,I2,I3,valid] = compute_avg_iter5(log10(rel_res),log10(rel_err),log10(rel_err_est));    
%     results = [ur1 ur2 Ir Ie Iest I1 I2 I3 Iconv valid];
       
    
    
% %     figure(2);
%     semilogy(rel_err,'black') %'-k');
%     hold on;
%    semilogy(rel_res,'g') %:k');
%    hold on;
%    semilogy(rel_err_est,'r') 
%    
    %plot(rel_err_est,'--k');
% %     plot(ed,'-k');
% %     plot(err.^2,'--k');
%      set(gca,'yscale','log');
%     xlabel('iteration');
%     ylabel('relative error');
%     title('Convergence plot');
%     legend('true relative error','relative error estimate','relative residual');
%     hold off;

%oh dont change this to 5 and 6; let it remain 5 and 5; its logic !
STORE1(5)=[];STORE1(5)=[];
STORE2(5)=[];STORE2(5)=[];

if(all(STORE1)==1 && all(STORE2)==1)
IterationsAtEachToleranceGMRES =  abs(STORE1-STORE2)./min(STORE1,STORE2);
else
    STORE1=zeros(1,6);
    STORE2=zeros(1,6);
    IterationsAtEachToleranceGMRES =  abs(STORE1-STORE2);
end


end