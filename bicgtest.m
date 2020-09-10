
%OnePoint is real UR ; OnePoint2 depicts the UR for error
function [OnePoint,OnePoint2,IterElapsed,IterationsAtEachToleranceBiCG] = bicgtest(A,b1,c1,x0,INV_A,d)

dim=size(A,1); TR_A = A';
i=1; STORE1=zeros(1,6); STORE2=zeros(1,6); 

r0=b1-A*x0;
p0=r0;
s0=c1;
q0=s0;

R=[];P=[];X=[]; 

x_final=A\b1;
f = sqrt(abs(x_final'*A*x_final));
g = sqrt(abs(b1'*A*b1));
Norm_b = norm(b1);

y0=zeros(dim(1),1);
%x0=randn(dim(1),1);     


while (  norm(r0)/Norm_b > (10^-6)  && i<dim(1)  )    
    
    const1 = A*p0;
    alpha(i)=s0'*r0/(q0'*const1);
    
    x1=x0 + alpha(i)*p0;
    y1=y0 + alpha(i)*q0;

    r1=r0 - alpha(i)*const1;
    s1=s0 - alpha(i)*TR_A*q0;
    
    n1= s1'*r1/(s0'*r0);
    
    p1= r1 + n1*p0;
    q1= s1 + n1*q0;
    
    if (i<=d+1)
    R=[R r0];
    P=[P p0];
    X=[X x0];
    end
    
    if (i>=d+2)
        R = R(:,2:d+1);
        P = P(:,2:d+1);
        X = X(:,2:d+1);
        R=[R r0];
        P=[P p0];
        X=[X x0];
    end
    
    %error(i) =  norm(x_final - x0);
    A_norm(i) = sqrt(abs(r1'*INV_A*r1))/abs(f);
    l2_norm(i) = norm(x_final-x0)/norm(x_final);
    rel_residual(i) = norm(r0)/Norm_b;
    rel_residual2(i) = sqrt(abs(r1'*A*r1))/g;
    
    
    S(:,i) = alpha(i)*p0;
    
    if (i>=d+1)
    error_Anorm(i-d) = sqrt(abs(-alpha(i-d)*R(:,1)'*P(:,1) + (alpha(i-d)^2)*P(:,1)'*A*P(:,1) + R(:,2)'*sum(S(:,i-d:i)')'))/abs(f);
    error_l2_norm(i-d) = sqrt(abs( - alpha(i-d)*P(:,1)'*(sum(S(:,i-d:i)')') + (sum(S(:,i-d+1:i)')')'*(sum(S(:,i-d:i)')') + alpha(i-d)^2 * norm(P(:,1))^2 ))/norm(x_final);    
    %error_l2_norm(i-d) = sqrt(abs( norm(sum(S(:,i-d:i)')')^2 + (alpha(i-d)^2)*(norm(P(:,1))^2) - 2*alpha(i-d)*P(:,1)'*sum(S(:,i-d:i)')') )/norm(x_final);
    %avg_denl(i-d) = abs(1 - norm(R(:,1))*norm(x_final)/(norm(b1)*l2_norm(i-d)*norm(x_final)) ); % /  abs(1 - norm(error_l2_norm(i-d)/l2_norm(i-d)));
    avg_denl(i-d) = abs((norm(R(:,1))/norm(b1) - l2_norm(i-d)/norm(x_final))) /  min((l2_norm(i-d)/norm(x_final)),norm(R(:,1))/norm(b1)); %this is new UR
   % avg_denl(i-d) = abs(1 - forward_prob_condl(A,X(:,1)) )  /  abs(1 - norm(error_l2_norm(i-d)/l2_norm(i-d)));
    avg_error(i-d) = abs(error_l2_norm(i-d) - l2_norm(i-d)) / min(l2_norm(i-d),error_l2_norm(i-d)); %this is new UR for error
    end
    
    
    if (rel_residual(i) > 10^-1)
        STORE1(1) = i;
    else if (rel_residual(i) > 10^-2)
            STORE1(2) = i;
            else if (rel_residual(i) > 10^-3)
            STORE1(3) = i;
                 else if (rel_residual(i) > 10^-4)
            STORE1(4) = i;
                     else if (rel_residual(i) > 10^-5)
            STORE1(5) = i;
              else if (rel_residual(i) > 10^-6)
            STORE1(6) = i;
                     end
                end
        end
        end
        end
    end
            

if (i>=d+1)
 if (error_l2_norm(i-d)  > 10^-1)
        STORE2(1) = i-d;
    else if (error_l2_norm(i-d) > 10^-2)
            STORE2(2) = i-d;
            else if (error_l2_norm(i-d) > 10^-3)
            STORE2(3) = i-d;
                 else if (error_l2_norm(i-d) > 10^-4)
            STORE2(4) = i-d;
                     else if (error_l2_norm(i-d) > 10^-5)
            STORE2(5) = i-d;
                else if (error_l2_norm(i-d) > 10^-6)
            STORE2(6) = i-d;
                     end
                end
        end
 end
 end
 end
end
      
    p0=p1;
    q0=q1;

    r0=r1;
    s0=s1;
    
    x0=x1;
    y0=y1;
    
    i=i+1;
end

% error_Anorm(i-d:i-1)=[];
% rel_residual2(i-d:i-1)=[];
% l2_norm(i-d:i-1)=[];
% A_norm(i-d:i-1)=[];
% error_l2_norm(i-d:i-1-d)=[];

% semilogy(A_norm);
% hold on 
% semilogy(error_Anorm,'r');
% hold on
% semilogy(abs(rel_residual2),'g');
% figure(2)
% semilogy(abs(rel_residual),'g');
% hold on
% semilogy(abs(l2_norm),'black')
% hold on
% semilogy(abs(error_l2_norm),'m');

%uncertainty in residue 
OnePoint  = mean(avg_denl);
%uncertainty in error estimator 
OnePoint2 = mean(avg_error);

%oh dont change this to 5 and 6; let it remain 5 and 5; its logic !
STORE1(5)=[];STORE1(5)=[];
STORE2(5)=[];STORE2(5)=[];

%chosing relative difference of iterations saved rather than just
%iterations saved
IterElapsed = mean(abs(STORE1-STORE2));
if(all(STORE1)==1 && all(STORE2)==1)
IterationsAtEachToleranceBiCG = abs(STORE1-STORE2)./min(STORE1,STORE2);
else
    STORE1=zeros(1,6);
    STORE2=zeros(1,6);
    IterationsAtEachToleranceBiCG = abs(STORE1-STORE2);
end


end