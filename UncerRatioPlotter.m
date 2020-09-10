%uncertainity ratio plotter

function[KAPPA,KAPPA_B,UR,URG,ITERELAPSED_BICG,TOL_1,TOL_2,TOL_3,TOL_4,UR_FOR_ERROR_BICG,URGforError_GMRES,TOLG_1,TOLG_2,TOLG_3,TOLG_4] = UncerRatioPlotter(n)

d=10; it = 1; m=1;
for j=1:1
    
    j
    
     
         [A,S] = crazyA(n,10^6);
         %A = MatGen(10^6,n);
% 
%            A = rand(n);
% % %         
%           [U S V] = svd(A);
% % %         S(1,1) = 10^6; S(n,n) = 1;
% %           S = diag(10.^linspace(0,6,n));
%           for i = 1:n-1
%               S(i,i) = 10^6;
%           end
%           S(n,n) = 1;
%           A = U*S*V';
%         
        
        
        INV_A = inv(A);
        NORM_A = norm(A);
        NORM_INVA = norm(INV_A);
        %C(it) = cond(A)
        %it=it+1;
      
       
       for iter = 1:10000
         if (mod(iter,100)==0)
             iter
         end
           
%         [A,S] = crazyA(n,10^2);
%        INV_A = inv(A);
%        NORM_A = norm(A);
%        NORM_INVA = norm(INV_A);
        
        b = randn(n,1);
       % b = ones(n,1)/1000;
        x = A\b;
        
        temp =  NORM_A*norm(x)/norm(b);
        temp2 = NORM_INVA*norm(b)/norm(x);
        
        if (temp > 10 || temp2 > 10) %check this when doing for gmres 
        %if (temp2 > 10^2)
        kappa(m) = temp; %assume kappa as kappaf
        kappab(m) = temp2;
        [val1,val2,val3,val4] = bicgtest(A,b,b,zeros(n,1),INV_A,d);
        ur(m) = val1;
        urforerrorbicg(m) = val2;
        iterelapsed_bicg(m) = val3;
        tol_1(m) = val4(1); tol_2(m) = val4(2); tol_3(m) = val4(3); tol_4(m) = val4(4); %tol_5(m) = val4(5); tol_6(m) = val4(6);
        [valg1, valg2, valg3] = gmres11(A,b,10^-7,zeros(n,1),x,d);
        urG(m) = valg1; urGforErrorgmres(m) = valg2; 
        tolg_1(m) = valg3(1); tolg_2(m) = valg3(2); tolg_3(m) = valg3(3); tolg_4(m) = valg3(4); %tolg_5(m) = valg3(5); tolg_6(m) = valg3(6);
       % Lur(m) = sqrt(2/pi) * (kappa(m)/(n-d)) * (SofSv(1,n-d-1,S,n,d,1) + SofSv(n-d,n,S,n,d,2));
       % Lur(m) = (sqrt(sum(sum(S)))/max(max(S)))*kappa(m); %its not even a bound its E(UR)
        m=m+1;
        
        end
        
       end
       
       
    
%        ur = mean(ur);
%        urG = mean(urG);
%        kappa = mean(kappa);
%        
    if (j==1)
   % LUR = Lur;
    UR = ur;
    URG = urG;
    URGforError_GMRES = urGforErrorgmres;
    KAPPA = kappa;
    KAPPA_B = kappab;
    ITERELAPSED_BICG = iterelapsed_bicg;
    UR_FOR_ERROR_BICG = urforerrorbicg;
    TOL_1 = tol_1; TOL_2 = tol_2; TOL_3 = tol_3; TOL_4 = tol_4; %TOL_5 = tol_5; TOL_6 = tol_6;
    TOLG_1 = tolg_1; TOLG_2 = tolg_2; TOLG_3 = tolg_3; TOLG_4 = tolg_4; %TOLG_5 = tolg_5; TOLG_6 = tolg_6;
    else 
     %   LUR = [LUR,Lur];
        UR = [UR,ur];
        URG = [URG,urG];
        URGforError_GMRES = [URGforError_GMRES,urGforErrorgmres];
        KAPPA = [KAPPA,kappa];
        KAPPA_B = [KAPPA_B,kappab];
        ITERELAPSED_BICG = [ITERELAPSED_BICG,iterelapsed_bicg];
        UR_FOR_ERROR_BICG = [UR_FOR_ERROR_BICG,urforerrorbicg];
        TOL_1 = [TOL_1,tol_1]; TOL_2 = [TOL_2,tol_2]; TOL_3 = [TOL_3,tol_3]; TOL_4 = [TOL_4,tol_4]; %TOL_5 = [TOL_5,tol_5]; TOL_6 = [TOL_6,tol_6];
        TOLG_1 = [TOLG_1,tolg_1]; TOLG_2 = [TOLG_2,tolg_2]; TOLG_3 = [TOLG_3,tolg_3]; TOLG_4 = [TOLG_4,tolg_4]; %TOLG_5 = [TOLG_5,tolg_5]; TOL_6 = [TOLG_6,tolg_6];      
        end

    
end

%loglog(KAPPA,UR,'*')
%hist(KAPPA)

end