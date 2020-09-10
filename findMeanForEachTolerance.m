function [ urbins,tolbins ] = findMeanForEachTolerance(UR,tol)

%Bin1 = [100,1000), Bin2 = [1000,10000), Bin3 = [10000,100000)
urBin1 = []; urBin2 = []; urBin3 = []; TOLBin1=[]; TOLBin2=[]; TOLBin3=[]; 
m1=1; m2=1; m3=1;
for j=1:size(UR,2)
    if (UR(j) >= 100 && UR(j) < 1000)
       TOLBin1(m1) = tol(j);
       urBin1(m1) = UR(j);
       m1 = m1 + 1;
    else if (UR(j) >=1000 && UR(j) < 10^4)
       TOLBin2(m2) = tol(j);
       urBin2(m2) = UR(j);  
       m2 = m2 + 1;
        else 
            TOLBin3(m3) = tol(j);
            urBin3(m3) = UR(j);
            m3 = m3 + 1;
        end
    end
end

urbins = [mean(urBin1),mean(urBin2),mean(urBin3)];
tolbins = [mean(TolBin1),mean(TolBin2),mean(TolBin3)];

            
end
