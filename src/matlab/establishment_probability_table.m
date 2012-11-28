% Jeffrey Barrick <jbarrick@msu.edu>
% Copyright (C) 2008-2009.

function out = establishment_probability_table(T,No,Step,Max,filename)

row = 1;
pr_gt_zero=0;
pr_gt_one=0;
for s=0.0:Step:Max
    if (s==0)
        Pe = 0;
    elseif (pr_gt_one)
        Pe = 1.0;
    else
        Ne = No*log(2)*T;
        Np = quad(@(x) establishment_probability(x, T, s, No), 0, T, 1.0e-6) / T;
        Pe = Np / Ne;
        
        if (Pe > 0)
            pr_gt_zero = 1;
        end;
            
        %it's possible for this to fail once s becomes too large
        %give it a chance of 1.0 if it fails for this reason
        if ((Pe==0) && (pr_gt_zero==1))
            Pe = 1.0;
            pr_gt_one=1;
        end;  
        
        if (Pe>=1.0)
            Pe=1.0;
            pr_gt_one=1;
        end;
    end
    
    s
    Pe
    M(row,:) = [s Pe];
    row=row+1;
end;
dlmwrite(filename, M,'\t');