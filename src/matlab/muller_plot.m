clear F N popnow poptime Graph tsz
%counter i j k 1;
N=1*10^7;
T=600;
U=1*10^-6
S=0.025
%S=0.000001

figure
hold on
axis([1 T 0 N])
set(gca, 'ytick', [])
set(gca, 'fontsize', 16)
N=[N];
tot=sum(N);
tot2=sum(N);
F=[1];
popnow=N;
piptime=popnow;

rndstate = 1;
rndnstate = 1;

%rndstate = rand('state');
%rndnstate = randn('state');


rand('state', rndstate);
randn('state', rndnstate);

for t=1:T
    time(1,t) = t;
    t
    avfit = sum(N.*F)/tot;
    lambda = N.*F./avfit;
    fn(1,t) = avfit;
    N=poissrnd(lambda);
    sizeN=size(N);
    i=1;
    L=1;
    counter=0;
    
    %waiting time until next mutation
    k=int32(exprnd(1/U))+1;
    
    while k<tot2
        while sum(N(:,1:L))<k
            if N(1,L) == 0
                N(:,L) = [];
                F(:,L) = [];
                L=L-1;
            end
            L=L+1;
        end
        s = exprnd(S);
        sz = N(1,L);
        test = sum(N(1,1:(L)));
        j = double(k-sum(N(1,1:(L-1))));
        gens = size(N);
        N = [N(1,1:(L-1)), (j-1), 1, (sz-j), N(1,(L+1):gens(1,2))];
        F = [F(1,1:(L)), F(L)*(1+s), F(1,(L):gens(1,2))];
        k = k + int32(exprnd(1/U));
        L=1;
        tot2=sum(N);
    end
    
    Graph = cumsum(N);
    Graph = [Graph' Graph'];
    Graph = [0 0; Graph];
    sizeG = size(Graph);
    C = [F' F'; 0 0];
    x = 1:1:sizeG(1,1);
    x = x./x;
    x = [x'*t x'*(t+1)];
    pcolor(x, Graph, C);
    shading flat;
end

colorbar
set(colorbar, 'box', 'off');
set(colorbar, 'FontName', 'Arial');
colormap(gray)
a = get(gca, 'Clim');
a(1,1) = 1;
set(gca, 'Clim', a);