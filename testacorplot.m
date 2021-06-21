
%% read data
B = importdata('IAcT_Nruns1000_N13to22.txt');

% last row contains reference solution
[M,POW] = size(B);
RELERROR = zeros(POW,1);
N = 2.^([13:13+POW-1]);

% acor simulation data
A = B(1:M-1,:);
% acor reference value
c = B(M,1);


%% compute relative error
for k = 1:POW
    
    RELERROR(k) = sqrt( mean((A(:,k) - c).^2)) / mean(A(:,k));
    
end

%% regression and double log plot

logx = log(N);
logy = log(RELERROR);
COVM = cov(logx,logy)/var(logx);
SLPE = COVM(1,2);


figure(1)
loglog(N,RELERROR,'b*--') 
hold on
loglog(N,7.5.*N.^SLPE,'r')
hold off
legend1 = sprintf('relative error (M=%1.0f)',M-1)
legend2 = sprintf('regression slope %f', SLPE)
legend(legend1,legend2)