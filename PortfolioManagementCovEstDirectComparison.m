%% Compare covariance approximation of a number of stocks
% Stocks chosen randomly from this list https://www.nasdaq.com/screening/companies-by-name.aspx?letter=O
%% Read in stocks
StockData = hist_stock_data('01012017','01012018','PortManCovDirCompTickerSymbols.txt','frequency','d');

N_S = length(StockData);
N_t = length(StockData(1).AdjClose);
Cs = zeros(N_t,N_S); % asset costs

for s = 1:N_S
    Cs(:,s) = StockData(s).AdjClose;
end

figure
for s = 1:N_S
    hold on
    plot(Cs(:,s))
end

figure
for s = 1:N_S
    hold on
    plot(Cs(:,s) / max(Cs(:,s)))
end

%% Use all samples as a 'true' covariance
S_hat = zeros(N_S, N_S);
for t = 1:N_t
    S_hat = S_hat + Cs(t,:)'*Cs(t,:);
end
S_hat_S = S_hat / (N_t +1);
S_hat_P = S_hat / N_t;

%% Now approximate the 'true' covariance  with subsamples and compare methods
p_ts = [2:2:50];
final_t = p_ts(end)+1;
for p_t = p_ts
    S_hat = zeros(N_S, N_S);
    for t = final_t-p_t:final_t
        S_hat = S_hat + Cs(t,:)'*Cs(t,:);
    end
end

