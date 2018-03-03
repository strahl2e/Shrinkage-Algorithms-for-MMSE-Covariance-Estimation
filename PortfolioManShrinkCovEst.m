%% Portfolio management with covariance shrinkage estimators

%% 1. Read in stock data for large n, large p
StockData = hist_stock_data('01012017','01012018','PortManCovShrinkEstTickerSymbols.txt','frequency','d');

%% 2. Try wealth management for small n, increasing n, comparing use of
% sample covariance vs methods in paper

N_S = length(StockData);
N_t = length(StockData(1).AdjClose);
nrs = zeros(N_t-1,N_S);
mus = zeros(N_S);
sds = zeros(N_S);

for s = 1:N_S
    t = 2:length(StockData(s).AdjClose);
    t_m_1 = t-1;

    % Daily net return r_{1t} and r_{2t}, r_{pt} =  ( C_{i,t} - C_{i,t-1} ) / C_{i,t-1}
    % Daily log return = log(1 + r_{pt})
    C = StockData(s).AdjClose;
    nr = (C(t) ./ C(t_m_1)) - 1;
    nrs(:,s) = nr;
    %lr_SP = log(C(t)) - log(C(t_m_1)); %For computational efficiency use log - log
    % Standardised returns: r_tilde_{it} = (r_{it} - mu_tilde_i)/sigma_tilde_i
    
    mus(s) = mean(nr);
    sds(s) = sqrt(var(nr));
end

p_s = 50; %Number of previous samples
start_t = 50;
t_offset = start_t-1;
gam = 20;
T = length(start_t:size(nrs,1));
ws_t = zeros(T,1);
wl_t = zeros(T,1);

for c_t = start_t:size(nrs,1)
	% Calculate relevant means
    current_mus = mean(nrs(max(1,c_t-p_s):c_t,:))';
    current_samp_cov = cov(nrs(max(1,c_t-p_s):c_t,:));
    
    % Compute weights using Expected Utility Maximization
    t = c_t - t_offset;
    samp_cov_inv = inv(current_samp_cov);
    one_tran_S = ones(N_S,1)'* samp_cov_inv;
    lambda_t = (gam - one_tran_S*current_mus) / (one_tran_S*ones(N_S,1));
    ws_t(t) = (1/gam)*samp_cov_inv*(current_mus + lambda_t*ones(N_S,1));
    wl_t(t) = max(min(1,ws_t(t)), 0);

    
end