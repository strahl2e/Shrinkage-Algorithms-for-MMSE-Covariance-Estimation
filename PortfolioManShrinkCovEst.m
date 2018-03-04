%% Portfolio management with covariance shrinkage estimators

%% 1. Read in stock data for large n, large p
StockData = hist_stock_data('01012017','01012018','PortManCovShrinkEstTickerSymbols.txt','frequency','d');

%% 2. Try wealth management for small n, increasing n, comparing use of
% sample covariance vs methods in paper

N_S = length(StockData);
N_t = length(StockData(1).AdjClose);

Cs = zeros(N_t,N_S);
nrs = zeros(N_t-1,N_S);
mus = zeros(N_S);
sds = zeros(N_S);

for s = 1:N_S
    t = 2:length(StockData(s).AdjClose);
    t_m_1 = t-1;

    % Daily net return r_{1t} and r_{2t}, r_{pt} =  ( C_{i,t} - C_{i,t-1} ) / C_{i,t-1}
    % Daily log return = log(1 + r_{pt})
    Cs(:,s) = StockData(s).AdjClose;
    nr = (Cs(t,s) ./ Cs(t_m_1,s)) - 1;
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
ws_t = zeros(T,N_S); % weights
wl_t = zeros(T,N_S);

Ws_t = zeros(T+1,1); % wealth
Wl_t = zeros(T+1,1);
Ws_t(1) = 1; % Start with 1 unit (million euros) of wealth
Wl_t(1) = 1;

for c_t = start_t:size(nrs,1)
	% Calculate relevant statistics of net returns
    current_mus = mean(nrs(max(1,c_t-p_s):c_t,:))';
    current_samp_cov = cov(nrs(max(1,c_t-p_s):c_t,:));
    
    % Compute weights using Expected Utility Maximization
    t = c_t - t_offset;
    samp_cov_inv = inv(current_samp_cov);
    one_tran_S = ones(N_S,1)'* samp_cov_inv;
    lambda_t = (gam - one_tran_S*current_mus) / (one_tran_S*ones(N_S,1));
    ws_t(t,:) = (1/gam)*samp_cov_inv*(current_mus + lambda_t*ones(N_S,1));
    wl_t(t,:) = max(min(1,ws_t(t,:)), 0);
    
    % Compute wealth for next time step, units of one million.
    %mil_units = 1000000;
    ns_t = (ws_t(t,:).*((Ws_t(t)*1000000) * Cs(c_t,:).^-1))'; %Current short investments 
    nl_t = (wl_t(t,:).*((Wl_t(t)*1000000) * Cs(c_t,:).^-1))'; 
    
    
    Ws_t(t+1) = (ns_t'*Cs(c_t+1,:)')/1000000; %Wealth investing in next time-step
    Wl_t(t+1) = (nl_t'*Cs(c_t+1,:)')/1000000;
end

% Plot portfolio weights
figure
plot(ws_t(:,1))
hold on
plot(wl_t(:,1))
hold off

% Plot portfolio returns
Rs_t = sum(ws_t(1:end-1,:) .* nrs((start_t+1):end,:),2);
Rl_t = sum(wl_t(1:end-1,:) .* nrs((start_t+1):end,:),2);

mean(Rs_t)
mean(Rl_t)

plot(Rs_t,'k')
hold on
plot(Rl_t,'r')
hold off

% Wealth process starting with w=1 at t=50 start.
fig=figure; 
hax=axes; 
plot(Ws_t,'k');
hold on
plot(Wl_t,'r');
hold on
plot(W_SP_t,'b');
hold on
plot(W_ND_t,'g');
hold on
line([R_t_min_idx+1 R_t_min_idx+1], get(hax,'YLim'), 'Color', [0 1 0], 'LineStyle',':')
hold off

% delete(findall(0,'Type','figure'))