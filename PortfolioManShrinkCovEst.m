%% Portfolio management with covariance shrinkage estimators

%% 1. Read in stock data for large n, large p
StockData = hist_stock_data('01012017','01012018','PortManCovShrinkEstTickerSymbols.txt','frequency','d');

%% 2. Try wealth management for small n, increasing n, comparing use of
% sample covariance vs methods in paper

N_S = length(StockData); % number of assets
N_t = length(StockData(1).AdjClose); % number of closing values

Cs = zeros(N_t,N_S); % asset costs
nrs = zeros(N_t-1,N_S); %n net returns
mus = zeros(N_S); % means
sds = zeros(N_S); % s.d.s

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

p_ss = 2:3:50; % number of previous samples (n) to use for computing weights
N_p_s = length(p_ss);

start_t = 50;
t_offset = start_t-1;
gam = 20;
T = length(start_t:size(nrs,1));
sc_names = {'SC','PC','F','LW','FBLW','OAS'}
N_sc = length(sc_names);

ws_t = zeros(T,N_S,N_p_s,N_sc); % weights
wl_t = zeros(T,N_S,N_p_s,N_sc);

Ws_t = zeros(T+1,N_p_s,N_sc); % wealth
Wl_t = zeros(T+1,N_p_s,N_sc);

for p_s_i = 1:N_p_s
    fprintf(2,'Computing weights and wealth with p=%i \n',p_ss(p_s_i));
    %p_s = 50; %Number of previous samples
    Ws_t(1,p_s_i,sc) = 1; % Start with 1 unit (million euros) of wealth
    Wl_t(1,p_s_i,sc) = 1;

    for c_t = start_t:size(nrs,1)
        % Calculate relevant statistics of net returns
        current_mus = mean(nrs(max(1,c_t-p_ss(p_s_i)):c_t,:))';
        current_samp_cov = cov(nrs(max(1,c_t-p_ss(p_s_i)):c_t,:));
        S_hat = zeros(N_S,N_S);
        for j=max(1,c_t-p_ss(p_s_i)):c_t
            S_hat = S_hat + nrs(j,:)' * nrs(j,:);
        end
        S_hat = S_hat ./ length(max(1,c_t-p_ss(p_s_i)):c_t);
        F_hat = (trace(S_hat) / N_S) * eye(N_S);
        
        S_hats = zeros(N_S, N_S, 6);
        S_hats(:,:,1) = current_samp_cov;
        S_hats(:,:,2) = S_hat;
        S_hats(:,:,3) = F_hat;
        
        sum_sq_fro_diffs = 0;
        for j=max(1,c_t-p_ss(p_s_i)):c_t
            sum_sq_fro_diffs = sum_sq_fro_diffs + norm(nrs(j,:)'*nrs(j,:) - S_hat, 'fro')^2;
        end
        n = length(max(1,c_t-p_ss(p_s_i)):c_t);
        rho_LW =  min(1, sum_sq_fro_diffs / (n^2*(trace(S_hat^2) - (trace(S_hat)^2/N_S))));
        rho_RBLW = min(1,(((n-2)/n) * trace(S_hat^2) + trace(S_hat)^2) / ((n+2)*(trace(S_hat^2) - (trace(S_hat)^2/N_S))));

        OAS_numer = (1-2/N_S)*trace(S_hat^2) + trace(S_hat)^2;
        OAS_denom = (n+1-2/N_S)*(trace(S_hat^2) - ((trace(S_hat)^2)/N_S));
        rho_OAS = min(1, OAS_numer / OAS_denom);

        S_hat_LW = (1-rho_LW)*S_hat + rho_LW * F_hat;
        S_hat_RBLW = (1-rho_RBLW)*S_hat + rho_RBLW * F_hat;
        S_hat_OAS = (1-rho_OAS)*S_hat + rho_OAS * F_hat;

        S_hats(:,:,4) = S_hat_LW;
        S_hats(:,:,5) = S_hat_RBLW;
        S_hats(:,:,6) = S_hat_OAS;
        
        % Compute weights using Expected Utility Maximization
        t = c_t - t_offset;
        for sc = 1:6
            samp_cov_inv = inv(S_hats(:,:,sc));
            one_tran_S = ones(N_S,1)'* samp_cov_inv;
            lambda_t = (gam - one_tran_S*current_mus) / (one_tran_S*ones(N_S,1));
            ws_t(t,:,p_s_i,sc) = (1/gam)*samp_cov_inv*(current_mus + lambda_t*ones(N_S,1));
            wl_t(t,:,p_s_i,sc) = max(min(1,ws_t(t,:,p_s_i,sc)), 0);

            % Compute wealth for next time step, units of one million.
            %mil_units = 1000000;
            ns_t = (ws_t(t,:,p_s_i,sc).*((Ws_t(t,p_s_i,sc)*1000000) * Cs(c_t,:).^-1))'; %Current short investments 
            nl_t = (wl_t(t,:,p_s_i,sc).*((Wl_t(t,p_s_i,sc)*1000000) * Cs(c_t,:).^-1))'; 

            Ws_t(t+1,p_s_i,sc) = (ns_t'*Cs(c_t+1,:)')/1000000; %Wealth investing in next time-step
            Wl_t(t+1,p_s_i,sc) = (nl_t'*Cs(c_t+1,:)')/1000000;
        end
    end
end

% Plot portfolio weights
figure
plot_num = 1;
for sc = 1:N_sc
    leg_txt = cell(length(p_ss)*2,1);
    lt_i = 1;
    for p_s_i = 1:N_p_s
        hold on
        subplot(3,2,plot_num)
        plot(ws_t(:,1,p_s_i,sc))
        leg_txt{lt_i} = strcat("S, p = ", num2str(p_ss(p_s_i)));
        lt_i = lt_i + 1;
    end
    for p_s_i = 1:N_p_s
        hold on
        subplot(3,2,plot_num)
        plot(wl_t(:,1,p_s_i,sc))
        leg_txt{lt_i} = strcat("L, p = ", num2str(p_ss(p_s_i)));
        lt_i = lt_i + 1;
    end
    hold off
    title(strcat('Portfolio weights',", (", sc_names{sc},")"))
    plot_num = plot_num + 1;
end
pw_leg = legend(leg_txt)


% Plot portfolio returns
Rs_t = zeros(length(ws_t)-1,N_p_s);
Rl_t = zeros(length(wl_t)-1,N_p_s);
for p_s_i = 1:N_p_s
    Rs_t(:,p_s_i) = sum(ws_t(1:end-1,:,p_s_i) .* nrs((start_t+1):end,:),2);
    Rl_t(:,p_s_i) = sum(wl_t(1:end-1,:,p_s_i) .* nrs((start_t+1):end,:),2);
end
mean(Rs_t)  
mean(Rl_t)
figure
leg_txt = cell(length(p_ss)*2,1);
lt_i = 1;
for p_s_i = 1:N_p_s
    hold on
    plot(Rs_t(:,p_s_i))
    leg_txt{lt_i} = strcat('S, p=',num2str(p_ss(p_s_i)));
    lt_i = lt_i +1;
end
for p_s_i = 1:N_p_s
    hold on
    plot(Rl_t(:,p_s_i))
    leg_txt{lt_i} = strcat('L, p=',num2str(p_ss(p_s_i)));
    lt_i = lt_i +1;
end
hold off
title('Returns')
legend(leg_txt)

[r c] = find(Rs_t == min(min(Rs_t)))

% Wealth process starting with w=1 at t=50 start.
fig=figure; 
leg_txt = cell(length(p_ss)*2,1);
lt_i = 1;
for p_s_i = 1:N_p_s
    hold on
    plot(Ws_t(:,p_s_i));
    leg_txt{lt_i} = strcat('S, p=',num2str(p_ss(p_s_i)));
    lt_i = lt_i +1;
end
for p_s_i = 1:N_p_s
    hold on
    plot(Wl_t(:,p_s_i));
    leg_txt{lt_i} = strcat('L, p=',num2str(p_ss(p_s_i)));
    lt_i = lt_i +1;
end
hold off
legend(leg_txt)

% delete(findall(0,'Type','figure'))