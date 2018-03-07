%% Compare covariance approximation of a number of stocks
% Stocks chosen randomly from this list https://www.nasdaq.com/screening/companies-by-name.aspx?letter=O
%% Read in stocks
StockData = hist_stock_data('01012017','01012018','PortManCovDirCompTickerSymbols.txt','frequency','d');

N_S = length(StockData);
N_t = length(StockData(1).AdjClose);
Cs = zeros(N_t,N_S); % asset costs
N_nr = N_t-1;
nrs = zeros(N_nr,N_S);

for s = 1:N_S
    Cs(:,s) = StockData(s).AdjClose;
    t = 2:length(Cs(:,s));
    t_m_1 = t-1;

    % Daily net return r_{1t} and r_{2t}, r_{pt} =  ( C_{i,t} - C_{i,t-1} ) / C_{i,t-1}
    % Daily log return = log(1 + r_{pt})
    
    nr = (Cs(t,s) ./ Cs(t_m_1,s)) - 1;
    nrs(:,s) = nr;
    %lr_SP = log(C(t)) - log(C(t_m_1)); %For computational efficiency use log - log
    % Standardised returns: r_tilde_{it} = (r_{it} - mu_tilde_i)/sigma_tilde_i
end

figure
for s = 1:N_S
    hold on
    plot(nrs(:,s))
end
figure
for s = 1:N_S
    hold on
    plot(Cs(:,s) / max(Cs(:,s)))
end

%% Use all samples as a 'true' covariance
S_true = zeros(N_S, N_S);
for t = 1:N_nr
    S_true = S_true + nrs(t,:)'*nrs(t,:);
end
S_S = S_true / (N_nr +1);
S_P = S_true / N_nr;  % We use this for the MSE as this is used in the paper.
S_M_S = cov(nrs);
S_M_P = cov(nrs,1);

imagesc(S_S);
imagesc(S_P);
colorbar
title("Covariance (1/N) for all 251 daily closing prices")
imagesc(S_M_S);
imagesc(S_M_P);

%% Now approximate the 'true' covariance  with subsamples and compare methods
p_ts = [2:2:50];
N_ps = length(p_ts);
S_hat_names = {'S (1/(N+1))}', 'S (1/N)', 'S M(1/(N+1))', 'S M(1/N)', 'F hat', 'LW', 'RBLW', 'OAS'}
N_Sh = length(S_hat_names);

final_ts = (p_ts(end)+1):N_t;
N_fts = length(final_ts);
MSEs = zeros(N_ps, N_Sh, N_fts);
rhos = zeros(N_ps, N_Sh, N_fts);

for final_t_i = 1:N_fts
    final_t = final_ts(final_t_i);
    S_hat_approx_err = zeros(N_S, N_S, N_ps, N_Sh);
    for p_t_i = 1:N_ps
        p_t = p_ts(p_t_i);
        S_hat = zeros(N_S, N_S);
        start_sample = final_t-p_t;
        end_sample = final_t-1;
        for t = start_sample:end_sample
            S_hat = S_hat + nrs(t,:)'*nrs(t,:);
        end
        S_S_hat = S_hat / (N_t +1);
        S_P_hat = S_hat / N_t;
        S_M_S_hat = cov(nrs(start_sample:end_sample,:));
        S_M_P_hat = cov(nrs(start_sample:end_sample,:),1);

        F_hat = (trace(S_P_hat) / N_S) * eye(N_S);

        S_hats = zeros(N_S, N_S, N_Sh);
        S_hats(:,:,1) = S_S_hat;
        S_hats(:,:,2) = S_P_hat;
        S_hats(:,:,3) = S_M_S_hat;
        S_hats(:,:,4) = S_M_P_hat;
        S_hats(:,:,5) = F_hat;

        sum_sq_fro_diffs = 0;
        for j=start_sample:end_sample
            sum_sq_fro_diffs = sum_sq_fro_diffs + norm(nrs(j,:)'*nrs(j,:) - S_P_hat, 'fro')^2;
        end
        n = (end_sample - start_sample) + 1;
        rho_LW =  min(1, sum_sq_fro_diffs / (n^2*(trace(S_P_hat^2) - (trace(S_P_hat)^2/N_S))));
        rho_RBLW = min(1,(((n-2)/n) * trace(S_P_hat^2) + trace(S_P_hat)^2) / ((n+2)*(trace(S_P_hat^2) - (trace(S_P_hat)^2/N_S))));

        OAS_numer = (1-2/N_S)*trace(S_P_hat^2) + trace(S_P_hat)^2;
        OAS_denom = (n+1-2/N_S)*(trace(S_P_hat^2) - ((trace(S_P_hat)^2)/N_S));
        rho_OAS = min(1, OAS_numer / OAS_denom);

        S_hat_LW = (1-rho_LW)*S_P_hat + rho_LW * F_hat;
        S_hat_RBLW = (1-rho_RBLW)*S_P_hat + rho_RBLW * F_hat;
        S_hat_OAS = (1-rho_OAS)*S_P_hat + rho_OAS * F_hat;

        S_hats(:,:,6) = S_hat_LW;
        S_hats(:,:,7) = S_hat_RBLW;
        S_hats(:,:,8) = S_hat_OAS;

        rhos(p_t_i, 1, final_t_i) = 0;
        rhos(p_t_i, 2, final_t_i) = 0;
        rhos(p_t_i, 3, final_t_i) = 0;
        rhos(p_t_i, 4, final_t_i) = 0;
        rhos(p_t_i, 5, final_t_i) = 1;
        rhos(p_t_i, 6, final_t_i) = rho_LW;
        rhos(p_t_i, 7, final_t_i) = rho_RBLW;
        rhos(p_t_i, 8, final_t_i) = rho_OAS;
        
        % Check MSE 
        for sh = 1:N_Sh
            
            MSEs(p_t_i, sh, final_t_i) = norm(S_hats(:,:,sh) - S_P, 'fro')^2;
        end
    end
end

MSE_avgs = sum(MSEs,3) / N_fts;
rho_avgs = sum(rhos,3) / N_fts;

figure
for sh = 1:N_Sh
    hold on
    plot(MSE_avgs(:,sh));
end
legend(S_hat_names)
title("Avg MSEs for increasing numbers of subsamples")
xlabel("p (#previous time-step samples)")
ylabel("MSE (Frobenius norm)")

figure
for sh = [1 2 5 6 7 8]
    hold on
    plot(MSE_avgs(:,sh));
end
legend(S_hat_names([1 2 5 6 7 8]))
title("Avg MSEs for increasing numbers of subsamples")
xlabel("p (#previous time-step samples)")
ylabel("MSE (Frobenius norm)")
set(gca,'XTick',0:5:25)
set(gca,'XTickLabel',0:10:50)

figure
p_range = 1:3;
for sh = [1 2 5 6 7 8]
    hold on
    plot(MSE_avgs(p_range,sh));
end
legend(S_hat_names([1 2 5 6 7 8]))
title("Avg MSEs for increasing numbers of subsamples")
xlabel("p (#previous time-step samples)")
ylabel("MSE (Frobenius norm)")
set(gca,'XTick',1:3)
set(gca,'XTickLabel',2:2:6)

figure
for sh = [1 2 5 6 7 8]
    hold on
    plot(rho_avgs(p_range,sh));
end
legend(S_hat_names([1 2 5 6 7 8]))
title("Avg \rho values for increasing numbers of subsamples")
xlabel("p (#previous time-step samples)")
ylabel("\rho (0= S hat, 1 = F hat)")
set(gca,'XTick',1:3)
set(gca,'XTickLabel',2:2:6)
    