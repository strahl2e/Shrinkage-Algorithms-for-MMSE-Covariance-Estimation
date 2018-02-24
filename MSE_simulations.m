%% Replicate the MSE simulations from the paper

%% Setup general
p=100;

S_i_j_diffs = zeros(p);
for i=1:p
    S_i_j_diffs(i, i:p) = abs(i - (i:p));
    S_i_j_diffs(i:p, i) = S_i_j_diffs(i, i:p);
end
figure
imagesc(S_i_j_diffs)
title('S_i_j_diffs')


%% Setup AR1 covariance matrices

S_AR1_r_0_1 = 0.1.^S_i_j_diffs
figure
imagesc(S_AR1_r_0_1)
title('S_AR1_r_0_1')
S_AR1_r_0_5 = 0.5.^S_i_j_diffs
figure
imagesc(S_AR1_r_0_5)
title('S_AR1_r_0_5')
S_AR1_r_0_9 = 0.9.^S_i_j_diffs
figure
imagesc(S_AR1_r_0_9)
title('S_AR1_r_0_9')

%% setup FBM covariance matrices
% Note: definition in paper is wrong! 
% See MATLAB code example here for a proper implementation https://se.mathworks.com/matlabcentral/fileexchange/19797-simulation-of-fractional-gaussian-noise--exact-
% See this paper for a correct definition: https://arxiv.org/pdf/1709.06115.pdf
FBM_matrix = @(h, S_i_j_diffs) 0.5*(abs(S_i_j_diffs + 1).^(2*h) - 2*(S_i_j_diffs.^(2*h)) + abs(S_i_j_diffs - 1).^(2*h));

S_FBM_h_0_6 = FBM_matrix(0.6, S_i_j_diffs);
figure
imagesc(S_FBM_h_0_6)
title('S_FBM_h_0_6')

S_FBM_h_0_7 = FBM_matrix(0.7, S_i_j_diffs);
figure
imagesc(S_FBM_h_0_7)
title('S_FBM_h_0_7')

S_FBM_h_0_8 = FBM_matrix(0.8, S_i_j_diffs);
figure
imagesc(S_FBM_h_0_8)
title('S_FBM_h_0_8')

%% Create samples for testing AR(1) covariance matrix

iters = 10
ns = 6:2:30;

avg_rho_LWs = zeros(length(ns),1);
avg_rho_RBLWs = zeros(length(ns),1);
avg_rho_OASs = zeros(length(ns),1);

avg_MSE_LWs = zeros(length(ns), 1);
avg_MSE_RBLWs = zeros(length(ns), 1);
avg_MSE_OASs = zeros(length(ns), 1);

for n = ns
    sum_rhos_LW = 0;
    sum_rhos_RBLW = 0;
    sum_rhos_OAS = 0;
    
    sum_MSEs_LW = 0;
    sum_MSEs_RBLW = 0;
    sum_MSEs_OAS = 0;
    
    for i = 1:iters
        % Create n samples
        X = mvnrnd(zeros(p,1), S_AR1_r_0_1, n);
        
        % Approximate covariance matrix using LW, RBLW, OAS and true oracle
        S_hat = cov(X);
        F_hat = (trace(S_hat) / p) * eye(p);
        
        sum_sq_fro_diffs = 0;
        for j=1:n
            sum_sq_fro_diffs = sum_sq_fro_diffs + norm(X(j,:)'*X(j,:) - S_hat, 'fro')^2;
        end
        rho_LW =  min(1, sum_sq_fro_diffs / (n^2*(trace(S_hat^2) - (trace(S_hat)^2/p))));
        rho_RBLW = min(1,(((n-2)/n) * trace(S_hat^2) + trace(S_hat)^2) / ((n+2)*(trace(S_hat^2) - (trace(S_hat)^2/p))));
        rho_OAS = min(1, (((1-2)/p)*trace(S_hat^2) + trace(S_hat)^2) / (((n+1-2)/p)*(trace(S_hat^2) - trace(S_hat)^2/p)));

        sum_rhos_LW = sum_rhos_LW + rho_LW;
        sum_rhos_RBLW = sum_rhos_RBLW + rho_RBLW;
        sum_rhos_OAS = sum_rhos_OAS + rho_OAS;
        
        S_hat_LW = (1-rho_LW)*S_hat + rho_LW * F_hat;
        S_hat_RBLW = (1-rho_RBLW)*S_hat + rho_RBLW * F_hat;
        S_hat_OAS = (1-rho_OAS)*S_hat + rho_OAS * F_hat;
        
        % Record MSE
        LW_Fro_sq_sum = 0;
        RBLW_Fro_sq_sum = 0;
        OAS_Fro_sq_sum = 0;
        for j=1:n
            LW_Fro_sq_sum = LW_Fro_sq_sum + norm(S_hat_LW - S_AR1_r_0_1, 'fro')^2;
            RBLW_Fro_sq_sum = RBLW_Fro_sq_sum + norm(S_hat_RBLW - S_AR1_r_0_1, 'fro')^2;
            OAS_Fro_sq_sum = OAS_Fro_sq_sum + norm(S_hat_OAS - S_AR1_r_0_1, 'fro')^2;
        end
        sum_MSEs_LW = LW_Fro_sq_sum / n;
        sum_MSEs_RBLW = RBLW_Fro_sq_sum / n;
        sum_MSEs_OAS = OAS_Fro_sq_sum / n;
    end
% Repeat 'iter' times
end

% delete(findall(0,'Type','figure'))