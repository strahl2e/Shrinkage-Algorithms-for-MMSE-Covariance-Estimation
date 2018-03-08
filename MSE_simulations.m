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
title('S AR1 r=0.1')
colorbar
S_AR1_r_0_5 = 0.5.^S_i_j_diffs
figure
imagesc(S_AR1_r_0_5)
title('S_AR1_r_0_5')
S_AR1_r_0_9 = 0.9.^S_i_j_diffs
figure
imagesc(S_AR1_r_0_9)
title('S AR1 r=0.9')
colorbar

%% setup FBM covariance matrices
% Note: definition in paper is wrong! 
% See MATLAB code example here for a proper implementation https://se.mathworks.com/matlabcentral/fileexchange/19797-simulation-of-fractional-gaussian-noise--exact-
% See this paper for a correct definition: https://arxiv.org/pdf/1709.06115.pdf
FBM_matrix = @(h, S_i_j_diffs) 0.5*(abs(S_i_j_diffs + 1).^(2*h) - 2*(S_i_j_diffs.^(2*h)) + abs(S_i_j_diffs - 1).^(2*h));

S_FBM_h_0_6 = FBM_matrix(0.6, S_i_j_diffs);
figure
imagesc(S_FBM_h_0_6)
title('S FBM h=0.6')
colorbar

S_FBM_h_0_7 = FBM_matrix(0.7, S_i_j_diffs);
figure
imagesc(S_FBM_h_0_7)
title('S_FBM_h_0_7')

S_FBM_h_0_8 = FBM_matrix(0.8, S_i_j_diffs);
figure
imagesc(S_FBM_h_0_8)
title('S FBM h=0.8')
colorbar

%% Create samples for testing AR(1) covariance matrix

tst_cov_mats = {S_AR1_r_0_1, S_AR1_r_0_5, S_AR1_r_0_9, S_FBM_h_0_6, S_FBM_h_0_7, S_FBM_h_0_8}
tst_cov_mat_names = {'S_{AR1} r = 0.1', 'S_{AR1} r=0.5', 'S_{AR1} r = 0.9', 'S_{FBM} h = 0.6', 'S_{FBM} h = 0.7', 'S_{FBM} h = 0.8'}
tst_cov_mat_filenames = {'S_AR1_r_0_1', 'S_AR1_r_0_5', 'S_AR1_r_0_9', 'S_FBM_h_0_6', 'S_FBM_h_0_7', 'S_FBM_h_0_8'}
iters = 5000;
ns = 6:2:30;

avg_rho_LWs = zeros(length(tst_cov_mats), length(ns),1);
avg_rho_RBLWs = zeros(length(tst_cov_mats), length(ns),1);
avg_rho_OASs = zeros(length(tst_cov_mats), length(ns),1);
avg_rho_Ora = zeros(length(tst_cov_mats), length(ns),1);

avg_MSE_LWs = zeros(length(tst_cov_mats), length(ns), 1);
avg_MSE_RBLWs = zeros(length(tst_cov_mats), length(ns), 1);
avg_MSE_OASs = zeros(length(tst_cov_mats), length(ns), 1);
avg_MSE_Ora = zeros(length(tst_cov_mats), length(ns), 1);

for tst_mat_i = 1:length(tst_cov_mats)
    for n = ns
        fprintf(2,'Evaluating covariance matrix %s with n=%d \n',tst_cov_mat_names{tst_mat_i},n);
        sum_rhos_LW = 0;
        sum_rhos_RBLW = 0;
        sum_rhos_OAS = 0;
        sum_rhos_Ora = 0;

        sum_MSEs_LW = 0;
        sum_MSEs_RBLW = 0;
        sum_MSEs_OAS = 0;
        sum_MSEs_Ora = 0;

        for i = 1:iters
            if mod(i,500) == 0
                fprintf(2,'Iteration %d of %d \n', i, iters);
            end
            % Create n samples using chosen covariance matrix
            TrueCov = tst_cov_mats{tst_mat_i};
            X = mvnrnd(zeros(p,1), TrueCov, n);

            % Approximate covariance matrix using LW, RBLW, OAS and true oracle
            % S_hat = cov(X,1); MATLAB removes column means!
            S_hat = zeros(p,p);
            for j=1:n
                S_hat = S_hat + X(j,:)' * X(j,:);
            end
            S_hat = S_hat ./ n;
            
            F_hat = (trace(S_hat) / p) * eye(p);

            sum_sq_fro_diffs = 0;
            for j=1:n
                sum_sq_fro_diffs = sum_sq_fro_diffs + norm(X(j,:)'*X(j,:) - S_hat, 'fro')^2;
            end
            rho_LW =  min(1, sum_sq_fro_diffs / (n^2*(trace(S_hat^2) - (trace(S_hat)^2/p))));
            rho_RBLW = min(1,(((n-2)/n) * trace(S_hat^2) + trace(S_hat)^2) / ((n+2)*(trace(S_hat^2) - (trace(S_hat)^2/p))));
            
            % Equation 23 is incorrect! Check equation 25 for the correct
            % formulation!
            OAS_numer = (1-2/p)*trace(S_hat^2) + trace(S_hat)^2;
            OAS_denom = (n+1-2/p)*(trace(S_hat^2) - ((trace(S_hat)^2)/p));
            rho_OAS = min(1, OAS_numer / OAS_denom);
            
            Ora_numer = (1 - 2/p)*trace(TrueCov^2) + trace(TrueCov)^2;
            Ora_denom = (n + 1 - 2/p)*trace(TrueCov^2) + (1 - n/p)*trace(TrueCov)^2;
            rho_Ora = Ora_numer / Ora_denom;
            
            S_hat_LW = (1-rho_LW)*S_hat + rho_LW * F_hat;
            S_hat_RBLW = (1-rho_RBLW)*S_hat + rho_RBLW * F_hat;
            S_hat_OAS = (1-rho_OAS)*S_hat + rho_OAS * F_hat;
            S_hat_Ora = (1-rho_Ora)*S_hat + rho_Ora * F_hat;

            % Record MSEs
            sum_MSEs_LW = sum_MSEs_LW + norm(S_hat_LW - TrueCov, 'fro')^2;
            sum_MSEs_RBLW = sum_MSEs_RBLW + norm(S_hat_RBLW - TrueCov, 'fro')^2;
            sum_MSEs_OAS = sum_MSEs_OAS + norm(S_hat_OAS - TrueCov, 'fro')^2;
            sum_MSEs_Ora = sum_MSEs_Ora + norm(S_hat_Ora - TrueCov, 'fro')^2;

            % Track rho values
            sum_rhos_LW = sum_rhos_LW + rho_LW;
            sum_rhos_RBLW = sum_rhos_RBLW + rho_RBLW;
            sum_rhos_OAS = sum_rhos_OAS + rho_OAS;
            sum_rhos_Ora = sum_rhos_Ora + rho_Ora;
        end % Repeat 'iter' times
        % Record average rho values and MSEs
        i_n = find(ns == n);
        avg_rho_LWs(tst_mat_i,i_n,:) = sum_rhos_LW / iters;
        avg_rho_RBLWs(tst_mat_i,i_n,:) = sum_rhos_RBLW / iters;
        avg_rho_OASs(tst_mat_i,i_n,:) = sum_rhos_OAS / iters;
        avg_rho_Ora(tst_mat_i,i_n,:) = sum_rhos_Ora / iters;

        avg_MSE_LWs(tst_mat_i,i_n,:) = sum_MSEs_LW / iters;
        avg_MSE_RBLWs(tst_mat_i,i_n,:) = sum_MSEs_RBLW / iters;
        avg_MSE_OASs(tst_mat_i,i_n,:) = sum_MSEs_OAS / iters;
        avg_MSE_Ora(tst_mat_i,i_n,:) = sum_MSEs_Ora / iters;
    end
    figure
    hold on
    plot(avg_rho_LWs(tst_mat_i,:,:))
    plot(avg_rho_RBLWs(tst_mat_i,:,:))
    plot(avg_rho_OASs(tst_mat_i,:,:))
    plot(avg_rho_Ora(tst_mat_i,:,:))
    title(['Covariance ', tst_cov_mat_names{tst_mat_i} ])
    legend('LW', 'RBLW', 'OAS', 'Oracle');
    hold off
end

for tst_mat_i = 1:length(tst_cov_mats)
    figure
    hold on
    plot(avg_rho_LWs(tst_mat_i,:,:),'-o')
    plot(avg_rho_RBLWs(tst_mat_i,:,:),'-+')
    plot(avg_rho_OASs(tst_mat_i,:,:),'-s')
    plot(avg_rho_Ora(tst_mat_i,:,:),'-d')
    title(['Values for \rho (Covariance ', tst_cov_mat_names{tst_mat_i}, ')' ])
    legend('LW', 'RBLW', 'OAS', 'Oracle');
    hold off
    print(['RhoValue_', tst_cov_mat_filenames{tst_mat_i}],'-dpng')
    
    figure
    hold on 
    plot(avg_MSE_LWs(tst_mat_i,:,:),'-o')
    plot(avg_MSE_RBLWs(tst_mat_i,:,:),'-+')
    plot(avg_MSE_OASs(tst_mat_i,:,:),'-s')
    plot(avg_MSE_Ora(tst_mat_i,:,:),'-d')
    title(['MSE (Covariance ', tst_cov_mat_names{tst_mat_i}, ')' ])
    legend('LW', 'RBLW', 'OAS', 'Oracle');
    hold off
    print(['MSE_', tst_cov_mat_filenames{tst_mat_i}],'-dpng')
end

% delete(findall(0,'Type','figure'))