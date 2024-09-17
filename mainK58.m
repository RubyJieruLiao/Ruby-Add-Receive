
tau0 = 6667; 
patient = 'P1';
connectivity_path = strcat('../data/connectivity_', patient, '/');
K = load(strcat(connectivity_path, 'weights.txt')); 
K = normal(K); 
[N, ~] = size(K); 


results = cell(1, 85);
node_indices = cell(1, 85);
results{1} = K;  
node_indices{1} = 1:N; 

for i = 1:84
    K_modified = K;  
    K_modified(i, :) = []; 
    K_modified(:, i) = [];  
    results{i+1} = K_modified;  
    node_indices{i+1} = setdiff(1:N, i); 
end

format long
x0_range1 = -2.3;
x0_EZ_range1 = -1.5:0.01:-1
disp(x0_EZ_range1);

StabilityMatrix = NaN(1, N); 


EZ_values = [];
for orig_EZ = [23, 27, 72, 76, 78]
    fprintf('Processing original EZ = %d\n', orig_EZ);

    for k_idx = 59  
        K_current = results{k_idx};
        [N_current, ~] = size(K_current); 

        current_indices = node_indices{k_idx}; 
        if ismember(orig_EZ, current_indices)
            EZ = find(current_indices == orig_EZ); 
            EZ_values = [EZ_values, orig_EZ]; 
            failed_x0_EZ = NaN; 
            last_real_x0_EZ = NaN; 
            found_transition = false;

            for ix0_EZ = 1:length(x0_EZ_range1)
                x0 = x0_range1 + zeros(N_current, 1);
                x0(EZ) = x0_EZ_range1(ix0_EZ);
                Z0 = 3 + zeros(N_current, 1); 
                one_dim_epileptor_fun = @(z) oneDepileptor(z, x0, K_current, tau0); 
                opt = optimset('TolFun', 1e-14, 'TolX', 1e-14);
                [z_fixed_numerical, fval, exitflag, output] = fsolve(one_dim_epileptor_fun, Z0, opt);

                fprintf('Exit flag: %d\n', exitflag);


                if exitflag < 0
                    warning('Failed to find fixed point, skipping current x0(EZ) value.');
                    continue; 
                end


                C1 = CouplingMatrix(z_fixed_numerical, K_current, tau0);


                try
                    [covariance1, err] = con2cov(C1, false, 10000000, 10, 1);
                    stability1 = Stability(covariance1); 
                catch ME
    
                    warning('Error computing stability1: %s', ME.message);
                    stability1 = NaN;
                end


                if isreal(stability1) && ~isnan(stability1)
                

                    fprintf('stability1 = %f at x0(EZ) = %f\n', stability1, x0_EZ_range1(ix0_EZ));
                    last_real_x0_EZ = x0_EZ_range1(ix0_EZ);
                    continue;
                else
               
                    failed_x0_EZ = x0_EZ_range1(ix0_EZ);
                    found_transition = true;
                    warning('stability1 is not a real number or is NaN at x0(EZ) = %f', failed_x0_EZ);
     
                    fine_x0_EZ_range = linspace(last_real_x0_EZ, failed_x0_EZ, 11); 
                    for fine_x0_EZ = fine_x0_EZ_range(2:end-1) 
                        x0(EZ) = fine_x0_EZ; 

                  
                        Z0 = 3 + zeros(N_current, 1);
                        one_dim_epileptor_fun = @(z) oneDepileptor(z, x0, K_current, tau0);
                        [z_fixed_numerical, fval, exitflag, output] = fsolve(one_dim_epileptor_fun, Z0, opt);

                        if exitflag < 0
                            continue;
                        end

                        C1 = CouplingMatrix(z_fixed_numerical, K_current, tau0);

                        try
                            [covariance1, err] = con2cov(C1, false, 10000000, 10, 1);
                            stability1 = Stability(covariance1); 
                        catch
                            stability1 = NaN;
                        end

                        %if isreal(stability1) && ~isnan(stability1)
                        if ~isnan(stability1)
                            fprintf('Fine stability1 = %f at x0(EZ) = %f\n', stability1, fine_x0_EZ);
                            last_real_x0_EZ = fine_x0_EZ;
                        else
                            fprintf('Fine stability1 is not real or is NaN at x0(EZ) = %f\n', fine_x0_EZ);
                            failed_x0_EZ = fine_x0_EZ;
                            break;
                        end
                    end
                    break;
                end
            end

     
            StabilityMatrix(EZ) = failed_x0_EZ;
        end
    end
end


EZ_range_min = min(EZ_values);
EZ_range_max = max(EZ_values); 

filename = sprintf('StabilityMatrix_K%d_EZ%d_to_%d.mat', k_idx-1, EZ_range_min, EZ_range_max);


save(filename, 'StabilityMatrix');
disp(['Saved ', filename, ' successfully.']);
