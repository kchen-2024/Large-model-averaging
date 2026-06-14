clear all
%% U.S. Crime dataset
data = readtable('UScrime.csv');
X = table2array(data(:, 2:end-1));  % Predictive variables (15)
y = data.y;                         % Response variable
n = size(X, 1);                     % Sample size (47)
var_names = data.Properties.VariableNames(2:end-1); 
%% Motor Trend Car dataset
% data = readtable('mtcars_dataset.csv');
% X = table2array(data(:, 3:end));    % Predictive variables (10)
% y = data.mpg;                       % Response variable 
% n = size(X, 1);                     % Sample size (32)
% var_names = data.Properties.VariableNames(3:end); 

%% Forward selection method (based on Mallows Cp)
full_model = fitlm(X, y);
MSE_full = full_model.MSE;  % MSE for the full model
selected = [];        
remaining = 1:size(X, 2);      
priority_order = {};  
cp_values = [];        
fprintf('Mallows Cp variable priority ranking:\n');
fprintf('%-5s %-10s %-8s\n', 'Steps', 'Add variables', 'Cp values');
fprintf('----------------------------\n');

for step = 1:size(X, 2)
    best_cp = inf;
    best_idx = 0;
    %% Add each remaining variable
    for j = 1:length(remaining)
        idx = [selected, remaining(j)];
        X_sub = X(:, idx);
        %% Fit the model and calculate Cp
        model = fitlm(X_sub, y);
        SSE_k = sum(model.Residuals.Raw.^2);
        p_k = size(X_sub, 2)+1;
        cp_val = (SSE_k / MSE_full) - (n - 2*p_k);      
        %% Update the best Cp
        if cp_val < best_cp
            best_cp = cp_val;
            best_idx = remaining(j);
        end
    end    
    %% Update the variable list
    selected = [selected, best_idx];
    remaining(remaining == best_idx) = [];
    priority_order{end+1} = var_names{best_idx};
    cp_values(end+1) = best_cp;
    fprintf('%-5d %-10s %-8.2f\n', step, var_names{best_idx}, best_cp);
end

[~, order] = ismember(priority_order, var_names);
X_ordered = X(:, order);
fprintf('Original variable order:\n');
disp(var_names);
fprintf('\nFinal variable priority ranking (from high to low):\n');
for i = 1:length(priority_order)
    fprintf('%d. %s\n', i, priority_order{i});
end

%% Create a sorted table
var_names_ordered=priority_order;
var_names_ordered{end+1}='y';
ordered_table = array2table([X_ordered,y], 'VariableNames', var_names_ordered);
% writetable(ordered_table, 'ordered_UScrime.csv')
% writetable(ordered_table, 'ordered_mtcars.csv')