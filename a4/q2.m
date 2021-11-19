clc
clear all

load iYL1228.mat

% addpath('D:\OneDrive - IIT-Madras(IC&SR)\sem8\BT5240\cobratoolbox')

%% Retreiving IDs
model = iYL1228;

formate_bigg_id = 'EX_for_e';
glucose_bigg_id = 'EX_glc__D_e';
galactose_bigg_id = 'EX_gal_e';
maltose_bigg_id = 'EX_malt_e';

formate_id = find(ismember(model.rxns,formate_bigg_id));
glucose_id = find(ismember(model.rxns,glucose_bigg_id));
galactose_id = find(ismember(model.rxns,galactose_bigg_id));
maltose_id = find(ismember(model.rxns,maltose_bigg_id));
%% Init cobra
initCobraToolbox false
%% Wild-type biomass
sol = optimizeCbModel(model);
f_wt = sol.f
%% Glucose as source
model = iYL1228;
model.lb(glucose_id) = -30;
model.lb(galactose_id) = 0;
model.lb(maltose_id) = 0;
sol = optimizeCbModel(model);
f_glc = sol.f
[min_v max_v] = fluxVariability(model,'optPercentage',10);
v_max_formate = max_v(formate_id)
[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleGeneDeletion(model);   
lethal_genes = model.genes(find(grRatio == 0))                                               
%% Galactose as source
model = iYL1228;
model.lb(glucose_id) = 0;
model.lb(galactose_id) = -30;
model.lb(maltose_id) = 0;
sol = optimizeCbModel(model);
f_gal = sol.f
[min_v max_v] = fluxVariability(model,'optPercentage',10);
v_max_formate = max_v(formate_id)
[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleGeneDeletion(model);
lethal_genes = model.genes(find(grRatio == 0))
%% Maltose as source
model = iYL1228;
model.lb(glucose_id) = 0;
model.lb(galactose_id) = 0;
model.lb(maltose_id) = -30;
sol = optimizeCbModel(model);
f_mal = sol.f
[min_v max_v] = fluxVariability(model,'optPercentage',10);
v_max_formate = max_v(formate_id)
[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleGeneDeletion(model);
lethal_genes = model.genes(find(grRatio == 0))
%% All three sources together
model = iYL1228;
model.lb(glucose_id) = -30;
model.lb(galactose_id) = -30;
model.lb(maltose_id) = -30;
sol = optimizeCbModel(model);
f_all = sol.f
[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleGeneDeletion(model);
lethal_genes = model.genes(find(grRatio == 0))
%% Perform FVA for default source
model = iYL1228;
sol = optimizeCbModel(model);
[min_v max_v] = fluxVariability(model,'optPercentage',0);
%% FSEOF

biomassRxn = find(model.c);
targetRxn = formate_id;
f  = 0.9;
num_iter = 20;
v_max_prod = abs(max_v(targetRxn)*f);   % 90 percent of the theoretical yield
v_target_vals = [];
v_bio_vals = [];
fseof.fluxes = zeros(length(model.rxns),num_iter);
fseof.overexp = zeros(length(model.rxns),1);
fseof.ko = zeros(length(model.rxns),1);

for i=1:num_iter
    n = i*v_max_prod/num_iter;
    model.lb(targetRxn) = n;
    model.ub(targetRxn) = n;
    v_target_vals = [v_target_vals n];
    sol = optimizeCbModel(model);
    v_bio_vals = [v_bio_vals sol.f];
    fseof.fluxes(:,i) = sol.v;
    for j=1:length(fseof.fluxes)
        if fseof.fluxes(j,1) > 0
            if i == 1   
                fseof.overexp(j,1)=1;
                fseof.ko(j,1)=1;
            else
                if (fseof.fluxes(j,i) > fseof.fluxes(j,i-1)) & fseof.overexp(j,1)
                    fseof.overexp(j,1)=1;
                else
                    fseof.overexp(j,1)=0;
                end
                if (fseof.fluxes(j,i) < fseof.fluxes(j,i-1)) & fseof.ko(j,1)
                    fseof.ko(j,1)=1;
                else
                    fseof.ko(j,1)=0;
                end
            end
            
        elseif fseof.fluxes(j,1) < 0
            if i == 1
                fseof.overexp(j,1)=1;
                fseof.ko(j,1)=1;
            else
                if (fseof.fluxes(j,i) < fseof.fluxes(j,i-1)) & fseof.overexp(j,1)
                    fseof.overexp(j,1)=1;
                else
                    fseof.overexp(j,1)=0;
                end
                if (fseof.fluxes(j,i) > fseof.fluxes(j,i-1)) & fseof.ko(j,1)
                    fseof.ko(j,1)=1;
                else
                    fseof.ko(j,1)=0;
                end
            end
            
        end
        
    end
end
%% Printing

print_format = '%s\t%s\t%s\t%s\n';
fprintf(print_format,'rxnNumber','rxnID','rxnName','Type of target');
for num=1:length(fseof.ko)
    if fseof.ko(num,1) == 1
        c1 = num2str(num);                                  
        c2 = char(model.rxns(num));                         
        c3 = char(model.rxnNames(num));
        c4 = 'Potential Knockout';
        fprintf(print_format,c1,c2,c3,c4);
    end
end
disp(' ')
disp('################################################################')
disp(' ')
fprintf(print_format,'rxnNumber','rxnID','rxnName','Type of target');
for num=1:length(fseof.overexp)
    if fseof.overexp(num,1) == 1
        c1 = num2str(num);                                  
        c2 = char(model.rxns(num));                         
        c3 = char(model.rxnNames(num));
        c4 = 'Over-Expression';
        fprintf(print_format,c1,c2,c3,c4);
    end
end
%% Some plots
reac_num = 1830; % knockouts: 130,140,781,823,1023 | over-expression: 665,797,1538,1544,1830
reac_flux = abs(fseof.fluxes(reac_num,:));
figure
% scatter(v_target_vals,reac_flux)
% xlabel('Formate flux')
% ylabel('abs('+string(model.rxnNames(reac_num))+' flux)')
scatter(v_bio_vals,v_target_vals)
xlabel('Biomass flux')
ylabel('Formate flux')
