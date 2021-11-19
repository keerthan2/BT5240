clc
clear all

c0 = [4.49064; 0.108367; 0.00261149; 0.0112817; 0.65939; 
    0.00770135; 0.00190919; 3.62057; 0.616118; 0.299109; 
    0.0021125; 0.00422702; 0.0738334; 0.33981; 0];

tspan = [0:500];
[t,c] = ode23(@(t,x) odemodel(t,x), tspan, c0);

figure 
plot(t,c)
legend('ATP','ADP','AMP','GLC','F6P','FBP','GAP','NAD','NADH','DPG','PEP','PYR','ACA','EtOH','P','Location','Best')
xlabel('Time (seconds)')
ylabel('Concentration of Metabolites (mmol/l)')
title('Time-course simulation results of all Metabolites with all parameters set to default value')
%% Vary V4
[~,c50] = ode23(@(t,x) odemodelV450(t,x), tspan, c0);
[~,c150] = ode23(@(t,x) odemodelV4150(t,x), tspan, c0);

cnadh = c(:,9);
cnadh50 = c50(:,9);
cnadh150 = c150(:,9);

figure
plot(t,cnadh,t,cnadh50,t,cnadh150)
legend('default(100%)','50%', '150%','Location','Best')
xlabel('Time (seconds)')
ylabel('Concentration of NADH (mmol/l)')
title('Time-course simulation results of NADH concentration by changing parameter V4')
%% Varying K4GAP
[~,c50] = ode23(@(t,x) odemodelK4GAP50(t,x), tspan, c0);
[~,c150] = ode23(@(t,x) odemodelK4GAP150(t,x), tspan, c0);

cnadh = c(:,9);
cnadh50 = c50(:,9);
cnadh150 = c150(:,9);

figure
plot(t,cnadh,t,cnadh50,t,cnadh150)
legend('default(100%)','50%', '150%','Location','Best')
xlabel('Time (seconds)')
ylabel('Concentration of NADH (mmol/l)')
title('Time-course simulation results of NADH concentration by changing parameter K4GAP')

%% Varying k8f
tspan_zoom = [0:100];
tspan2 = tspan_zoom;
[~,c50] = ode23(@(t,x) odemodelk8f50(t,x), tspan2, c0);
[~,c150] = ode23(@(t,x) odemodelk8f150(t,x), tspan2, c0);

cnadh = c(:,9);
cnadh50 = c50(:,9);
cnadh150 = c150(:,9);

figure
plot(tspan2,cnadh(1+tspan2),tspan2,cnadh50,tspan2,cnadh150)
legend('default(100%)','50%', '150%','Location','Best')
xlabel('Time (seconds)')
ylabel('Concentration of NADH (mmol/l)')
title('(Zoomed) Time-course simulation results of NADH concentration by changing parameter k8f')

%% Varying k8b
[~,c50] = ode23(@(t,x) odemodelk8b50(t,x), tspan, c0);
[~,c150] = ode23(@(t,x) odemodelk8b150(t,x), tspan, c0);

cnadh = c(:,9);
cnadh50 = c50(:,9);
cnadh150 = c150(:,9);

figure
plot(t,cnadh,t,cnadh50,t,cnadh150)
legend('default(100%)','50%', '150%','Location','Best')
xlabel('Time (seconds)')
ylabel('Concentration of NADH (mmol/l)')
title('Time-course simulation results of NADH concentration by changing parameter k8b')




