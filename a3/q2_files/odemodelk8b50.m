function xdot=odemodelk8b50(t,x)
% Define the parameter in the 'par' variable name    
    par = [0.5,0.1,0.063,1.5,0.0016,0.017,0.01,1.0,50.0,10.0,1.0,1.0,1.0,0.5,10.0,0.2,0.3,2.0,0.3,1.0,1.43E-4,10.0,10.0,0.05,0.011];
% Parameter: V1
	par_V1= par(1);
% Parameter: K1GLC
	par_K1GLC=par(2);
% Parameter: K1ATP
	par_K1ATP=par(3);
% Parameter:  V2
	par_V2=par(4);
% Parameter: K2
	par_K2=par(5);
% Parameter: k2
	par_k2=par(6);
% Parameter: K2ATP
	par_K2ATP=par(7);
% Parameter: k3f
	par_k3f=par(8);
% Parameter: k3b
	par_k3b=par(9);
% Parameter: V4
	par_V4=par(10);
% Parameter: K4GAP
	par_K4GAP=par(11);
% Parameter: K4NAD
	par_K4NAD=par(12);
% Parameter: k5f
	par_k5f=par(13);
% Parameter: k5b
	par_k5b=par(14);
% Parameter: V6
	par_V6=par(15);
% Parameter:  K6PEP
	par_K6PEP=par(16);
% Parameter: K6ADP
	par_K6ADP=par(17);
% Parameter: V7
	par_V7=par(18);
% Parameter: K7PYR
	par_K7PYR=par(19);
% Parameter: k8f
	par_k8f=par(20);
% Parameter: k8b
	par_k8b=par(21)*0.5;
% Parameter: k9f
	par_k9f=par(22);
% Parameter: k9b
	par_k9b=par(23);
% Parameter: k10
	par_k10=par(24);
% Parameter: flow
	par_flow=par(25);

% Reaction: ATPflow
	reaction_ATPflow= (3.5-x(1))*par_flow;

% Reaction: ADPflow
	reaction_ADPflow=(1.1-x(2))*par_flow;

% Reaction: NADHflow
	reaction_NADHflow=(0.24-x(9))*par_flow;

% Reaction: NADflow
	reaction_NADflow=(4-x(8))*par_flow;

% Reaction: GLCflow
	reaction_GLCflow=(50-x(4))*par_flow;

% Reaction: F6Pflow
	reaction_F6Pflow= x(5)*par_flow;

% Reaction: FBPflow
	reaction_FBPflow= x(6)*par_flow;

% Reaction: GAPflow
	reaction_GAPflow= x(7)*par_flow;

% Reaction: DPGflow
	reaction_DPGflow= x(10)*par_flow;

% Reaction: PEPflow
	reaction_PEPflow= x(11)*par_flow;

% Reaction: PYRflow
	reaction_PYRflow= x(12)*par_flow;

% Reaction: ACAflow
	reaction_ACAflow= x(13)*par_flow;

% Reaction: EtOHflow
	reaction_EtOHflow= x(14)*par_flow;

% Reaction: AMPflow
	reaction_AMPflow= x(3)*par_flow;

% Reaction: Pflow
	reaction_Pflow= x(15)*par_flow;

% Reaction: reaction_1
	reaction_reaction_1= par_V1*x(1)*x(4)/((par_K1GLC+x(4))*(par_K1ATP+x(1)));

% Reaction: reaction_2
	reaction_reaction_2= par_V2*x(1)*x(5)^2/((par_K2*(1+par_k2*(x(1)/x(3))^2)+x(5)^2)*(par_K2ATP+x(1)));

% Reaction: reaction_3
	reaction_reaction_3=(par_k3f*x(6)-par_k3b*x(7)^2);

% Reaction: reaction_4
	reaction_reaction_4=par_V4*x(8)*x(7)/((par_K4GAP+x(7))*(par_K4NAD+x(8)));

% Reaction: reaction_5
	reaction_reaction_5=(par_k5f*x(10)*x(2)-par_k5b*x(11)*x(1));

% Reaction: reaction_6
	reaction_reaction_6=par_V6*x(2)*x(11)/((par_K6PEP+x(11))*(par_K6ADP+x(2)));

% Reaction: reaction_7
	reaction_reaction_7=par_V7*x(12)/(par_K7PYR+x(12));

% Reaction: reaction_8
	reaction_reaction_8=(par_k8f*x(9)*x(13)-par_k8b*x(8)*x(14));

% Reaction: reaction_9
	reaction_reaction_9=(par_k9f*x(3)*x(1)-par_k9b*x(2)^2);

% Reaction: reaction_10
	reaction_reaction_10= par_k10*x(5);

	xdot=zeros(15,1);
	
% Species:ATP
	xdot(1) = ( 1.0 * reaction_ATPflow) + (-1.0 * reaction_reaction_1) + (-1.0 * reaction_reaction_2) + ( 1.0 * reaction_reaction_5) + ( 1.0 * reaction_reaction_6) + (-1.0 * reaction_reaction_9);
	
% Species: ADP
	xdot(2) = ( 1.0 * reaction_ADPflow) + ( 1.0 * reaction_reaction_1) + ( 1.0 * reaction_reaction_2) + (-1.0 * reaction_reaction_5) + (-1.0 * reaction_reaction_6) + ( 2.0 * reaction_reaction_9);
	
% Species: AMP
	xdot(3) = (-1.0 * reaction_AMPflow) + (-1.0 * reaction_reaction_9);
	
% Species:GLC
	xdot(4) = ( 1.0 * reaction_GLCflow) + (-1.0 * reaction_reaction_1);
	
% Species:F6P
	xdot(5) = (-1.0 * reaction_F6Pflow) + ( 1.0 * reaction_reaction_1) + (-1.0 * reaction_reaction_2) + (-1.0 * reaction_reaction_10);
	
% Species:FBP
	xdot(6) = (-1.0 * reaction_FBPflow) + ( 1.0 * reaction_reaction_2) + (-1.0 * reaction_reaction_3);
	
% Species:GAP
	xdot(7) = (-1.0 * reaction_GAPflow) + ( 2.0 * reaction_reaction_3) + (-1.0 * reaction_reaction_4);
	
% Species: NAD
	xdot(8) = ( 1.0 * reaction_NADflow) + (-1.0 * reaction_reaction_4) + ( 1.0 * reaction_reaction_8);
	
% Species: NADH
	xdot(9) = ( 1.0 * reaction_NADHflow) + ( 1.0 * reaction_reaction_4) + (-1.0 * reaction_reaction_8);
	
% Species: DPG
	xdot(10) = (-1.0 * reaction_DPGflow) + ( 1.0 * reaction_reaction_4) + (-1.0 * reaction_reaction_5);
	
% Species: PEP
	xdot(11) = (-1.0 * reaction_PEPflow) + ( 1.0 * reaction_reaction_5) + (-1.0 * reaction_reaction_6);
	
% Species: PYRw
	xdot(12) =(-1.0 * reaction_PYRflow) + ( 1.0 * reaction_reaction_6) + (-1.0 * reaction_reaction_7);
	
% Species: ACA
	xdot(13) = (-1.0 * reaction_ACAflow) + ( 1.0 * reaction_reaction_7) + (-1.0 * reaction_reaction_8);
	
% Species:EtOH
	xdot(14) = (-1.0 * reaction_EtOHflow) + ( 1.0 * reaction_reaction_8);
	
% Species:P
	xdot(15) = (-1.0 * reaction_Pflow) + ( 1.0 * reaction_reaction_10);
end