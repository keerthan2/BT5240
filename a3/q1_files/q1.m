clc
clear all

v0 = [2 100 2 2];
data = table2array(readtable('data.csv'));

v = fminsearch(@(v) loss_func(v,data),v0);

disp('Estimated Vm after minimizing the objective function: ')
fprintf('--> Vm1 = %f \n', v(1))
fprintf('--> Vm2 = %f \n', v(2))
fprintf('--> Vm3 = %f \n', v(3))
fprintf('--> Vm4 = %f \n', v(4))

[t, y] = pred(v); 

function loss = loss_func(v,data)
    [t,y] = pred(v);
    loss = norm((data-y),'fro');
end

function [t,y] = pred(v)
    cap0=10; cnap0=0.75; c20=300; ca0=0; cs0=0; cna0=0;
    k = [0.034;0.112;1.091;0.408;0.004;
        0.007;0.07;0.899;10.973;0.017;
        6.077;0.01;7.576;0.079;0.006;
        20.303;0.535;1.886;0.005;2.242;
        0.348;29.551;0.007;288.016;
        ];
    tspan = [0:20];
    y0 = [cap0 cs0 ca0 c20 cna0 cnap0];
    [t,y] = ode23(@(t,y) solver(t,y,v,k), tspan, y0);
end

function dydt = solver(t,y,v,k)
    dydt = zeros(6,1);
    
    cap=y(1); cs=y(2); ca=y(3); c2=y(4); cna=y(5); cnap=y(6);
    
    vm1 = v(1); vm2 = v(2); vm3 = v(3); vm4 = v(4);
    
    kmap=k(1); kms=k(2); kma=k(3); km2=k(4); kmna=k(5);
    kmnap=k(6); ki1s=k(7); ki1a=k(8); ki12=k(9); ki1nap=k(10);
    ki1ap=k(11); ki2ap=k(12); ki2a=k(13); ki22=k(14); ki2na=k(15);
    ki32=k(16); ki3s=k(17); ki3ap=k(18); ki3nap=k(19); ki4a=k(20);
    ki4ap=k(21); ki4s=k(22); ki4na=k(23); ki42=k(24);
    
    r1_den = (cap+kmap*(1+cs/ki1s+ca/ki1a+c2/ki12)+(cap.^2/ki1ap)).*(cna+kmna*(1+cnap/ki1nap));
    r1 = vm1*cap.*cna./r1_den;
    
    r2_den = (cs+kms*(1+cap/ki2ap+ca/ki2a+c2/ki22)).*(cnap+kmnap*(1+cna/ki2na));
    r2 = vm2*cs.*cnap./r2_den;
    
    r3_den = (ca+kma*(1+cs/ki3s+cap/ki3ap+c2/ki32)).*(cna+kmna*(1+cnap/ki3nap));
    r3 = vm3*ca.*cna./r3_den;
    
    r4_den = (c2+km2*(1+cs/ki4s+ca/ki4a+cap/ki4ap)+(c2.^2/ki42)).*(cnap+kmnap*(1+cna/ki4na));
    r4 = vm4*c2.*cnap./r4_den;
    
    dydt(1,1) = -r1 + r2;
    dydt(2,1) = r1 - r2;
    dydt(3,1) = -r3 + r4;
    dydt(4,1) = r3 - r4;
    dydt(5,1) = -r1 + r2 - r3 + r4;
    dydt(6,1) = r1 - r2 + r3 - r4;

end