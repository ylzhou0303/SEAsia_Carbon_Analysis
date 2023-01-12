Em = [290:2:600];
Ex = [240:5:450];

contour(Ex, Em, ones(length(Em), length(Ex)));
freq = 3382;
Ram1 = 1./(1./Ex - freq/(10^7));
Ray1 = Ex;
Ray2 = 2.*Ex;
Ram2 = 2.*Ram1;

hold on
plot(Ex, Ram1, 'r-');
hold on
plot(Ex, Ram2, 'r--');hold on
plot(Ex, Ray1, 'b-'); hold on
plot (Ex, Ray2, 'b--')

legend('','Ram1', 'Ram2', 'Ray1', 'Ray2')
xlim([240 450])
ylim([290 600])
xlabel('Ex')
ylabel('Em')