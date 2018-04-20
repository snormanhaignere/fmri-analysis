function tcout = hrf_double_gamma(tcin,a1,b1,d1,a2,b2,d2,s2)

t1 = max(tcin-d1,0);t2 = max(tcin-d2,0);
gam1 = gampdf(t1,a1,b1);
gam2 = gampdf(t2,a2,b2);
tcout = gam1-gam2*s2;



% tcout = (((tt/tau).^(n-1)).*exp(-tt/tau))/(tau*factorial(n-1));

