function [R1, C2, R3, R4, C5] = lp_reaMult (Ho, alpha, omega, ref)
  C5 = ref;
  
  k = (4*(Ho+1))/(alpha)^2;
  C2 = k*C5;
  
  raiz = sqrt(1-((4*(Ho+1))/(k*alpha^2)));
  R4 = (alpha/(2*omega*C5))*(1 + raiz); % Revisar, aqui falta +-
  R1 = R4/Ho;
  R3 = (1)/((omega^2)*(C5^2)*R4*k);
end
