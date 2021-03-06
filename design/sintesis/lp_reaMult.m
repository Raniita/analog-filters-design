function [R1, C2, R3, R4, C5] = lp_reaMult (Ho, alpha, omega, refC5, refC2)
  C5 = refC5;
  
  k_min = (4*(Ho+1))/(alpha)^2;
  %C2 = k * C5;
  C2 = refC2;
  k = C2/C5;
  
  if k < k_min
      C5 = 0;
      C2 = 0;
  end
  
  raiz = sqrt(1-(k_min*k^(-1)));
  %while(not(isreal(raiz)))
  %    k1 = k + 0.5;
  %    raiz = sqrt(1-(k*k1^(-1)));
  %end
  
  R4 = (alpha/(2*omega*C5))*(1 + raiz); 
  if(not(isreal(R4)))
    R4 = (alpha/(2*omega*C5))*(1 - raiz);
  end
  
  R1 = R4/Ho;
  R3 = ((omega^2)*(C5^2)*R4*round(k))^-1;
end