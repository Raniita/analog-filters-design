function [C1, R2, C3, C4, R5] = hp_reaMult (Ho, alpha, omega, ref)
  C1 = ref;
  C3 = ref;
  
  R5 = ((1)/(alpha*omega*C1))*(2*Ho+1);
  R2 = (alpha*Ho)/(omega*C1*(2*Ho+1));
  C4 = C1/Ho;
end
