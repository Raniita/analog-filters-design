function [R1, R2, C3, C4, R5] = bp_reaMult (Ho, alpha, omega, ref)
  C3 = ref;
  C4 = ref;
  
  Q = 1/alpha;
  if Ho > 2*Q^2
    R1 = (Q)/(Ho*omega*C3);
    R2 = (Q)/((2*Q-Ho)*omega*C3);
    R5 = (2*Q)/(omega*C3);
  else
    R1 = 0;
    R2 = 0;
    R5 = 0;
  end
  end
