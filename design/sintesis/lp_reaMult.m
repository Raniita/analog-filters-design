## Copyright (C) 2020 
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} lp_reaMult (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author:  <rani@raniita>
## Created: 2020-03-09

function [R1, C2, R3, R4, C5] = lp_reaMult (Ho, alpha, omega, ref)
  C5 = ref;
  
  k = (4*(Ho+1))/(alpha)^2;
  C2 = k*C5;
  
  raiz = sqrt(1-((4*(Ho+1))/(k*alpha^2)));
  R4 = (alpha/(2*omega*C5))*(1 + raiz); % Revisar, aqui falta +-
  R1 = R4/Ho;
  R3 = (1)/((omega^2)*(C5^2)*R4*k);
endfunction
