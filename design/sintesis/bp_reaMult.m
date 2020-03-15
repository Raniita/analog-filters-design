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
## @deftypefn {} {@var{retval} =} bp_reaMult (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author:  <rani@raniita>
## Created: 2020-03-09

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
  endif
endfunction
