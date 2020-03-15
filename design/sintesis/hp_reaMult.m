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
## @deftypefn {} {@var{retval} =} hp_reaMult (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author:  <rani@raniita>
## Created: 2020-03-09

function [C1, R2, C3, C4, R5] = hp_reaMult (Ho, alpha, omega, ref)
  C1 = ref;
  C3 = ref;
  
  R5 = ((1)/(alpha*omega*C1))*(2*Ho+1);
  R2 = (alpha*Ho)/(omega*C1*(2*Ho+1));
  C4 = C1/Ho;
endfunction
