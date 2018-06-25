function [y dy] = cholera_scalar_qoi(x) 
tmesh = [0 : 0.1 : 10];
ndim = 8;
[I] = runModelCholera(x, tmesh);
%y = trapz(tmesh, I);
y = I(end);

h = 1e-4;

if nargout == 2
   E = eye(ndim);
   for i = 1 : ndim
      ei = E(:, i);

      % f(x + h * e_i)
      [I] = runModelCholera(x + h * ei, tmesh);
      %yi = trapz(tmesh, I); 
      yi = I(end);

      dy(i) = (yi - y) / h;
   end
end


   

