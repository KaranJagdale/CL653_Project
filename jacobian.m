function jac = jacobian(x)
    jac = [0 1 0 0 0 0 0;
           0 0 0 -x(7) 0 0 -x(4);
           0 0 0 1 0 0 0;
           0 x(7) 0 0 0 0 x(2);
           0 0 0 0 0 1 0;
           0 0 0 0 0 0 0;
           0 0 0 0 0 0 0];
end