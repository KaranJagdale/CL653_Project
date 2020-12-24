function derx = dyn(t,x)
    
    derx = [x(2);
            -x(4)*x(7);
            x(4);
            x(2)*x(7);
            x(6);
            0;
            0];
    
end