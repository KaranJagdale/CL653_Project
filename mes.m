function measurement = mes(x)
    measurement = [sqrt(x(1)^2 + x(3)^2 + x(5)^2);
        atan2(x(3),x(1));
        atan2(x(5),sqrt(x(1)^2 + x(3)^2))];
    if measurement(2)<0
        measurement(2) = measurement(2) + 2*pi;
    end
end