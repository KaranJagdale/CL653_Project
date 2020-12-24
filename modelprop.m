function propval = modelprop(x0, t1, t2)
    [tU,XpredU] = ode45(@(t,x)dyn(t,x), [t1 t2], x0);
    propval = XpredU(size(XpredU,1),:)';
end