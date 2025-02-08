% Gallagher
function result = c_GGG(T,B)


result = (GGG_ENTROPY2(T+0.0001,B)-GGG_ENTROPY2(T,B))/0.0001*T;
end
