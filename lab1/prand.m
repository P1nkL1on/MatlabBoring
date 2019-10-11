function val = prand(lambda)
% генератор пуассоновской случайной величины
val = 0;
while val < lambda
    val = val-log(rand);
end
return;
