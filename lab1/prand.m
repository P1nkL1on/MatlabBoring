function val = prand(lambda)
% ��������� ������������� ��������� ��������
val = 0;
while val < lambda
    val = val-log(rand);
end
return;
