function bin = digito(l1,l2,l3,l4,c1,c2,c3)

bin = zeros(size(l1,1),4);
for i = 1:size(l1,1)
decimal = bi2de([l1(i) l2(i) l3(i) l4(i) c1(i) c2(i) c3(i)]);
verifica = sum([l1(i) l2(i) l3(i) l4(i) c1(i) c2(i) c3(i)]);

switch (decimal)
    case 0 
        bin(1,:) = [0 0 0 0]; % 0   % nada foi selecionado
    case 17
        bin(i,:) = [1 0 0 0]; % 1   % linha 1 e coluna 1
    case 33
        bin(i,:) = [0 1 0 0]; % 2   % linha 1 e coluna 2
    case 65
        bin(i,:) = [1 1 0 0]; % 3   % linha 1 e coluna 3
    case 18
        bin(i,:) = [0 0 1 0]; % 4   % linha 2 e coluna 1
    case 34
        bin(i,:) = [1 0 1 0]; % 5   % linha 2 e coluna 2
    case 66
        bin(i,:) = [0 1 1 0]; % 6   % linha 2 e coluna 3
    case 20
        bin(i,:) = [1 1 1 0]; % 7   % linha 3 e coluna 1
    case 36
        bin(i,:) = [0 0 0 1]; % 8   % linha 3 e coluna 2
    case 68
        bin(i,:) = [1 0 0 1]; % 9   % linha 3 e coluna 3
    case 24 
        bin(i,:) = [0 1 0 1]; % 10  % linha 4 e coluna 2
    case 40
        bin(i,:) = [0 0 1 1]; %*    % linha 4 e coluna 1
    case 72
        bin(i,:) = [0 0 1 1]; % #   % linha 4 e coluna 3
   otherwise
        switch (verifica) % Erro
            case 3
                bin(i,:) = [1 0 1 1];
            case 2
                bin(i,:) = [0 1 1 1];
            otherwise
                bin(i,:) = [1 1 1 1];
        end
    end
end
