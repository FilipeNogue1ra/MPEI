% teste_bloom.m
% Teste de validacao do Filtro de Bloom

clc; clear; close all;

fprintf('=== Teste Isolado do Filtro de Bloom ===\n\n');

m = 8000;
k = 3;
filtro = FiltroBloom(m, k);

n_treino = 1000;
dados_treino = cell(n_treino, 1);
for i = 1:n_treino
    dados_treino{i} = sprintf('VAR_%d', randi(1000000));
    filtro.Inserir(dados_treino{i});
end

erros_negativos = 0;
for i = 1:n_treino
    if ~filtro.Verificar(dados_treino{i})
        erros_negativos = erros_negativos + 1;
    end
end
fprintf('Falsos Negativos Encontrados (deve ser 0): %d\n', erros_negativos);

n_teste = 10000;
falsos_positivos = 0;
for i = 1:n_teste
    elemento_teste = sprintf('TEST_%d', randi(1000000));
    while ismember(elemento_teste, dados_treino)
        elemento_teste = sprintf('TEST_%d', randi(1000000));
    end
    
    if filtro.Verificar(elemento_teste)
        falsos_positivos = falsos_positivos + 1;
    end
end

taxa_empirica = falsos_positivos / n_teste;
taxa_teorica = (1 - exp(-k * n_treino / m))^k;

fprintf('Falsos Positivos Detectados: %d em %d testes\n', falsos_positivos, n_teste);
fprintf('Taxa Empirica: %.4f (%.2f%%)\n', taxa_empirica, taxa_empirica * 100);
fprintf('Taxa Teorica:  %.4f (%.2f%%)\n', taxa_teorica, taxa_teorica * 100);
fprintf('Diferenca:     %.4f\n', abs(taxa_empirica - taxa_teorica));
