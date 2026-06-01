% main_oncologia.m
% Script principal de Integracao (MinHash)

clc; clear; close all;
addpath('utils');

fprintf('=== Pipeline de Triagem Genomica Computacional ===\n\n');

%% 1. Carregamento e Preparacao de Dados
dados_mat = 'dados_clinvar_processados.mat';
if ~isfile(dados_mat)
    [dados_treino, dados_teste, resumo] = preparacao_dados_clinvar('clinvar.csv', dados_mat, 5000);
else
    load(dados_mat, 'dados_treino', 'dados_teste', 'resumo');
end

fprintf('Total de dados de treino: %d\n', height(dados_treino));
fprintf('Total de dados de teste: %d\n\n', height(dados_teste));

%% 2. Indexacao MinHash (Treino)
numHashes = 150;
shingleSize = 5;
minhash = MinHash(numHashes, shingleSize, 42);

n_treino = height(dados_treino);
assinaturas_treino = zeros(n_treino, numHashes, 'uint64');

for i = 1:n_treino
    perfilStr = char(dados_treino.Perfil(i));
    assinaturas_treino(i, :) = minhash.gerarAssinatura(perfilStr);
end

%% 3. Consulta MinHash (Teste)
n_teste_demo = min(5, height(dados_teste));
fprintf('Demonstracao com os primeiros %d pacientes de teste:\n\n', n_teste_demo);

for i = 1:n_teste_demo
    paciente_atual = dados_teste(i, :);
    perfil_novo = char(paciente_atual.Perfil);
    
    fprintf('--------------------------------------------------\n');
    fprintf('Paciente de Teste #%d:\n', i);
    fprintf('Gene: %s | Variacao: %s> %s\n', paciente_atual.Gene, paciente_atual.REF, paciente_atual.ALT);
    fprintf('Classe Real: %s\n', paciente_atual.CLNSIG);
    
    assinatura_novo = minhash.gerarAssinatura(perfil_novo);
    
    similaridades = zeros(n_treino, 1);
    for j = 1:n_treino
        similaridades(j) = minhash.estimarJaccard(assinatura_novo, assinaturas_treino(j, :));
    end
    
    [valores_sim, indices_top] = maxk(similaridades, 3);
    
    fprintf('\nResultados MinHash (Top 3):\n');
    for k = 1:3
        idx_hist = indices_top(k);
        paciente_hist = dados_treino(idx_hist, :);
        fprintf('  %d. Similaridade: %.1f%% | Classe: %s | Gene: %s\n', ...
            k, valores_sim(k)*100, paciente_hist.CLNSIG, paciente_hist.Gene);
    end
    fprintf('--------------------------------------------------\n\n');
end
