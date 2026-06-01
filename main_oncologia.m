% main_oncologia.m
% Script principal de Integracao (Bloom -> MinHash -> Naive Bayes)

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

dados_treino(dados_treino.CLNSIG == "" | ismissing(dados_treino.CLNSIG), :) = [];
dados_teste(dados_teste.CLNSIG == "" | ismissing(dados_teste.CLNSIG), :) = [];

fprintf('Total de dados de treino: %d\n', height(dados_treino));
fprintf('Total de dados de teste: %d\n\n', height(dados_teste));

%% 2. Treino dos Modelos
fprintf('--- A Treinar Modelos (Bloom, MinHash, Naive Bayes) ---\n');

% 2.1 Filtro de Bloom
n_high_risk = sum(dados_treino.is_high_risk);
m_bloom = max(1000, round(8 * n_high_risk));
k_bloom = 5;
filtro_bloom = FiltroBloom(m_bloom, k_bloom);
chaves_alto_risco = dados_treino.variant_key(dados_treino.is_high_risk);
for i = 1:numel(chaves_alto_risco)
    filtro_bloom.Inserir(chaves_alto_risco(i));
end

% 2.2 MinHash
numHashes = 150;
shingleSize = 5;
minhash = MinHash(numHashes, shingleSize, 42);
n_treino = height(dados_treino);
assinaturas_treino = zeros(n_treino, numHashes, 'uint64');
for i = 1:n_treino
    assinaturas_treino(i, :) = minhash.gerarAssinatura(char(dados_treino.Perfil(i)));
end

% 2.3 Naive Bayes
nb = NaiveBayes(1);
nb.treinar(dados_treino, 'CLNSIG', {'Gene', 'CHROM', 'REF', 'ALT', 'Consequence'});

fprintf('Treino concluido!\n\n');

%% 3. Fase de Teste e Avaliacao
fprintf('--- A Iniciar Fase de Teste ---\n');

n_teste = height(dados_teste);
previsoes = strings(n_teste, 1);
classes_reais = string(dados_teste.CLNSIG);
bloom_alertas = false(n_teste, 1);

n_demo = min(5, n_teste);
for i = 1:n_teste
    paciente_atual = dados_teste(i, :);
    
    % 1. Filtro de Bloom
    bloom_alertas(i) = filtro_bloom.Verificar(paciente_atual.variant_key);
    
    % 2. Naive Bayes
    [classe_predita, confianca, ~] = nb.classificar(paciente_atual);
    previsoes(i) = classe_predita;
    
    if i <= n_demo
        fprintf('--------------------------------------------------\n');
        fprintf('Paciente de Teste #%d:\n', i);
        fprintf('Gene: %s | Variacao: %s>%s\n', paciente_atual.Gene, paciente_atual.REF, paciente_atual.ALT);
        fprintf('Classe Real: %s\n\n', paciente_atual.CLNSIG);
        
        if bloom_alertas(i)
            fprintf('[Filtro Bloom] ALERTA: Variante de risco detectada!\n');
        else
            fprintf('[Filtro Bloom] Variante nao catalogada como alto risco.\n');
        end
        
        assinatura_novo = minhash.gerarAssinatura(char(paciente_atual.Perfil));
        similaridades = zeros(n_treino, 1);
        for j = 1:n_treino
            similaridades(j) = minhash.estimarJaccard(assinatura_novo, assinaturas_treino(j, :));
        end
        [valores_sim, indices_top] = maxk(similaridades, 5);
        
        fprintf('[MinHash] Top 5 pacientes semelhantes:\n');
        for k = 1:5
            idx_hist = indices_top(k);
            paciente_hist = dados_treino(idx_hist, :);
            fprintf('  %d. Similaridade: %.1f%% | Classe: %s\n', ...
                k, valores_sim(k)*100, paciente_hist.CLNSIG);
        end
        
        fprintf('[Naive Bayes] Previsto: %s (Confianca: %.1f%%)\n', ...
            classe_predita, confianca * 100);
        fprintf('--------------------------------------------------\n\n');
    end
end

%% 4. Metricas Finais
accuracy = sum(previsoes == classes_reais) / n_teste;
fprintf('=== Resultados Globais ===\n');
fprintf('Accuracy Naive Bayes: %.2f%%\n', accuracy * 100);
fprintf('Variantes marcadas no Bloom Filter: %d / %d (%.1f%%)\n\n', ...
    sum(bloom_alertas), n_teste, (sum(bloom_alertas)/n_teste)*100);

figure('Name', 'Matriz de Confusao - Naive Bayes');
confusionchart(classes_reais, previsoes, 'Title', 'Matriz de Confusao - Classificador');
