% main_oncologia.m
% Pipeline principal de integracao (Bloom -> MinHash -> Naive Bayes)

clc; clear; close all;
addpath('utils');

fprintf('=== Pipeline de Triagem Genomica Computacional ===\n\n');

% 1. Carregamento e Preparacao de Dados
ficheiro_mat = 'dados_clinvar_processados.mat';
if ~isfile(ficheiro_mat)
    [dados_treino, dados_teste, resumo] = preparacao_dados_clinvar('clinvar.csv', ficheiro_mat, 5000);
else
    load(ficheiro_mat, 'dados_treino', 'dados_teste', 'resumo');
end

dados_treino(dados_treino.CLNSIG == "" | ismissing(dados_treino.CLNSIG), :) = [];
dados_teste(dados_teste.CLNSIG == "" | ismissing(dados_teste.CLNSIG), :) = [];

fprintf('Total de dados de treino: %d\n', height(dados_treino));
fprintf('Total de dados de teste: %d\n\n', height(dados_teste));

% 2. Treino dos Modelos
fprintf('--- A Treinar Modelos (Bloom, MinHash, Naive Bayes) ---\n');
t_treino_inicio = tic;

% Filtro de Bloom
n_alto_risco = sum(dados_treino.is_high_risk);
m_bloom = max(1000, round(8 * n_alto_risco));
k_bloom = 5;
filtro_bloom = FiltroBloom(m_bloom, k_bloom);
chaves_alto_risco = dados_treino.variant_key(dados_treino.is_high_risk);
for i = 1:numel(chaves_alto_risco)
    filtro_bloom.Inserir(chaves_alto_risco(i));
end

% MinHash
num_hashes = 150;
tamanho_shingle = 5;
minhash = MinHash(num_hashes, tamanho_shingle, 42);
n_treino = height(dados_treino);
assinaturas_treino = zeros(n_treino, num_hashes, 'uint64');
for i = 1:n_treino
    assinaturas_treino(i, :) = minhash.gerarAssinatura(char(dados_treino.Perfil(i)));
end

% Naive Bayes
nb = NaiveBayes(1);
nb.treinar(dados_treino, 'CLNSIG', {'Gene', 'CHROM', 'REF', 'ALT', 'Consequence'});

tempo_treino = toc(t_treino_inicio);
fprintf('Treino concluido em %.4f segundos!\n\n', tempo_treino);

% 3. Fase de Teste e Avaliacao
fprintf('--- A Iniciar Fase de Teste ---\n');

n_teste = height(dados_teste);
previsoes = strings(n_teste, 1);
classes_reais = string(dados_teste.CLNSIG);
alertas_bloom = false(n_teste, 1);

tempo_bloom = 0;
tempo_nb = 0;
tempo_minhash = 0;
contador_minhash = 0;

n_demo = min(5, n_teste);
for i = 1:n_teste
    paciente_atual = dados_teste(i, :);
    
    % Tempo Bloom
    t_bl = tic;
    alertas_bloom(i) = filtro_bloom.Verificar(paciente_atual.variant_key);
    tempo_bloom = tempo_bloom + toc(t_bl);
    
    % Tempo Naive Bayes
    t_n = tic;
    [classe_predita, confianca, ~] = nb.classificar(paciente_atual);
    tempo_nb = tempo_nb + toc(t_n);
    
    previsoes(i) = classe_predita;
    
    % MinHash (Calculado para todos, mas apenas contabilizamos tempo e mostramos nos demos)
    t_mh = tic;
    assinatura_novo = minhash.gerarAssinatura(char(paciente_atual.Perfil));
    similaridades = zeros(n_treino, 1);
    for j = 1:n_treino
        similaridades(j) = minhash.estimarJaccard(assinatura_novo, assinaturas_treino(j, :));
    end
    [valores_sim, indices_top] = maxk(similaridades, 5);
    t_mh_decorrido = toc(t_mh);
    
    if i <= n_demo
        tempo_minhash = tempo_minhash + t_mh_decorrido;
        contador_minhash = contador_minhash + 1;
        
        fprintf('--------------------------------------------------\n');
        fprintf('Paciente de Teste #%d:\n', i);
        fprintf('Gene: %s | Variacao: %s>%s\n', paciente_atual.Gene, paciente_atual.REF, paciente_atual.ALT);
        fprintf('Classe Real: %s\n\n', paciente_atual.CLNSIG);
        
        if alertas_bloom(i)
            fprintf('[Filtro Bloom] ALERTA: Variante de risco detectada!\n');
        else
            fprintf('[Filtro Bloom] Variante nao catalogada como alto risco.\n');
        end
        
        fprintf('[MinHash] Top 5 pacientes semelhantes:\n');
        for k = 1:5
            idx_historico = indices_top(k);
            paciente_historico = dados_treino(idx_historico, :);
            fprintf('  %d. Similaridade: %.1f%% | Classe: %s\n', ...
                k, valores_sim(k)*100, paciente_historico.CLNSIG);
        end
        
        fprintf('[Naive Bayes] Previsto: %s (Confianca: %.1f%%)\n', ...
            classe_predita, confianca * 100);
        fprintf('--------------------------------------------------\n\n');
    end
end

% 4. Metricas Finais
acuracia = sum(previsoes == classes_reais) / n_teste;
fprintf('=== Resultados Globais ===\n');
fprintf('Accuracy Naive Bayes: %.2f%%\n', acuracia * 100);
fprintf('Variantes marcadas no Bloom Filter: %d / %d (%.1f%%)\n\n', ...
    sum(alertas_bloom), n_teste, (sum(alertas_bloom)/n_teste)*100);

fprintf('=== Tempos Medios de Execucao (Eficiencia) ===\n');
fprintf('Tempo de Treino Total: %.4f segundos\n', tempo_treino);
fprintf('Tempo medio por consulta:\n');
fprintf('  - Filtro de Bloom (Triagem): %.6f ms\n', (tempo_bloom / n_teste) * 1000);
fprintf('  - Naive Bayes (Classificacao): %.6f ms\n', (tempo_nb / n_teste) * 1000);
if contador_minhash > 0
    fprintf('  - MinHash (Procura Top 5 em %d instâncias): %.6f ms\n', n_treino, (tempo_minhash / contador_minhash) * 1000);
end
fprintf('==============================================\n\n');

figure('Name', 'Matriz de Confusao - Naive Bayes');
confusionchart(classes_reais, previsoes, 'Title', 'Matriz de Confusao - Classificador');
