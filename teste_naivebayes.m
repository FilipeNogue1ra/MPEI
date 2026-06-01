% teste_naivebayes.m
clc; clear; close all;
fprintf('=== Teste Avançado do Classificador Naïve Bayes ===\n\n');

%% 1. Carregamento dos Dados Reais
ficheiro_mat = 'dados_clinvar_processados.mat';
if ~isfile(ficheiro_mat)
    error('Ficheiro de dados %s não encontrado. Por favor, execute o script principal primeiro.', ficheiro_mat);
end

load(ficheiro_mat, 'dados_treino', 'dados_teste');

dados_treino(dados_treino.CLNSIG == "" | ismissing(dados_treino.CLNSIG), :) = [];
dados_teste(dados_teste.CLNSIG == "" | ismissing(dados_teste.CLNSIG), :) = [];

dados_totais = [dados_treino; dados_teste];
n_total = height(dados_totais);
fprintf('Total de instâncias reais disponíveis para validação: %d\n\n', n_total);

%% 2. Curva de Aprendizagem (Accuracy vs. Percentagem de Treino)
proporcoes_treino = 0.5:0.1:0.9;
valores_acuracia = zeros(size(proporcoes_treino));

atributos = {'Gene', 'CHROM', 'REF', 'ALT', 'Consequence'};
alvo = 'CLNSIG';

classes_unicas = unique(string(dados_totais.(alvo)));
num_classes = numel(classes_unicas);

classes_reais_finais = [];
classes_preditas_finais = [];

for idx = 1:numel(proporcoes_treino)
    proporcao = proporcoes_treino(idx);
    
    rng(100 + idx);
    particao = cvpartition(n_total, 'HoldOut', 1 - proporcao);
    treino = dados_totais(training(particao), :);
    teste = dados_totais(test(particao), :);
    
    nb = NaiveBayes(1);
    nb.treinar(treino, alvo, atributos);
    
    n_teste = height(teste);
    predicoes = strings(n_teste, 1);
    for r = 1:n_teste
        [classe_pred, ~, ~] = nb.classificar(teste(r, :));
        predicoes(r) = classe_pred;
    end
    reais = string(teste.(alvo));
    
    acuracia = sum(predicoes == reais) / n_teste;
    valores_acuracia(idx) = acuracia;
    
    if idx == numel(proporcoes_treino)
        classes_reais_finais = reais;
        classes_preditas_finais = predicoes;
    end
end

figure('Position',[100 100 600 400]);
plot(proporcoes_treino * 100, valores_acuracia * 100, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('Percentagem de dados de treino (%)');
ylabel('Accuracy (%)');
title('Curva de Aprendizagem do Naïve Bayes (Dados Reais)');
grid on;
saveas(gcf, 'fig_nb_accuracy_vs_train_new.png');

figure('Position',[100 100 550 450]);
confusionchart(classes_reais_finais, classes_preditas_finais, ...
    'Title', sprintf('Matriz de Confusão (Partição de %.0f%% Treino)', proporcoes_treino(end)*100), ...
    'ColumnSummary', 'column-normalized', 'RowSummary', 'row-normalized');
saveas(gcf, 'fig_nb_confusion_new.png');

%% 3. Relatório de Métricas Finais
fprintf('--- Relatório de Métricas Finais (Naïve Bayes no Dataset Real) ---\n');
fprintf('Partição de Treino: %.0f%% | Partição de Teste: %.0f%%\n', proporcoes_treino(end)*100, (1-proporcoes_treino(end))*100);
fprintf('Accuracy Global: %.2f%%\n\n', valores_acuracia(end) * 100);

fprintf('Métricas Detalhadas por Classe:\n');
for c = 1:num_classes
    nome_classe = classes_unicas(c);
    
    verdadeiros_positivos = sum(classes_preditas_finais == nome_classe & classes_reais_finais == nome_classe);
    falsos_positivos = sum(classes_preditas_finais == nome_classe & classes_reais_finais ~= nome_classe);
    falsos_negativos = sum(classes_preditas_finais ~= nome_classe & classes_reais_finais == nome_classe);
    
    precisao = verdadeiros_positivos / (verdadeiros_positivos + falsos_positivos + eps);
    sensibilidade = verdadeiros_positivos / (verdadeiros_positivos + falsos_negativos + eps);
    f1 = 2 * precisao * sensibilidade / (precisao + sensibilidade + eps);
    
    fprintf(' Classe: %s\n', nome_classe);
    fprintf('   Precision: %.2f%%\n', precisao * 100);
    fprintf('   Recall (Sensibilidade): %.2f%%\n', sensibilidade * 100);
    fprintf('   F1-Score: %.4f\n', f1);
end
