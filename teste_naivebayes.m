% teste_naivebayes.m
% Teste avançado e rigoroso do classificador Naïve Bayes com dados reais (MPEI)

clc; clear; close all;
fprintf('=== Teste Avançado do Classificador Naïve Bayes ===\n\n');

%% 1. Carregamento dos Dados Reais
dados_mat = 'dados_clinvar_processados.mat';
if ~isfile(dados_mat)
    error('Ficheiro de dados %s não encontrado. Por favor, execute o script principal primeiro.', dados_mat);
end

load(dados_mat, 'dados_treino', 'dados_teste');

% Limpeza de dados com classes em falta
dados_treino(dados_treino.CLNSIG == "" | ismissing(dados_treino.CLNSIG), :) = [];
dados_teste(dados_teste.CLNSIG == "" | ismissing(dados_teste.CLNSIG), :) = [];

% Juntar os dados para podermos fazer partições controladas
dados_totais = [dados_treino; dados_teste];
n_total = height(dados_totais);
fprintf('Total de instâncias reais disponíveis para validação: %d\n\n', n_total);

%% 2. Curva de Aprendizagem (Accuracy vs. Percentagem de Treino)
propTrain = 0.5:0.1:0.9; % Variamos de 50% a 90% de treino
accVals = zeros(size(propTrain));
precisionVals = zeros(size(propTrain));
recallVals = zeros(size(propTrain));
F1Vals = zeros(size(propTrain));

features = {'Gene', 'CHROM', 'REF', 'ALT', 'Consequence'};
target = 'CLNSIG';

% Classes únicas
classes_unicas = unique(string(dados_totais.(target)));
num_classes = numel(classes_unicas);

% Para guardar a matriz de confusão da última iteração (90% treino)
confMatFinal = [];
classesReaisFinal = [];
classesPreditasFinal = [];

for idx = 1:numel(propTrain)
    p = propTrain(idx);
    
    % Divisão Treino/Teste
    rng(100 + idx); % Semente diferente para cada divisão para diversidade
    cv = cvpartition(n_total, 'HoldOut', 1 - p);
    tr = dados_totais(training(cv), :);
    te = dados_totais(test(cv), :);
    
    % Treinar Classificador
    nb = NaiveBayes(1); % Suavização de Laplace alpha = 1
    nb.treinar(tr, target, features);
    
    % Classificar conjunto de teste
    n_test = height(te);
    predicoes = strings(n_test, 1);
    for r = 1:n_test
        [classe_pred, ~, ~] = nb.classificar(te(r, :));
        predicoes(r) = classe_pred;
    end
    reais = string(te.(target));
    
    % Métricas Globais
    accuracy = sum(predicoes == reais) / n_test;
    accVals(idx) = accuracy;
    
    % Guardar dados da última partição para o gráfico de matriz de confusão
    if idx == numel(propTrain)
        classesReaisFinal = reais;
        classesPreditasFinal = predicoes;
    end
end

% ---------- Gráfico 5 – Curva de Aprendizagem (Accuracy vs. Percentagem de Treino) ----------
figure('Position',[100 100 600 400]);
plot(propTrain * 100, accVals * 100, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('Percentagem de dados de treino (%)');
ylabel('Accuracy (%)');
title('Curva de Aprendizagem do Naïve Bayes (Dados Reais)');
grid on;
saveas(gcf, 'fig_nb_accuracy_vs_train_new.png');

% ---------- Gráfico 6 – Matriz de Confusão (heatmap real) ----------
figure('Position',[100 100 550 450]);
confusionchart(classesReaisFinal, classesPreditasFinal, ...
    'Title', sprintf('Matriz de Confusão (Partição de %.0f%% Treino)', propTrain(end)*100), ...
    'ColumnSummary', 'column-normalized', 'RowSummary', 'row-normalized');
saveas(gcf, 'fig_nb_confusion_new.png');

%% 3. Relatório de Métricas Finais
fprintf('--- Relatório de Métricas Finais (Naïve Bayes no Dataset Real) ---\n');
fprintf('Partição de Treino: %.0f%% | Partição de Teste: %.0f%%\n', propTrain(end)*100, (1-propTrain(end))*100);
fprintf('Accuracy Global: %.2f%%\n\n', accVals(end) * 100);

% Métricas por classe
fprintf('Métricas Detalhadas por Classe:\n');
for c = 1:num_classes
    class_name = classes_unicas(c);
    
    % Binário para esta classe
    tp = sum(classesPreditasFinal == class_name & classesReaisFinal == class_name);
    fp = sum(classesPreditasFinal == class_name & classesReaisFinal ~= class_name);
    fn = sum(classesPreditasFinal ~= class_name & classesReaisFinal == class_name);
    
    precision = tp / (tp + fp + eps);
    recall = tp / (tp + fn + eps);
    f1 = 2 * precision * recall / (precision + recall + eps);
    
    fprintf(' Classe: %s\n', class_name);
    fprintf('   Precision: %.2f%%\n', precision * 100);
    fprintf('   Recall (Sensibilidade): %.2f%%\n', recall * 100);
    fprintf('   F1-Score: %.4f\n', f1);
end
