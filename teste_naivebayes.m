% teste_naivebayes.m
% Teste de validacao do Classificador Naive Bayes

clc; clear; close all;

fprintf('=== Teste Isolado do Naive Bayes ===\n\n');

Gene = ["SAMD11"; "SAMD11"; "EGFR"; "EGFR"; "TP53"; "TP53"];
Consequence = ["missense"; "synonymous"; "missense"; "missense"; "nonsense"; "synonymous"];
CLNSIG = ["Pathogenic"; "Benign"; "Pathogenic"; "Pathogenic"; "Pathogenic"; "Benign"];

treino = table(Gene, Consequence, CLNSIG);
disp('Dados de Treino Mock:');
disp(treino);

nb = NaiveBayes(1);
nb.treinar(treino, 'CLNSIG', {'Gene', 'Consequence'});

amostra1 = table("SAMD11", "missense", 'VariableNames', {'Gene', 'Consequence'});
[classe1, conf1, ~] = nb.classificar(amostra1);

fprintf('Previsao para (Gene=SAMD11, Consequence=missense):\n');
fprintf('  Classe Prevista: %s\n', classe1);
fprintf('  Confianca: %.2f%%\n\n', conf1 * 100);

amostra2 = table("BRCA1", "missense", 'VariableNames', {'Gene', 'Consequence'});
[classe2, conf2, ~] = nb.classificar(amostra2);

fprintf('Previsao para Novo Gene (Gene=BRCA1, Consequence=missense):\n');
fprintf('  Classe Prevista: %s\n', classe2);
fprintf('  Confianca: %.2f%%\n', conf2 * 100);
