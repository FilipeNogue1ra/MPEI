% teste_minhash.m
% Teste avançado e rigoroso do algoritmo MinHash (MPEI)

clc; clear; close all;
fprintf('=== Teste Avançado do MinHash ===\n\n');

%% 1. Geração de Dados para Teste Estatístico
% Vamos gerar pares de perfis genéticos sintéticos com diferentes níveis de sobreposição
rng(42); % Para reprodutibilidade
base_vocab = ["Gene:EGFR", "Gene:TP53", "Gene:BRCA1", "Gene:EGFR", "Chrom:1", "Chrom:17", "Chrom:13", ...
              "Ref:A", "Ref:G", "Ref:T", "Ref:C", "Alt:A", "Alt:G", "Alt:T", "Alt:C", ...
              "Consequence:missense", "Consequence:nonsense", "Consequence:synonymous", "Consequence:frameshift", ...
              "Pathogenic", "Benign", "Likely_Benign", "Uncertain_Significance"];

num_pairs = 60;
realJaccards = zeros(num_pairs, 1);
shingle_size = 3;

% Armazenar os pares de strings
strPairs = cell(num_pairs, 2);
for i = 1:num_pairs
    % Par 1: selecionar k elementos aleatórios
    len1 = randi([4, 8]);
    tokens1 = base_vocab(randperm(numel(base_vocab), len1));
    s1 = char(join(tokens1, '|'));
    
    % Par 2: partilhar alguns elementos para simular similaridade variada
    if i <= 10
        % Identicos
        s2 = s1;
    elseif i <= 30
        % Alta similaridade (partilha a maioria dos tokens)
        num_share = min(len1 - 1, randi([3, len1]));
        tokens2 = [tokens1(1:num_share), base_vocab(randperm(numel(base_vocab), randi([1, 3])))];
        s2 = char(join(tokens2, '|'));
    elseif i <= 50
        % Baixa similaridade
        len2 = randi([4, 8]);
        tokens2 = base_vocab(randperm(numel(base_vocab), len2));
        s2 = char(join(tokens2, '|'));
    else
        % Completamente disjuntos (vocabulários diferentes)
        s1 = 'Gene:EGFR|Chrom:1|Ref:A';
        s2 = 'Gene:BRCA1|Chrom:17|Ref:T|Pathogenic';
    end
    
    strPairs{i, 1} = s1;
    strPairs{i, 2} = s2;
    
    % Calcular Jaccard real usando shingles
    mh_temp = MinHash(100, shingle_size, 42);
    sh1 = mh_temp.obterShingles(s1);
    sh2 = mh_temp.obterShingles(s2);
    realJaccards(i) = length(intersect(sh1, sh2)) / length(union(sh1, sh2));
end

%% 2. Curva de Erro Médio Absoluto vs. Número de Hashes (N)
Ns = 50:50:500;
meanAbsError = zeros(size(Ns));

for idx = 1:numel(Ns)
    mh = MinHash(Ns(idx), shingle_size, 42);
    errors = zeros(num_pairs, 1);
    
    for i = 1:num_pairs
        sig1 = mh.gerarAssinatura(strPairs{i, 1});
        sig2 = mh.gerarAssinatura(strPairs{i, 2});
        estJ = mh.estimarJaccard(sig1, sig2);
        errors(i) = abs(realJaccards(i) - estJ);
    end
    meanAbsError(idx) = mean(errors);
end

% Gráfico 3 – Erro Absoluto Médio vs N
figure('Position',[100 100 600 400]);
plot(Ns, meanAbsError, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('Número de hashes N');
ylabel('Erro Absoluto Médio (MAE)');
title('Erro Médio do MinHash vs. Número de Hashes N');
grid on;
saveas(gcf, 'fig_minhash_error_vs_N_new.png');


%% 3. Correlação: Real vs. Estimado (Scatter Plot para N = 150)
N_selected = 150; % Valor do pipeline principal
mh_select = MinHash(N_selected, shingle_size, 42);
estimatedJaccards = zeros(num_pairs, 1);

for i = 1:num_pairs
    sig1 = mh_select.gerarAssinatura(strPairs{i, 1});
    sig2 = mh_select.gerarAssinatura(strPairs{i, 2});
    estimatedJaccards(i) = mh_select.estimarJaccard(sig1, sig2);
end

% Calcular coeficiente de determinação R² para mostrar no console
correlation = corrcoef(realJaccards, estimatedJaccards);
r_squared = correlation(1,2)^2;
fprintf('--- Avaliação Estatística (N = %d, k = %d) ---\n', N_selected, shingle_size);
fprintf('Coeficiente de Determinação (R^2): %.4f\n', r_squared);
fprintf('Erro Absoluto Médio Global: %.4f\n\n', mean(abs(realJaccards - estimatedJaccards)));

% Gráfico 4 – Comparação de desempenho do MinHash (Scatter Plot)
figure('Position',[100 100 600 400]);
scatter(realJaccards, estimatedJaccards, 40, 'd', 'filled', 'MarkerFaceColor', 'b');
hold on;
% Linha de identidade y = x
plot([0 1], [0 1], 'r--', 'LineWidth', 1.5);
hold off;
xlabel('Similaridade de Jaccard Real');
ylabel('Estimativa por MinHash');
title(sprintf('Correlação: Jaccard Real vs. Estimado (N = %d)', N_selected));
legend('Casos de Teste', 'Ideal (y = x)', 'Location', 'northwest');
grid on;
saveas(gcf, 'fig_minhash_k_compare_new.png'); % Substitui fig_minhash_k_compare por este scatter plot de validação
