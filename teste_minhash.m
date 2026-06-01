% teste_minhash.m
clc; clear; close all;
fprintf('=== Teste Avançado do MinHash ===\n\n');

%% 1. Geração de Dados para Teste Estatístico
rng(42);
vocabulario_base = ["Gene:EGFR", "Gene:TP53", "Gene:BRCA1", "Gene:EGFR", "Chrom:1", "Chrom:17", "Chrom:13", ...
              "Ref:A", "Ref:G", "Ref:T", "Ref:C", "Alt:A", "Alt:G", "Alt:T", "Alt:C", ...
              "Consequence:missense", "Consequence:nonsense", "Consequence:synonymous", "Consequence:frameshift", ...
              "Pathogenic", "Benign", "Likely_Benign", "Uncertain_Significance"];

num_pares = 60;
jaccards_reais = zeros(num_pares, 1);
tamanho_shingle = 3;

pares_strings = cell(num_pares, 2);
for i = 1:num_pares
    tam1 = randi([4, 8]);
    tokens1 = vocabulario_base(randperm(numel(vocabulario_base), tam1));
    texto1 = char(join(tokens1, '|'));
    
    if i <= 10
        texto2 = texto1;
    elseif i <= 30
        num_partilhado = min(tam1 - 1, randi([3, tam1]));
        tokens2 = [tokens1(1:num_partilhado), vocabulario_base(randperm(numel(vocabulario_base), randi([1, 3])))];
        texto2 = char(join(tokens2, '|'));
    elseif i <= 50
        tam2 = randi([4, 8]);
        tokens2 = vocabulario_base(randperm(numel(vocabulario_base), tam2));
        texto2 = char(join(tokens2, '|'));
    else
        texto1 = 'Gene:EGFR|Chrom:1|Ref:A';
        texto2 = 'Gene:BRCA1|Chrom:17|Ref:T|Pathogenic';
    end
    
    pares_strings{i, 1} = texto1;
    pares_strings{i, 2} = texto2;
    
    minhash_temporario = MinHash(100, tamanho_shingle, 42);
    shingles1 = minhash_temporario.obterShingles(texto1);
    shingles2 = minhash_temporario.obterShingles(texto2);
    jaccards_reais(i) = length(intersect(shingles1, shingles2)) / length(union(shingles1, shingles2));
end

%% 2. Curva de Erro Médio Absoluto vs. Número de Hashes (N)
valores_N = 50:50:500;
erro_medio_absoluto = zeros(size(valores_N));

for idx = 1:numel(valores_N)
    mh = MinHash(valores_N(idx), tamanho_shingle, 42);
    erros = zeros(num_pares, 1);
    
    for i = 1:num_pares
        assinatura1 = mh.gerarAssinatura(pares_strings{i, 1});
        assinatura2 = mh.gerarAssinatura(pares_strings{i, 2});
        jaccard_estimado = mh.estimarJaccard(assinatura1, assinatura2);
        erros(i) = abs(jaccards_reais(i) - jaccard_estimado);
    end
    erro_medio_absoluto(idx) = mean(erros);
end

figure('Position',[100 100 600 400]);
plot(valores_N, erro_medio_absoluto, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('Número de hashes N');
ylabel('Erro Absoluto Médio (MAE)');
title('Erro Médio do MinHash vs. Número de Hashes N');
grid on;
saveas(gcf, 'fig_minhash_error_vs_N_new.png');

%% 3. Correlação: Real vs. Estimado (Scatter Plot para N = 150)
N_selecionado = 150;
minhash_selecionado = MinHash(N_selecionado, tamanho_shingle, 42);
jaccards_estimados = zeros(num_pares, 1);

for i = 1:num_pares
    assinatura1 = minhash_selecionado.gerarAssinatura(pares_strings{i, 1});
    assinatura2 = minhash_selecionado.gerarAssinatura(pares_strings{i, 2});
    jaccards_estimados(i) = minhash_selecionado.estimarJaccard(assinatura1, assinatura2);
end

correlacao = corrcoef(jaccards_reais, jaccards_estimados);
r_quadrado = correlacao(1,2)^2;
fprintf('--- Avaliação Estatística (N = %d, k = %d) ---\n', N_selecionado, tamanho_shingle);
fprintf('Coeficiente de Determinação (R^2): %.4f\n', r_quadrado);
fprintf('Erro Absoluto Médio Global: %.4f\n\n', mean(abs(jaccards_reais - jaccards_estimados)));

figure('Position',[100 100 600 400]);
scatter(jaccards_reais, jaccards_estimados, 40, 'd', 'filled', 'MarkerFaceColor', 'b');
hold on;
plot([0 1], [0 1], 'r--', 'LineWidth', 1.5);
hold off;
xlabel('Similaridade de Jaccard Real');
ylabel('Estimativa por MinHash');
title(sprintf('Correlação: Jaccard Real vs. Estimado (N = %d)', N_selecionado));
legend('Casos de Teste', 'Ideal (y = x)', 'Location', 'northwest');
grid on;
saveas(gcf, 'fig_minhash_k_compare_new.png');
