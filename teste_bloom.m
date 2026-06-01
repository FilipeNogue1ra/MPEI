% teste_bloom.m
% Teste avançado e rigoroso do Filtro de Bloom com validação teórica e empírica

clc; clear; close all;
fprintf('=== Teste Avançado do Filtro de Bloom (MPEI) ===\n\n');

%% 1. Validação de Falsos Negativos (Propriedade Fundamental)
% Testamos se o filtro garante 100% de acerto para elementos inseridos
m_test = 5000;
k_test = 3;
n_test_keys = 1000;
bf_verify = FiltroBloom(m_test, k_test);

inserted_keys = cell(n_test_keys, 1);
for i = 1:n_test_keys
    inserted_keys{i} = sprintf('VAR_CHR%d_%d_%d', randi(22), randi(1e7), i);
    bf_verify.Inserir(inserted_keys{i});
end

% Verificar todos os inseridos
false_negatives = 0;
for i = 1:n_test_keys
    if ~bf_verify.Verificar(inserted_keys{i})
        false_negatives = false_negatives + 1;
    end
end

fprintf('--- Validação de Corretude ---\n');
fprintf('Chaves inseridas: %d\n', n_test_keys);
fprintf('Falsos Negativos detetados: %d (Esperado: 0)\n', false_negatives);
assert(false_negatives == 0, 'Erro crítico: O Filtro de Bloom gerou um falso negativo!');
fprintf('Sucesso: Propriedade de 0 Falsos Negativos verificada com rigor.\n\n');


%% 2. Curva de Falsos Positivos vs. Tamanho do Filtro (m)
% Avaliamos a taxa de falsos positivos com k fixo e m variável
ks = 3;
ms = 2000:500:10000;            % Amostragem mais densa de m
nTreino = 1000;                 % Número de elementos a inserir
nTeste = 10000;                 % Consultas negativas para teste estatístico

fpRateM = zeros(size(ms));
fpTheoreticalM = zeros(size(ms));

for i = 1:numel(ms)
    bf = FiltroBloom(ms(i), ks);
    
    % Inserir chaves
    dados_treino = cell(nTreino, 1);
    for t = 1:nTreino
        dados_treino{t} = sprintf('VAR_%d', randi(1e7));
        bf.Inserir(dados_treino{t});
    end
    
    % Consultar chaves que sabemos não estarem no filtro
    fp = 0;
    for t = 1:nTeste
        elem = sprintf('TEST_%d', randi(1e7));
        while ismember(elem, dados_treino)
            elem = sprintf('TEST_%d', randi(1e7));
        end
        if bf.Verificar(elem)
            fp = fp + 1;
        end
    end
    
    fpRateM(i) = fp / nTeste;
    fpTheoreticalM(i) = (1 - exp(-ks * nTreino / ms(i)))^ks;
end

% Gráfico 1 – Taxa de FP vs. tamanho m
figure('Position',[100 100 600 400]);
plot(ms, fpRateM * 100, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hold on;
plot(ms, fpTheoreticalM * 100, '--r', 'LineWidth', 1.5);
hold off;
xlabel('Tamanho do vetor m');
ylabel('Taxa de falsos positivos (%)');
title('Taxa de FP vs. Tamanho do Vetor m (k = 3)');
legend('Empírico (Simulado)', 'Teórico', 'Location', 'northeast');
grid on;
saveas(gcf, 'fig_bloom_fp_vs_m_new.png');


%% 3. Curva de Falsos Positivos vs. Número de Hashes (k) - Procura do Ótimo
% Fixamos m = 8000 e n = 1000, variando k para evidenciar o k ótimo
mFixed = 8000;
ksVals = 1:10;                   % Variamos k até 10
fpRateK = zeros(size(ksVals));
fpTheoreticalK = zeros(size(ksVals));

for i = 1:numel(ksVals)
    bf = FiltroBloom(mFixed, ksVals(i));
    
    % Inserir chaves
    dados_treino = cell(nTreino, 1);
    for t = 1:nTreino
        dados_treino{t} = sprintf('VAR_%d', randi(1e7));
        bf.Inserir(dados_treino{t});
    end
    
    % Consultar chaves inexistentes
    fp = 0;
    for t = 1:nTeste
        elem = sprintf('TEST_%d', randi(1e7));
        while ismember(elem, dados_treino)
            elem = sprintf('TEST_%d', randi(1e7));
        end
        if bf.Verificar(elem)
            fp = fp + 1;
        end
    end
    
    fpRateK(i) = fp / nTeste;
    fpTheoreticalK(i) = (1 - exp(-ksVals(i) * nTreino / mFixed))^ksVals(i);
end

% Encontrar o k ótimo teórico e empírico
[minFP_teorico, idx_t] = min(fpTheoreticalK);
[minFP_empirico, idx_e] = min(fpRateK);

fprintf('--- Otimização do Parâmetro k ---\n');
fprintf('Tamanho do filtro m: %d, Elementos inseridos n: %d\n', mFixed, nTreino);
fprintf('k Ótimo Teórico: %d (FP Teórico: %.4f%%)\n', ksVals(idx_t), minFP_teorico * 100);
fprintf('k Ótimo Empírico: %d (FP Empírico: %.4f%%)\n\n', ksVals(idx_e), minFP_empirico * 100);

% Gráfico 2 – Taxa de FP vs. número de funções k
figure('Position',[100 100 600 400]);
plot(ksVals, fpRateK * 100, '-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
hold on;
plot(ksVals, fpTheoreticalK * 100, '--r', 'LineWidth', 1.5);
% Destacar o ponto mínimo
plot(ksVals(idx_e), fpRateK(idx_e) * 100, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
xlabel('Número de funções de hash k');
ylabel('Taxa de falsos positivos (%)');
title('Taxa de FP vs. Número de Funções de Hash k (m = 8000)');
legend('Empírico (Simulado)', 'Teórico', 'k Ótimo Empírico', 'Location', 'northeast');
grid on;
saveas(gcf, 'fig_bloom_fp_vs_k_new.png');
