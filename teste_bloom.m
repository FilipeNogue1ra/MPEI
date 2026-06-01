% teste_bloom.m
clc; clear; close all;
fprintf('=== Teste Avançado do Filtro de Bloom (MPEI) ===\n\n');

%% 1. Validação de Falsos Negativos
m_teste = 5000;
k_teste = 3;
n_chaves_teste = 1000;
bf_verificacao = FiltroBloom(m_teste, k_teste);

chaves_inseridas = cell(n_chaves_teste, 1);
for i = 1:n_chaves_teste
    chaves_inseridas{i} = sprintf('VAR_CHR%d_%d_%d', randi(22), randi(1e7), i);
    bf_verificacao.Inserir(chaves_inseridas{i});
end

falsos_negativos = 0;
for i = 1:n_chaves_teste
    if ~bf_verificacao.Verificar(chaves_inseridas{i})
        falsos_negativos = falsos_negativos + 1;
    end
end

fprintf('--- Validação de Corretude ---\n');
fprintf('Chaves inseridas: %d\n', n_chaves_teste);
fprintf('Falsos Negativos detetados: %d (Esperado: 0)\n', falsos_negativos);
assert(falsos_negativos == 0, 'Erro crítico: O Filtro de Bloom gerou um falso negativo!');
fprintf('Sucesso: Propriedade de 0 Falsos Negativos verificada com rigor.\n\n');

%% 2. Curva de Falsos Positivos vs. Tamanho do Filtro (m)
ks = 3;
ms = 2000:500:10000;
n_treino = 1000;
n_teste = 10000;

taxa_fp_m = zeros(size(ms));
teorico_fp_m = zeros(size(ms));

for i = 1:numel(ms)
    bf = FiltroBloom(ms(i), ks);
    
    dados_treino = cell(n_treino, 1);
    for t = 1:n_treino
        dados_treino{t} = sprintf('VAR_%d', randi(1e7));
        bf.Inserir(dados_treino{t});
    end
    
    falsos_positivos = 0;
    for t = 1:n_teste
        elem = sprintf('TEST_%d', randi(1e7));
        while ismember(elem, dados_treino)
            elem = sprintf('TEST_%d', randi(1e7));
        end
        if bf.Verificar(elem)
            falsos_positivos = falsos_positivos + 1;
        end
    end
    
    taxa_fp_m(i) = falsos_positivos / n_teste;
    teorico_fp_m(i) = (1 - exp(-ks * n_treino / ms(i)))^ks;
end

figure('Position',[100 100 600 400]);
plot(ms, taxa_fp_m * 100, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hold on;
plot(ms, teorico_fp_m * 100, '--r', 'LineWidth', 1.5);
hold off;
xlabel('Tamanho do vetor m');
ylabel('Taxa de falsos positivos (%)');
title('Taxa de FP vs. Tamanho do Vetor m (k = 3)');
legend('Empírico (Simulado)', 'Teórico', 'Location', 'northeast');
grid on;
saveas(gcf, 'fig_bloom_fp_vs_m_new.png');

%% 3. Curva de Falsos Positivos vs. Número de Hashes (k)
m_fixo = 8000;
valores_k = 1:10;
taxa_fp_k = zeros(size(valores_k));
teorico_fp_k = zeros(size(valores_k));

for i = 1:numel(valores_k)
    bf = FiltroBloom(m_fixo, valores_k(i));
    
    dados_treino = cell(n_treino, 1);
    for t = 1:n_treino
        dados_treino{t} = sprintf('VAR_%d', randi(1e7));
        bf.Inserir(dados_treino{t});
    end
    
    falsos_positivos = 0;
    for t = 1:n_teste
        elem = sprintf('TEST_%d', randi(1e7));
        while ismember(elem, dados_treino)
            elem = sprintf('TEST_%d', randi(1e7));
        end
        if bf.Verificar(elem)
            falsos_positivos = falsos_positivos + 1;
        end
    end
    
    taxa_fp_k(i) = falsos_positivos / n_teste;
    teorico_fp_k(i) = (1 - exp(-valores_k(i) * n_treino / m_fixo))^valores_k(i);
end

[min_fp_teorico, idx_teorico] = min(teorico_fp_k);
[min_fp_empirico, idx_empirico] = min(taxa_fp_k);

fprintf('--- Otimização do Parâmetro k ---\n');
fprintf('Tamanho do filtro m: %d, Elementos inseridos n: %d\n', m_fixo, n_treino);
fprintf('k Ótimo Teórico: %d (FP Teórico: %.4f%%)\n', valores_k(idx_teorico), min_fp_teorico * 100);
fprintf('k Ótimo Empírico: %d (FP Empírico: %.4f%%)\n\n', valores_k(idx_empirico), min_fp_empirico * 100);

figure('Position',[100 100 600 400]);
plot(valores_k, taxa_fp_k * 100, '-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
hold on;
plot(valores_k, teorico_fp_k * 100, '--r', 'LineWidth', 1.5);
plot(valores_k(idx_empirico), taxa_fp_k(idx_empirico) * 100, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
xlabel('Número de funções de hash k');
ylabel('Taxa de falsos positivos (%)');
title('Taxa de FP vs. Número de Funções de Hash k (m = 8000)');
legend('Empírico (Simulado)', 'Teórico', 'k Ótimo Empírico', 'Location', 'northeast');
grid on;
saveas(gcf, 'fig_bloom_fp_vs_k_new.png');
