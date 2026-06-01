% teste_minhash.m
% Teste de validacao do MinHash vs Jaccard Real

clc; clear; close all;

fprintf('=== Teste Isolado do Algoritmo MinHash ===\n\n');

str1 = 'Gene:SAMD11|Chrom:1|Ref:G|Alt:A|Pathogenic';
str2 = 'Gene:SAMD11|Chrom:1|Ref:G|Alt:T|Benign';

k_shingle = 3;
n_hashes = 200;

minhash = MinHash(n_hashes, k_shingle, 42);

shingles_1 = minhash.obterShingles(str1);
shingles_2 = minhash.obterShingles(str2);

jaccard_real = length(intersect(shingles_1, shingles_2)) / length(union(shingles_1, shingles_2));

assinatura1 = minhash.gerarAssinatura(str1);
assinatura2 = minhash.gerarAssinatura(str2);

jaccard_estimado = minhash.estimarJaccard(assinatura1, assinatura2);

fprintf('String 1: %s\n', str1);
fprintf('String 2: %s\n', str2);
fprintf('\n');
fprintf('Numero de Hashes: %d\n', n_hashes);
fprintf('Tamanho do Shingle: %d\n', k_shingle);
fprintf('\n');
fprintf('Jaccard Real/Exato: %.4f\n', jaccard_real);
fprintf('Jaccard Estimado (MinHash): %.4f\n', jaccard_estimado);
fprintf('Erro Absoluto: %.4f\n', abs(jaccard_real - jaccard_estimado));
