% Modulo de Naive Bayes
% Classifica a categoria de risco de uma mutação

classdef NaiveBayes < handle

    properties
        classes
        nomes_atributos
        probabilidades_priori
        valores_atributos
        verosimilhancas
        alpha = 1
    end

    methods
        function obj = NaiveBayes(alpha)
            if nargin >= 1 && ~isempty(alpha)
                obj.alpha = alpha;
            end
        end

        function treinar(obj, tabela, nome_alvo, nomes_atributos)
            if nargin < 4
                error('NaiveBayes requer a tabela, o nome da classe-alvo e as features.');
            end

            obj.nomes_atributos = string(nomes_atributos);
            nome_alvo = string(nome_alvo);

            classes_unicas = unique(string(tabela.(nome_alvo)));
            obj.classes = classes_unicas;
            total_linhas = height(tabela);
            obj.probabilidades_priori = containers.Map('KeyType', 'char', 'ValueType', 'double');
            obj.valores_atributos = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.verosimilhancas = containers.Map('KeyType', 'char', 'ValueType', 'any');

            for c = 1:numel(classes_unicas)
                valor_classe = classes_unicas(c);
                mascara_classe = string(tabela.(nome_alvo)) == valor_classe;
                tabela_classe = tabela(mascara_classe, :);
                obj.probabilidades_priori(char(valor_classe)) = (height(tabela_classe) + obj.alpha) / (total_linhas + obj.alpha * numel(classes_unicas));

                verosimilhancas_classe = containers.Map('KeyType', 'char', 'ValueType', 'any');
                valores_atributos_classe = containers.Map('KeyType', 'char', 'ValueType', 'any');

                for f = 1:numel(obj.nomes_atributos)
                    nome_atributo = obj.nomes_atributos(f);
                    valores = unique(string(tabela.(nome_atributo)));
                    valores(valores == "") = [];
                    valores_atributos_classe(char(nome_atributo)) = valores;

                    contagens = containers.Map('KeyType', 'char', 'ValueType', 'double');
                    for v = 1:numel(valores)
                        valor = valores(v);
                        contagem = sum(string(tabela_classe.(nome_atributo)) == valor);
                        numerador = contagem + obj.alpha;
                        denominador = height(tabela_classe) + obj.alpha * numel(valores);
                        contagens(char(valor)) = numerador / denominador;
                    end
                    verosimilhancas_classe(char(nome_atributo)) = contagens;
                end

                obj.verosimilhancas(char(valor_classe)) = verosimilhancas_classe;
                obj.valores_atributos(char(valor_classe)) = valores_atributos_classe;
            end
        end

        function [classe_predita, confianca, pontuacoes] = classificar(obj, amostra)
            if isempty(obj.classes)
                error('NaiveBayes tem de ser treinado antes de classificar.');
            end

            pontuacoes = zeros(1, numel(obj.classes));
            amostra = garantir_colunas_como_strings(amostra);

            for c = 1:numel(obj.classes)
                valor_classe = obj.classes(c);
                pontuacao = log(obj.probabilidades_priori(char(valor_classe)));
                verosimilhancas_classe = obj.verosimilhancas(char(valor_classe));
                valores_atributos_classe = obj.valores_atributos(char(valor_classe));

                for f = 1:numel(obj.nomes_atributos)
                    nome_atributo = obj.nomes_atributos(f);
                    mapa_atributo = verosimilhancas_classe(char(nome_atributo));
                    valor = string(amostra.(nome_atributo));
                    if isKey(mapa_atributo, char(valor))
                        pontuacao = pontuacao + log(mapa_atributo(char(valor)));
                    else
                        valores_conhecidos = valores_atributos_classe(char(nome_atributo));
                        alternativa = obj.alpha / (obj.alpha * max(1, numel(valores_conhecidos)) + 1);
                        pontuacao = pontuacao + log(alternativa);
                    end
                end

                pontuacoes(c) = pontuacao;
            end

            max_pontuacao = max(pontuacoes);
            expoentes = exp(pontuacoes - max_pontuacao);
            probabilidades = expoentes / sum(expoentes);
            [confianca, idx] = max(probabilidades);
            classe_predita = obj.classes(idx);
        end
    end
end

function amostra = garantir_colunas_como_strings(amostra)
    variaveis = amostra.Properties.VariableNames;
    for i = 1:numel(variaveis)
        if ~isstring(amostra.(variaveis{i}))
            amostra.(variaveis{i}) = string(amostra.(variaveis{i}));
        end
    end
end