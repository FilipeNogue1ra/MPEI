% Modulo de Naive Bayes
% Classifica a categoria de risco de uma mutação

classdef NaiveBayes < handle

    properties
        classes
        featureNames
        classPriors
        featureValues
        likelihoods
        alpha = 1
    end

    methods
        function obj = NaiveBayes(alpha)
            if nargin >= 1 && ~isempty(alpha)
                obj.alpha = alpha;
            end
        end

        function treinar(obj, tabela, targetName, featureNames)
            if nargin < 4
                error('NaiveBayes requer a tabela, o nome da classe-alvo e as features.');
            end

            obj.featureNames = string(featureNames);
            targetName = string(targetName);

            classes = unique(string(tabela.(targetName)));
            obj.classes = classes;
            totalRows = height(tabela);
            obj.classPriors = containers.Map('KeyType', 'char', 'ValueType', 'double');
            obj.featureValues = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.likelihoods = containers.Map('KeyType', 'char', 'ValueType', 'any');

            for c = 1:numel(classes)
                classValue = classes(c);
                classMask = string(tabela.(targetName)) == classValue;
                classTable = tabela(classMask, :);
                obj.classPriors(char(classValue)) = (height(classTable) + obj.alpha) / (totalRows + obj.alpha * numel(classes));

                classLikelihoods = containers.Map('KeyType', 'char', 'ValueType', 'any');
                classFeatureValues = containers.Map('KeyType', 'char', 'ValueType', 'any');

                for f = 1:numel(obj.featureNames)
                    featureName = obj.featureNames(f);
                    values = unique(string(tabela.(featureName)));
                    values(values == "") = [];
                    classFeatureValues(char(featureName)) = values;

                    counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
                    for v = 1:numel(values)
                        value = values(v);
                        count = sum(string(classTable.(featureName)) == value);
                        numerator = count + obj.alpha;
                        denominator = height(classTable) + obj.alpha * numel(values);
                        counts(char(value)) = numerator / denominator;
                    end
                    classLikelihoods(char(featureName)) = counts;
                end

                obj.likelihoods(char(classValue)) = classLikelihoods;
                obj.featureValues(char(classValue)) = classFeatureValues;
            end
        end

        function [classePredita, confianca, pontuacoes] = classificar(obj, amostra)
            if isempty(obj.classes)
                error('NaiveBayes tem de ser treinado antes de classificar.');
            end

            pontuacoes = zeros(1, numel(obj.classes));
            amostra = ensureRowAsStrings(amostra);

            for c = 1:numel(obj.classes)
                classValue = obj.classes(c);
                score = log(obj.classPriors(char(classValue)));
                classLikelihoods = obj.likelihoods(char(classValue));
                classFeatureValues = obj.featureValues(char(classValue));

                for f = 1:numel(obj.featureNames)
                    featureName = obj.featureNames(f);
                    featureMap = classLikelihoods(char(featureName));
                    value = string(amostra.(featureName));
                    if isKey(featureMap, char(value))
                        score = score + log(featureMap(char(value)));
                    else
                        knownValues = classFeatureValues(char(featureName));
                        fallback = obj.alpha / (obj.alpha * max(1, numel(knownValues)) + 1);
                        score = score + log(fallback);
                    end
                end

                pontuacoes(c) = score;
            end

            maxScore = max(pontuacoes);
            expoentes = exp(pontuacoes - maxScore);
            probabilidades = expoentes / sum(expoentes);
            [confianca, idx] = max(probabilidades);
            classePredita = obj.classes(idx);
        end
    end
end

function amostra = ensureRowAsStrings(amostra)
    vars = amostra.Properties.VariableNames;
    for i = 1:numel(vars)
        if ~isstring(amostra.(vars{i}))
            amostra.(vars{i}) = string(amostra.(vars{i}));
        end
    end
end