% Modulo de MinHash
% Implementação do algoritmo de deteção de mutações similares.
% Divide perfis genéticos/mutações em "shingles" (k-gramas),
% gera matrizes de assinaturas usando funções de hash lineares
% e estima a Similaridade de Jaccard.

classdef MinHash < handle
    
    properties
        numHashes   % Número de funções de hash (tamanho da assinatura)
        shingleSize % Tamanho do shingle (k-gramas, ex: 2 ou 3 caracteres)
        p           % Número primo grande para as funções de hash
        A           % Array de coeficientes 'a' para h(x) = (ax + b) mod p
        B           % Array de coeficientes 'b' para h(x) = (ax + b) mod p
    end
    
    methods
        % 1. Inicializar o objeto MinHash
        function obj = MinHash(n_hashes, k_shingle, seed)
            if nargin < 3
                seed = [];
            end

            if ~isempty(seed)
                rng(seed);
            end

            if nargin < 2
                error('MinHash requer pelo menos o número de hashes e o tamanho do shingle.');
            end

            obj.numHashes = n_hashes;
            obj.shingleSize = k_shingle;
            obj.p = uint64(4294967291); % Número primo grande para evitar colisões
            
            % Gerar coeficientes aleatórios para as funções de dispersão
            % h(x) = (A*x + B) mod p
            limiteSuperior = double(obj.p - 1);
            obj.A = uint64(randi([1, limiteSuperior], 1, n_hashes));
            obj.B = uint64(randi([0, limiteSuperior], 1, n_hashes));
        end
        
        % 2. Extrair Shingles e converter para IDs numéricos
        function ids = obterShingles(obj, str)
            str = char(str); % Garantir que é um vetor de caracteres
            tamanho = length(str);

            if isempty(str)
                ids = uint64(0);
                return;
            end
            
            % Se a string for mais curta que o shingle, usamos a string toda
            if tamanho < obj.shingleSize
                ids = obj.string2id(str);
                return;
            end
            
            numShingles = tamanho - obj.shingleSize + 1;
            ids = zeros(1, numShingles, 'uint64');
            
            % Deslizar a janela para criar os k-gramas
            for i = 1:numShingles
                shingle = str(i:(i + obj.shingleSize - 1));
                % Converter o texto do shingle num número inteiro identificador
                ids(i) = obj.string2id(shingle);
            end
            
            % Remover shingles duplicados na mesma string (Conjunto único)
            ids = unique(ids);
        end
        
        % 3. Gerar a Assinatura MinHash (A Matriz/Vetor Reduzido)
        function assinatura = gerarAssinatura(obj, str)
            % Obter os IDs dos shingles desta string
            shingles_ids = obj.obterShingles(str);
            assinatura = zeros(1, obj.numHashes, 'uint64');
            
            % Para cada função de hash...
            for i = 1:obj.numHashes
                % Aplicar a função h(x) = (a*x + b) mod p a todos os shingles
                valores_hash = mod(obj.A(i) .* shingles_ids + obj.B(i), obj.p);
                
                % O MINHASH: guardar apenas o valor mínimo!
                assinatura(i) = min(valores_hash);
            end
        end
        
        % 4. Estimar a Similaridade de Jaccard entre duas assinaturas
        function similaridade = estimarJaccard(obj, assinatura1, assinatura2)
            % A similaridade é a fração de posições onde as assinaturas coincidem
            coincidencias = sum(assinatura1 == assinatura2);
            similaridade = coincidencias / obj.numHashes;
        end
    end
    
    methods (Access = private)
        % Função auxiliar privada para converter um texto shingle num ID numérico único
        function id = string2id(~, str)
            % Uma adaptação simples do djb2 para gerar inteiros
            id = uint64(5381);
            for i = 1:length(str)
                id = mod(id * 33 + uint64(double(str(i))), uint64(4294967295));
            end
        end
    end
end