% Modulo de MinHash
classdef MinHash < handle
    
    properties
        num_hashes
        tamanho_shingle
        p
        A
        B
    end
    
    methods
        function obj = MinHash(n_hashes, k_shingle, semente)
            if nargin < 3
                semente = [];
            end

            if ~isempty(semente)
                rng(semente);
            end

            if nargin < 2
                error('MinHash requer pelo menos o número de hashes e o tamanho do shingle.');
            end

            obj.num_hashes = n_hashes;
            obj.tamanho_shingle = k_shingle;
            obj.p = uint64(4294967291);
            
            limite_superior = double(obj.p - 1);
            obj.A = uint64(randi([1, limite_superior], 1, n_hashes));
            obj.B = uint64(randi([0, limite_superior], 1, n_hashes));
        end
        
        function ids = obterShingles(obj, texto)
            texto = char(texto);
            tamanho = length(texto);

            if isempty(texto)
                ids = uint64(0);
                return;
            end
            
            if tamanho < obj.tamanho_shingle
                ids = obj.texto2id(texto);
                return;
            end
            
            num_shingles = tamanho - obj.tamanho_shingle + 1;
            ids = zeros(1, num_shingles, 'uint64');
            
            for i = 1:num_shingles
                shingle = texto(i:(i + obj.tamanho_shingle - 1));
                ids(i) = obj.texto2id(shingle);
            end
            
            ids = unique(ids);
        end
        
        function assinatura = gerarAssinatura(obj, texto)
            ids_shingles = obj.obterShingles(texto);
            assinatura = zeros(1, obj.num_hashes, 'uint64');
            
            for i = 1:obj.num_hashes
                valores_hash = mod(obj.A(i) .* ids_shingles + obj.B(i), obj.p);
                assinatura(i) = min(valores_hash);
            end
        end
        
        function similaridade = estimarJaccard(obj, assinatura1, assinatura2)
            coincidencias = sum(assinatura1 == assinatura2);
            similaridade = coincidencias / obj.num_hashes;
        end
    end
    
    methods (Access = private)
        % Adaptação do djb2 para gerar IDs inteiros para os shingles
        function id = texto2id(~, texto)
            id = uint64(5381);
            for i = 1:length(texto)
                id = mod(id * 33 + uint64(double(texto(i))), uint64(4294967295));
            end
        end
    end
end