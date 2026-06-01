% Modulo de Filtro de Bloom
classdef FiltroBloom < handle
    
    properties
        m
        k
        vetor
    end
    
    methods
        function obj = FiltroBloom(tamanho_m, num_hashes_k)
            if nargin < 2
                error('FiltroBloom requer o tamanho do vetor e o número de hashes.');
            end

            if ~isscalar(tamanho_m) || ~isscalar(num_hashes_k) || tamanho_m <= 0 || num_hashes_k <= 0
                error('Os parâmetros do FiltroBloom têm de ser escalares positivos.');
            end

            obj.m = tamanho_m;
            obj.k = num_hashes_k;
            obj.vetor = false(1, tamanho_m);
        end
        
        function Inserir(obj, elemento)
            elemento_char = char(elemento); 
            for i = 1:obj.k
                pos = obj.string2hash(elemento_char, i);
                obj.vetor(pos) = true;
            end
        end
        
        function resultado = Verificar(obj, elemento)
            elemento_char = char(elemento);
            resultado = true;
            
            for i = 1:obj.k
                pos = obj.string2hash(elemento_char, i);
                if obj.vetor(pos) == false
                    resultado = false;
                    return; 
                end
            end
        end
    end
    
    methods (Access = private)
        % djb2 hash function with seed
        function h = string2hash(obj, str, seed)
            h_val = 5381 + seed; 
            for c = 1:length(str)
                h_val = mod(h_val * 33 + double(str(c)), obj.m);
            end
            h = h_val + 1; 
        end
    end
end