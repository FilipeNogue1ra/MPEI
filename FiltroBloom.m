% Modulo de Filtro de Bloom
% Verifica a pertença de uma mutação a um grupo de mutações de risco

classdef FiltroBloom < handle
    
    properties
        m       % Tamanho do vetorr
        k       % Número de funções de hash
        vetor   % Vetor de bits (booleano)
    end
    
    methods
        % 1. Inicializar o vetor
        function obj = FiltroBloom(tamanho_m, num_hashes_k)
            if nargin < 2
                error('FiltroBloom requer o tamanho do vetor e o número de hashes.');
            end

            if ~isscalar(tamanho_m) || ~isscalar(num_hashes_k) || tamanho_m <= 0 || num_hashes_k <= 0
                error('Os parâmetros do FiltroBloom têm de ser escalares positivos.');
            end

            % cria um filtro vazio de tamanho m com k funções
            obj.m = tamanho_m;
            obj.k = num_hashes_k;
            obj.vetor = false(1, tamanho_m); % Vetor de zeros 
        end
        
        % 2. Inserir uma string (ex: "EGFR_L858R")
        function Inserir(obj, elemento)
            % Converte para char para podermos iterar sobre as letras
            elemento_char = char(elemento); 
            
            % Calcula as k posições de hash e coloca a 'true'
            for i = 1:obj.k
                pos = obj.string2hash(elemento_char, i);
                obj.vetor(pos) = true;
            end
        end
        
        % 3. Verificar a pertença de uma string
        function resultado = Verificar(obj, elemento)
            elemento_char = char(elemento);
            resultado = true; % Assumimos que pertence até provar o contrário
            
            % Verifica as k posições
            for i = 1:obj.k
                pos = obj.string2hash(elemento_char, i);
                if obj.vetor(pos) == false
                    resultado = false; % Se = 0, NÃO pertence de certeza
                    return; 
                end
            end
        end
    end
    
    methods (Access = private)
        % Função de dispersão (Hash) baseada no algoritmo djb2 
        function h = string2hash(obj, str, seed)
            % A semente (seed) varia de 1 até k 
            h_val = 5381 + seed; 
            
            for c = 1:length(str)
                % Algoritmo: hash = hash * 33 + char
                h_val = mod(h_val * 33 + double(str(c)), obj.m);
            end
            h = h_val + 1; 
        end
    end
end