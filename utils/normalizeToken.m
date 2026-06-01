function token = normalizeToken(value)
% normalizeToken Normaliza um valor textual para token consistente
%   Convém usar esta utilidade a partir de vários módulos para garantir
%   comportamento idêntico ao transformar campos em tokens.

    token = lower(strtrim(string(value)));
    token = replace(token, [" ", ",", ";", "\t"], "_");
    token = regexprep(token, "[^a-z0-9_:\\-\\.\\|]", "_");
    token = regexprep(token, "_+", "_");
    token = strip(token, "_");
end
