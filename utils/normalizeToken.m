function token = normalizeToken(valor)
% Normaliza um valor de string para token consistente

    token = lower(strtrim(string(valor)));
    token = replace(token, [" ", ",", ";", "\t"], "_");
    token = regexprep(token, "[^a-z0-9_:\\-\\.\\|]", "_");
    token = regexprep(token, "_+", "_");
    token = strip(token, "_");
end
