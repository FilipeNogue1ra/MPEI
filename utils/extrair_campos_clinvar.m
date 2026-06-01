% Extrai Gene e Consequence do campo INFO.
% O campo INFO é uma string com formato "key=value;key=value;...".
% Utilizado pelo modulo de preparação de dados

function [gene, consequence] = extrair_campos_clinvar(infoValue)

    gene = "";
    consequence = "";

    if ismissing(infoValue) || strlength(strtrim(string(infoValue))) == 0
        return;
    end

    entries = split(string(infoValue), ";");
    for i = 1:numel(entries)
        entry = strtrim(entries(i));
        if strlength(entry) == 0
            continue;
        end

        parts = split(entry, "=");
        key = upper(strtrim(parts(1)));

        if key == "GENEINFO" && numel(parts) > 1 && strlength(gene) == 0
            rawGene = strtrim(split(string(strjoin(parts(2:end), "=")), "|"));
            if ~isempty(rawGene)
                gene = normalizeToken(rawGene(1));
            end
        elseif key == "MC" && numel(parts) > 1 && strlength(consequence) == 0
            consequence = normalizeToken(strjoin(parts(2:end), "="));
        elseif key == "CLNDN" && numel(parts) > 1 && strlength(consequence) == 0
            consequence = normalizeToken(strjoin(parts(2:end), "="));
        end
    end
end