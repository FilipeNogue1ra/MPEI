% Extrai Gene e Consequence do campo INFO.
function [gene, consequencia] = extrair_campos_clinvar(valor_info)

    gene = "";
    consequencia = "";

    if ismissing(valor_info) || strlength(strtrim(string(valor_info))) == 0
        return;
    end

    entradas = split(string(valor_info), ";");
    for i = 1:numel(entradas)
        entrada = strtrim(entradas(i));
        if strlength(entrada) == 0
            continue;
        end

        partes = split(entrada, "=");
        chave = upper(strtrim(partes(1)));

        if chave == "GENEINFO" && numel(partes) > 1 && strlength(gene) == 0
            gene_bruto = strtrim(split(string(strjoin(partes(2:end), "=")), "|"));
            if ~isempty(gene_bruto)
                gene = normalizeToken(gene_bruto(1));
            end
        elseif chave == "MC" && numel(partes) > 1 && strlength(consequencia) == 0
            consequencia = normalizeToken(strjoin(partes(2:end), "="));
        elseif chave == "CLNDN" && numel(partes) > 1 && strlength(consequencia) == 0
            consequencia = normalizeToken(strjoin(partes(2:end), "="));
        end
    end
end