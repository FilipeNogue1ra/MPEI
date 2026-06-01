% Constrói o perfil textual da variante.
% O perfil segue a estrutura: Gene, Cromossoma, REF, ALT e Consequence.

function variantText = build_clinvar_variant_string(row)

    if ~istable(row) || height(row) ~= 1
        error('build_clinvar_variant_string espera uma tabela com uma única linha.');
    end

    variantText = strjoin(collectTokens(row), " ");
end

function tokens = collectTokens(row)
    tokens = strings(0, 1);
    variableNames = string(row.Properties.VariableNames);

    geneCol = findColumnName(variableNames, ["GENE", "Gene"]);
    chromCol = findColumnName(variableNames, ["CHROM", "Chrom", "chrom"]);
    refCol = findColumnName(variableNames, ["REF", "Ref", "ref"]);
    altCol = findColumnName(variableNames, ["ALT", "Alt", "alt"]);
    infoCol = findColumnName(variableNames, ["INFO", "Info", "info"]);

    gene = "";
    consequence = "";

    if strlength(geneCol) > 0
        gene = normalizeToken(row.(geneCol));
    end

    if strlength(infoCol) > 0
        [geneFromInfo, consequenceFromInfo] = extrair_campos_clinvar(row.(infoCol));
        if strlength(gene) == 0
            gene = geneFromInfo;
        end
        consequence = consequenceFromInfo;
    end

    if strlength(gene) > 0
        tokens(end + 1, 1) = "gene:" + gene; %#ok<AGROW>
    end
    if strlength(chromCol) > 0
        chromValue = normalizeToken(row.(chromCol));
        if strlength(chromValue) > 0
            tokens(end + 1, 1) = "chrom:" + chromValue; %#ok<AGROW>
        end
    end
    if strlength(refCol) > 0
        refValue = normalizeToken(row.(refCol));
        if strlength(refValue) > 0
            tokens(end + 1, 1) = "ref:" + refValue; %#ok<AGROW>
        end
    end
    if strlength(altCol) > 0
        altValue = normalizeToken(row.(altCol));
        if strlength(altValue) > 0
            tokens(end + 1, 1) = "alt:" + altValue; %#ok<AGROW>
        end
    end
    if strlength(consequence) > 0
        tokens(end + 1, 1) = "conseq:" + consequence; %#ok<AGROW>
    end

    tokens = unique(tokens(tokens ~= ""));
end

function fieldName = findColumnName(variableNames, desiredNames)
    fieldName = "";
    variableNames = string(variableNames);
    desiredNames = string(desiredNames);
    for i = 1:numel(desiredNames)
        idx = find(lower(variableNames) == lower(desiredNames(i)), 1);
        if ~isempty(idx)
            fieldName = variableNames(idx);
            return;
        end
    end
end

% normalizeToken now implemented in utils/normalizeToken.m