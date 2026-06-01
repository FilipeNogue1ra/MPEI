% Constrói o perfil textual da variante.
function texto_variante = build_clinvar_variant_string(linha)

    if ~istable(linha) || height(linha) ~= 1
        error('build_clinvar_variant_string espera uma tabela com uma única linha.');
    end

    texto_variante = strjoin(obter_tokens(linha), " ");
end

function tokens = obter_tokens(linha)
    tokens = strings(0, 1);
    nomes_variaveis = string(linha.Properties.VariableNames);

    col_gene = procurar_nome_coluna(nomes_variaveis, ["GENE", "Gene"]);
    col_chrom = procurar_nome_coluna(nomes_variaveis, ["CHROM", "Chrom", "chrom"]);
    col_ref = procurar_nome_coluna(nomes_variaveis, ["REF", "Ref", "ref"]);
    col_alt = procurar_nome_coluna(nomes_variaveis, ["ALT", "Alt", "alt"]);
    col_info = procurar_nome_coluna(nomes_variaveis, ["INFO", "Info", "info"]);

    gene = "";
    consequencia = "";

    if strlength(col_gene) > 0
        gene = normalizeToken(linha.(col_gene));
    end

    if strlength(col_info) > 0
        [gene_da_info, consequencia_da_info] = extrair_campos_clinvar(linha.(col_info));
        if strlength(gene) == 0
            gene = gene_da_info;
        end
        consequencia = consequencia_da_info;
    end

    if strlength(gene) > 0
        tokens(end + 1, 1) = "gene:" + gene;
    end
    if strlength(col_chrom) > 0
        valor_chrom = normalizeToken(linha.(col_chrom));
        if strlength(valor_chrom) > 0
            tokens(end + 1, 1) = "chrom:" + valor_chrom;
        end
    end
    if strlength(col_ref) > 0
        valor_ref = normalizeToken(linha.(col_ref));
        if strlength(valor_ref) > 0
            tokens(end + 1, 1) = "ref:" + valor_ref;
        end
    end
    if strlength(col_alt) > 0
        valor_alt = normalizeToken(linha.(col_alt));
        if strlength(valor_alt) > 0
            tokens(end + 1, 1) = "alt:" + valor_alt;
        end
    end
    if strlength(consequencia) > 0
        tokens(end + 1, 1) = "conseq:" + consequencia;
    end

    tokens = unique(tokens(tokens ~= ""));
end

function nome_campo = procurar_nome_coluna(nomes_variaveis, nomes_desejados)
    nome_campo = "";
    nomes_variaveis = string(nomes_variaveis);
    nomes_desejados = string(nomes_desejados);
    for i = 1:numel(nomes_desejados)
        idx = find(lower(nomes_variaveis) == lower(nomes_desejados(i)), 1);
        if ~isempty(idx)
            nome_campo = nomes_variaveis(idx);
            return;
        end
    end
end