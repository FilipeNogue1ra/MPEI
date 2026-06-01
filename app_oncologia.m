% app_oncologia.m
% Interface Grafica do MinHash

function app_oncologia()
    addpath('utils');
    
    dados_treino = [];
    dados_teste = [];
    assinaturas_treino = [];
    minhash = [];
    
    fig = uifigure('Name', 'Sistema de Triagem Clinico-Oncologica - MPEI', ...
        'Position', [100 100 800 620]);
    
    % Componentes UI
    pnlControl = uipanel(fig, 'Title', '1. Carregamento de Dados e Treino', 'Position', [20 500 760 100]);
    lblStatus = uilabel(pnlControl, 'Text', 'Estado: Aguardando carregamento...', 'Position', [20 40 400 22], 'FontWeight', 'bold');
    btnLoad = uibutton(pnlControl, 'Text', 'Carregar ClinVar e Treinar MinHash', ...
        'Position', [450 30 280 35], 'ButtonPushedFcn', @(btn, event) carregarETreinar(), ...
        'BackgroundColor', [0.1 0.5 0.8], 'FontColor', [1 1 1], 'FontWeight', 'bold');
        
    pnlPaciente = uipanel(fig, 'Title', '2. Selecionar Paciente de Teste', 'Position', [20 280 760 200]);
    uilabel(pnlPaciente, 'Text', 'Escolha um Paciente do Conjunto de Teste:', 'Position', [20 145 300 22]);
    ddPacientes = uidropdown(pnlPaciente, 'Position', [20 115 720 25], 'Items', {'(Carregue os dados primeiro)'}, 'Enable', 'off', ...
        'ValueChangedFcn', @(dd, event) atualizarDetalhesPaciente());
    lblDetalhes = uilabel(pnlPaciente, 'Text', 'Detalhes do paciente selecionado aparecerao aqui.', 'Position', [20 15 720 85], 'WordWrap', 'on', 'FontAngle', 'italic');
    
    pnlResultados = uipanel(fig, 'Title', '3. Resultados da Comparacao MinHash', 'Position', [20 20 760 240]);
    btnComparar = uibutton(pnlResultados, 'Text', 'Executar Comparacao MinHash', ...
        'Position', [20 175 240 30], 'Enable', 'off', 'ButtonPushedFcn', @(btn, event) executarMinHash(), ...
        'BackgroundColor', [0.1 0.7 0.4], 'FontColor', [1 1 1], 'FontWeight', 'bold');
    tblResultados = uitable(pnlResultados, 'Position', [20 20 720 150]);
    tblResultados.ColumnName = {'Rank', 'Similaridade Jaccard', 'Gene', 'Mutacao (REF>ALT)', 'Significancia Clinica'};
    tblResultados.ColumnWidth = {60, 140, 120, 140, 'auto'};
    
    % Callbacks
    function carregarETreinar()
        btnLoad.Enable = 'off';
        lblStatus.Text = 'A ler clinvar.csv e a preparar dados...';
        drawnow;
        
        dados_mat = 'dados_clinvar_processados.mat';
        try
            if ~isfile(dados_mat)
                [dados_treino, dados_teste, ~] = preparacao_dados_clinvar('clinvar.csv', dados_mat, 5000);
            else
                load(dados_mat, 'dados_treino', 'dados_teste');
            end
            
            lblStatus.Text = 'A gerar assinaturas MinHash...';
            drawnow;
            
            numHashes = 150;
            shingleSize = 5;
            minhash = MinHash(numHashes, shingleSize, 42);
            
            n_treino = height(dados_treino);
            assinaturas_treino = zeros(n_treino, numHashes, 'uint64');
            for i = 1:n_treino
                assinaturas_treino(i, :) = minhash.gerarAssinatura(char(dados_treino.Perfil(i)));
            end
            
            lblStatus.Text = sprintf('Dados Carregados! Treino Concluido (%d amostras)', n_treino);
            
            n_teste = height(dados_teste);
            itens_dropdown = cell(n_teste, 1);
            for i = 1:n_teste
                itens_dropdown{i} = sprintf('Paciente %d: Gene %s (Mut: %s>%s)', ...
                    i, dados_teste.Gene(i), dados_teste.REF(i), dados_teste.ALT(i));
            end
            ddPacientes.Items = itens_dropdown;
            ddPacientes.Enable = 'on';
            btnComparar.Enable = 'on';
            atualizarDetalhesPaciente();
            
        catch ME
            lblStatus.Text = ['Erro: ' ME.message];
            btnLoad.Enable = 'on';
        end
    end

    function atualizarDetalhesPaciente()
        if isempty(dados_teste), return; end
        
        idx_selecionado = ddPacientes.Value;
        idx_partes = sscanf(idx_selecionado, 'Paciente %d:');
        idx_paciente = idx_partes(1);
        
        paciente = dados_teste(idx_paciente, :);
        lblDetalhes.Text = sprintf('<b>Gene:</b> %s | <b>Cromossoma:</b> %s | <b>Mutacao:</b> %s > %s | <b>Consequencia:</b> %s<br><b>Significancia Clinica Real:</b> %s<br><b>Perfil:</b> %s', ...
            paciente.Gene, paciente.CHROM, paciente.REF, paciente.ALT, paciente.Consequence, paciente.CLNSIG, paciente.Perfil);
        lblDetalhes.Interpreter = 'html'; 
    end

    function executarMinHash()
        if isempty(dados_teste) || isempty(minhash), return; end
        
        idx_selecionado = ddPacientes.Value;
        idx_partes = sscanf(idx_selecionado, 'Paciente %d:');
        idx_paciente = idx_partes(1);
        
        paciente = dados_teste(idx_paciente, :);
        perfil_novo = char(paciente.Perfil);
        
        assinatura_novo = minhash.gerarAssinatura(perfil_novo);
        
        n_treino = height(dados_treino);
        similaridades = zeros(n_treino, 1);
        for j = 1:n_treino
            similaridades(j) = minhash.estimarJaccard(assinatura_novo, assinaturas_treino(j, :));
        end
        
        [valores_sim, indices_top] = maxk(similaridades, 3);
        
        data = cell(3, 5);
        for k = 1:3
            idx_hist = indices_top(k);
            paciente_hist = dados_treino(idx_hist, :);
            
            data{k, 1} = sprintf('Top %d', k);
            data{k, 2} = sprintf('%.1f%%', valores_sim(k) * 100);
            data{k, 3} = char(paciente_hist.Gene);
            data{k, 4} = sprintf('%s > %s', paciente_hist.REF, paciente_hist.ALT);
            data{k, 5} = char(paciente_hist.CLNSIG);
        end
        tblResultados.Data = data;
    end
end
