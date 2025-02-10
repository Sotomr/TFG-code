function EQMOmod()
    %=============================================================
    % SCRIPT PRINCIPAL: EQMOmod
    % Se ejecuta la optimización multiobjetivo con la lógica de EQMO
    % y se muestran los gráficos del cronograma y la curva de satisfacción.
    %=============================================================
    fprintf('>> EQMOmod\n');
    rng(145763);  % Semilla para reproducibilidad

    % Parámetros de entrada (puedes ajustar estos valores)
    porcentajeOpt    = 40;   % [0..100] (porcentaje de optimización)
    porcentajeGroup  = 50;   % [0..100] (porcentaje para agrupamiento)
    numIteracions    = 300;  % Número máximo de iteraciones
    maxNoImprovement = 10;   % Número máximo de iteraciones sin mejora

    % Llamada al algoritmo mejorado (EQMO)
    [costTotalGeneral, satisfFinal, bestIter] = EQMO(porcentajeOpt, porcentajeGroup, numIteracions, maxNoImprovement);
    
    % Recuperar las variables generadas dentro de EQMO (guardadas en la base)
    EqFinal         = evalin('base', 'EqFinal');
    tEntregaFinal   = evalin('base', 'tEntregaFinal');
    horaPref        = evalin('base', 'horaPref');
    dist            = evalin('base', 'dist');
    vyD             = evalin('base', 'vyD');
    numD            = evalin('base', 'numD');
    numP            = evalin('base', 'numP');
    K               = evalin('base', 'K');
    CapBat0         = evalin('base', 'CapBat0');
    baseConsumD     = evalin('base', 'baseConsumD');
    TRecarregaMax   = evalin('base', 'TRecarregaMax');
    UmbralRecarrega = evalin('base', 'UmbralRecarrega');
    preuD           = evalin('base', 'preuD');
    distMax         = evalin('base', 'distMax');
    costMant0       = evalin('base', 'costMant0');
    costE           = evalin('base', 'costE');
    costAss         = evalin('base', 'costAss');
    pesos           = evalin('base', 'pesos');
    ncicles         = evalin('base', 'ncicles');
    alfa            = evalin('base', 'alfa');
    costBat         = evalin('base', 'costBat');
    tvol            = evalin('base', 'tvol');
    mD              = evalin('base', 'mD');
    horesVolAny     = evalin('base', 'horesVolAny');
    costFixDron     = evalin('base', 'costFixDron');
    costExternalitzacio = evalin('base', 'costExternalitzacio');
    
    % Mostrar resultados en la ventana de comandos
    fprintf('------------------------------------------\n');
    fprintf('Cost Total General: %.2f euros\n', costTotalGeneral);
    fprintf('Satisfacción Final  : %.2f %%\n', satisfFinal);
    fprintf('Mejor iteración     : %d\n', bestIter);
    fprintf('------------------------------------------\n');
    
    %---------------------------------------------------------
    % Graficar el cronograma de entregas y la curva de satisfacción
    %---------------------------------------------------------
    deltaT_min = 5;  % Discretización en minutos para el cronograma
    ConstruirCronograma(EqFinal, dist, vyD, numD, numP, K, deltaT_min, CapBat0, baseConsumD, TRecarregaMax, UmbralRecarrega);
    GraficarSatisfaccio(EqFinal, tEntregaFinal, horaPref, K);
    
    % Mostrar resultados finales en pantalla
    MostrarResultatsFinals(EqFinal, tEntregaFinal, bestIter, ...
        horaPref, numP, numD, preuD, CapBat0, distMax, ...
        costMant0, costE, costAss, pesos, dist, ...
        ncicles, alfa, costBat, tvol, mD, vyD, ...
        horesVolAny, K, costFixDron, costExternalitzacio);
end

%==========================================================================
% FUNCIÓN EQMO (Lógica mejorada de optimización multiobjetivo)
%==========================================================================
function [costTotalGeneral, satisfFinal, bestIter] = EQMO(porcentajeOpt, porcentajeGroup, numIteracions, maxNoImprovement)
    fprintf('>> DistribucioMultiobjectiu4c_Modified2\n');
    rng(145763);
    
    numP = 50;  
    numD = 15;
    K    = 12;   % Horas totales de la jornada
    
    preuD = round(linspace(3000, 10000, numD));
    CapBat0     = zeros(1, numD);
    vidaUtilD   = zeros(1, numD);
    distMax     = zeros(1, numD);
    costMant0   = zeros(1, numD);
    costE       = zeros(1, numD);
    costAss     = zeros(1, numD);
    baseConsumD = zeros(1, numD);
    
    for d = 1:numD
        if preuD(d) < 5000
            CapBat0(d)   = round(30 + (800 - 30)*rand);
            vidaUtilD(d) = round(300 + (600 - 300)*rand);
            distMax(d)   = round(30 + (50 - 30)*rand);
            costMant0(d) = round(2 + (4 - 2)*rand);
            costE(d)     = 0.0002 + (0.0003 - 0.0002)*rand; 
            costAss(d)   = round(200 + (600 - 200)*rand);
            baseConsumD(d)= 15; 
        elseif preuD(d) < 7500
            CapBat0(d)   = round(800 + (1200 - 800)*rand);
            vidaUtilD(d) = round(600 + (800 - 600)*rand);
            distMax(d)   = round(50 + (70 - 50)*rand);
            costMant0(d) = round(3 + (5 - 3)*rand);
            costE(d)     = 0.0002 + (0.00025 - 0.0002)*rand;
            costAss(d)   = round(400 + (700 - 400)*rand);
            baseConsumD(d)= 12;
        else
            CapBat0(d)   = round(1200 + (2000 - 1200)*rand);
            vidaUtilD(d) = round(800 + (1000 - 800)*rand);
            distMax(d)   = round(70 + (90 - 70)*rand);
            costMant0(d) = round(4 + (6 - 4)*rand);
            costE(d)     = 0.00018 + (0.00022 - 0.00018)*rand;
            costAss(d)   = round(600 + (800 - 600)*rand);
            baseConsumD(d)= 10;
        end
    end

    horaPref = randi([0, K*60], 1, numP);
    pesos = round(1 + (20 - 1)*rand(1, numP));
    dist  = round(5 + (10 - 5)*rand(1, numP));

    ncicles = round(50 + (150 - 50)*rand(1, numD));
    alfa    = 0.01 + (0.015 - 0.01)*rand(1, numD);
    costBat = 0.05 + (0.15 - 0.05)*rand(1, numD);
    tvol    = zeros(1, numD);
    mD      = round(5 + (10 - 5)*rand(1, numD));
    vyD     = round(15 + (25 - 15)*rand(1, numD));
    horesVolAny = round(150 + (300 - 150)*rand(1, numD));
    BateriaRestant = CapBat0;

    costFixDron         = 0.5;
    costExternalitzacio = 1;
    preuE               = 0.2;
    TRecarregaMax       = 0.5;
    UmbralRecarrega     = 0.2;
    PesFactor           = 0.1;
    
    % --- MODIFICACIÓN: Ajuste en la lógica de agrupamiento ---
    if porcentajeGroup < 0
        grupTamany = 1;          % Cada paquete en grupo individual
    elseif porcentajeGroup > 100
        grupTamany = numP;       % Todos los paquetes en un solo grupo
    else
        grupTamany = max(1, round(numP * ((100 - porcentajeGroup) / 100)));
    end

    maxEstq         = 30;   
    maxEstq_Satisf  = 100;

    [EqFinal, tEntregaFinal, bestIter] = MultiObjetivo( ...
        porcentajeOpt, porcentajeGroup, numIteracions, maxNoImprovement, ...
        numP, numD, K, ...
        preuD, CapBat0, vidaUtilD, distMax, costMant0, costE, ...
        costAss, baseConsumD, pesos, dist, ...
        ncicles, alfa, costBat, tvol, mD, vyD, ...
        horesVolAny, BateriaRestant, costFixDron, ...
        costExternalitzacio, preuE, TRecarregaMax, ...
        UmbralRecarrega, PesFactor, grupTamany, ...
        maxEstq, maxEstq_Satisf, horaPref);

    S = zeros(1, numP);
    for p = 1:numP
        if ~isnan(tEntregaFinal(p))
            S(p) = CalcularSatisfaccioExacta(tEntregaFinal(p), horaPref(p), K);
        else
            S(p) = 0;
        end
    end

    costTRec_dummy = zeros(1,numD); 
    [~, ~, costTotalGeneralLocal] = CostFinal( ...
        EqFinal, preuD, vidaUtilD, CapBat0, distMax, costMant0, costE, ...
        pesos, dist, ncicles, alfa, costBat, tvol, mD, vyD, ...
        horesVolAny, costAss, K, find(all(EqFinal == 0,1)), ...
        costFixDron, costExternalitzacio, costTRec_dummy);

    paquetsAmbSatisfaccio = any(EqFinal == 1, 1);
    S_ambPaquets = S(paquetsAmbSatisfaccio);
    if isempty(S_ambPaquets)
        satisfFinal = 0;
    else
        satisfFinal = mean(S_ambPaquets);
    end

    costTotalGeneral = costTotalGeneralLocal;
    
    %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % Asignar variables clave a la base para poder usarlas en los gráficos
    assignin('base','EqFinal', EqFinal);
    assignin('base','tEntregaFinal', tEntregaFinal);
    assignin('base','horaPref', horaPref);
    assignin('base','dist', dist);
    assignin('base','vyD', vyD);
    assignin('base','numD', numD);
    assignin('base','numP', numP);
    assignin('base','K', K);
    assignin('base','CapBat0', CapBat0);
    assignin('base','baseConsumD', baseConsumD);
    assignin('base','TRecarregaMax', TRecarregaMax);
    assignin('base','UmbralRecarrega', UmbralRecarrega);
    assignin('base','preuD', preuD);
    assignin('base','distMax', distMax);
    assignin('base','costMant0', costMant0);
    assignin('base','costE', costE);
    assignin('base','costAss', costAss);
    assignin('base','pesos', pesos);
    assignin('base','ncicles', ncicles);
    assignin('base','alfa', alfa);
    assignin('base','costBat', costBat);
    assignin('base','tvol', tvol);
    assignin('base','mD', mD);
    assignin('base','horesVolAny', horesVolAny);
    assignin('base','costFixDron', costFixDron);
    assignin('base','costExternalitzacio', costExternalitzacio);
    %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end

%==========================================================================
function [EqFinal, tEntregaFinal, bestIter] = MultiObjetivo( ...
    porcentajeOpt, porcentajeGroup, numIteracions, maxNoImprovement, ...
    numP, numD, K, ...
    preuD, CapBat0, vidaUtilD, distMax, costMant0, costE, ...
    costAss, baseConsumD, pesos, dist, ...
    ncicles, alfa, costBat, tvol, mD, vyD, ...
    horesVolAny, BateriaRestant, costFixDron, ...
    costExternalitzacio, preuE, TRecarregaMax, ...
    UmbralRecarrega, PesFactor, grupTamany, ...
    maxEstq, maxEstq_Satisf, horaPref)
    
    costMaxSim = numP*costExternalitzacio*sum(dist);  % Límite máximo para externalización
    
    % Inicialización
    Eq = zeros(numD, numP);
    tEntrega = nan(1, numP);
    
    bestWeightedSum = -realmax; 
    bestEq = Eq;
    best_tEntrega = tEntrega;
    bestIter = 0;
    noImprovementCounter = 0;
    
    % División en grupos
    gruposTotal = ceil(numP / grupTamany);
    numGruposSatisf = round((porcentajeOpt / 100) * gruposTotal);
    numGruposCoste = gruposTotal - numGruposSatisf;
    
    grupos = {};
    ini = 1;
    while ini <= numP
        fin = min(ini + grupTamany - 1, numP);
        grupos{end+1} = ini:fin; 
        ini = fin + 1;
    end
    numGrupos = length(grupos);
    assert(numGrupos == gruposTotal, 'Error en la división de grupos.');
    
    % Aleatoriedad
    indicesGrupos = randperm(numGrupos);
    gruposSatisf = indicesGrupos(1:numGruposSatisf);
    gruposCoste  = indicesGrupos(numGruposSatisf+1:end);
    asignacionGrupos = zeros(1, numGrupos);
    asignacionGrupos(gruposSatisf) = 1; 
    
    % Historial para graficar
    costHistory = zeros(1, numIteracions);
    satisfHistory = zeros(1, numIteracions);
    
    for it = 1:numIteracions
        fprintf('\nITERACIÓN #%d de %d\n', it, numIteracions);
        
        if it > 1 && mod(it, ceil(numIteracions / 10)) == 0
            numReasignar = max(1, round(0.1 * numGrupos));
            gruposAReasignar = randperm(numGrupos, numReasignar);
            for g = gruposAReasignar
                if asignacionGrupos(g) == 1 && numGruposSatisf > 0
                    asignacionGrupos(g) = 0;  % cambiar a Coste
                    numGruposSatisf = numGruposSatisf - 1;
                    numGruposCoste = numGruposCoste + 1;
                elseif asignacionGrupos(g) == 0 && numGruposSatisf < gruposTotal
                    asignacionGrupos(g) = 1;  % cambiar a Satisf
                    numGruposSatisf = numGruposSatisf + 1;
                    numGruposCoste = numGruposCoste - 1;
                end
            end
            fprintf('  -> Reasignación de %d grupos.\n', numReasignar);
        end
        
        for g = 1:numGrupos
            paqsG = grupos{g};
            
            if asignacionGrupos(g) == 1
                % Satisfacción
                [EqGrupo, ~, ~, tEntregaGrupo] = NashOptimizado( ...
                    numP, numD, CapBat0, distMax, pesos, dist, tvol, ...
                    BateriaRestant, TRecarregaMax, UmbralRecarrega, ...
                    PesFactor, baseConsumD, maxEstq_Satisf, K, ...
                    vyD, horaPref, paqsG);
                
                for p = paqsG
                    Eq(:, p) = 0;
                    if any(EqGrupo(:, p) == 1)
                        dAsig = find(EqGrupo(:, p) == 1);
                        Eq(dAsig, p) = 1;
                        tEntrega(p) = tEntregaGrupo(p);
                    end
                end
                
            else
                % Coste
                [EqGrupo, ~, ~, tEntregaGrupo] = NashDrons( ...
                    numP, numD, Eq, preuD, vidaUtilD, CapBat0, distMax, ...
                    costMant0, costE, pesos, dist, ncicles, alfa, costBat, tvol, mD, vyD, ...
                    horesVolAny, costAss, K, costFixDron, costExternalitzacio, preuE, ...
                    grupTamany, maxEstq, BateriaRestant, TRecarregaMax, UmbralRecarrega, ...
                    PesFactor, baseConsumD, horaPref, paqsG);
                
                for p = paqsG
                    Eq(:, p) = 0;
                    if any(EqGrupo(:, p) == 1)
                        dAsig = find(EqGrupo(:, p) == 1);
                        Eq(dAsig, p) = 1;
                        tEntrega(p) = tEntregaGrupo(p);
                    end
                end
            end
        end
        
        [~, ~, cTotal] = CostFinal(Eq, preuD, vidaUtilD, CapBat0, ...
            distMax, costMant0, costE, pesos, dist, ncicles, alfa, ...
            costBat, tvol, mD, vyD, horesVolAny, costAss, K, ...
            find(all(Eq == 0,1)), costFixDron, costExternalitzacio, zeros(1, numD));
        
        % Satisfacción media
        SIter = zeros(1, numP);
        for p = 1:numP
            if ~isnan(tEntrega(p))
                SIter(p) = CalcularSatisfaccioExacta(tEntrega(p), horaPref(p), K);
            end
        end
        satIter = mean(SIter);
        
        costHistory(it) = cTotal;
        satisfHistory(it) = satIter;
     
        satNorm = satIter / 100;   % la satisfacción varía de 0 a 100
        costNorm = min(1, cTotal / costMaxSim); 
        alpha = (porcentajeOpt / 100);  % peso de la satisfacción
        
        % Maximizar la combinación => más satisfacción y menor coste
        weightedSum = alpha * satNorm + (1 - alpha) * (1 - costNorm);
        
        fprintf('Fin iter %d: Coste= %.2f, Satisf= %.2f, WeightedSum= %.4f\n', ...
            it, cTotal, satIter, weightedSum);
        
        if (weightedSum > bestWeightedSum)
            bestWeightedSum = weightedSum;
            bestEq = Eq;
            best_tEntrega = tEntrega;
            bestIter = it;
            noImprovementCounter = 0;
            fprintf('Mejor asignación encontrada en la iteración %d.\n', it);
        else
            noImprovementCounter = noImprovementCounter + 1;
            fprintf('Sin mejora significativa. Contador: %d\n', noImprovementCounter);
        end
        
        if noImprovementCounter >= maxNoImprovement
            fprintf('Convergencia alcanzada en la iteración %d.\n', it);
            break;
        end
    end
    
    % Mejor solución encontrada
    EqFinal = bestEq;
    tEntregaFinal = best_tEntrega;
    
    assignin('base', 'costHistory', costHistory(1:bestIter));
    assignin('base', 'satisfHistory', satisfHistory(1:bestIter));
end

%==========================================================================
function [EqGrupo, iter, costTRec, tEntregaGrupo] = NashDrons(...
    numP, numD, Eq, preuD, vidaUtilD, CapBat0, distMax, ...
    costMant0, costE, pesos, dist, ncicles, alfa, costBat, ...
    tvol, mD, vyD, horesVolAny, costAss, K, costFixDron, ...
    costExternalitzacio, preuE, grupTamany, maxEstq, ...
    BateriaRestant, TRecarregaMax, UmbralRecarrega, PesFactor, ...
    baseConsumD, horaPref, paqsAOptimizar)
    
    tEntregaGrupo = nan(1, numP);
    canvi = true;
    iter  = 0;
    bestCost = realmax;
    EqBest = Eq;
    costTRec = zeros(1,numD);
    
    while canvi && iter < length(paqsAOptimizar) * 2  
        canvi = false;
        iter = iter + 1;
        
        for p = paqsAOptimizar
            dronAct = find(Eq(:, p) == 1);
            if isempty(dronAct)
                dronAct = -1;
            end
            millorCost = bestCost;
            millorDron = dronAct;
            
            for d = 1:numD
                if dist(p) > distMax(d)
                    continue;  
                end
                
                % Calcular coste incremental
                distTotal = 2 * dist(p);
                tiempoVueloNeces = distTotal / vyD(d);
                consumoEnergia = baseConsumD(d) * distTotal * (1 + PesFactor * pesos(p));
                
                if BateriaRestant(d) < consumoEnergia
                    faltaEnergia = consumoEnergia - BateriaRestant(d);
                    tiempoRecarga = (faltaEnergia / CapBat0(d)) * TRecarregaMax;
                else
                    tiempoRecarga = 0;
                end
                
                % Restricciones
                if (tvol(d) + tiempoRecarga + tiempoVueloNeces) > K
                    continue;  
                end
                
                EqSim = Eq;
                EqSim(:, p) = 0;
                EqSim(d, p) = 1;
                
                % Calcular coste
                [~, ~, cTotalSim] = CostFinal(EqSim, preuD, vidaUtilD, ...
                    CapBat0, distMax, costMant0, costE, pesos, dist, ...
                    ncicles, alfa, costBat, tvol, mD, vyD, ...
                    horesVolAny, costAss, K, find(all(EqSim == 0,1)), ...
                    costFixDron, costExternalitzacio, costTRec);
                
                if cTotalSim < millorCost
                    millorCost = cTotalSim;
                    millorDron = d;
                end
            end
            
            if millorDron ~= dronAct && millorDron ~= -1
                Eq(:, p) = 0;
                Eq(millorDron, p) = 1;
                tEntregaGrupo(p) = (tvol(millorDron) + tiempoVueloNeces + tiempoRecarga) * 60;
                
                % Actualizar estado del dron
                if BateriaRestant(millorDron) < consumoEnergia
                    BateriaRestant(millorDron) = CapBat0(millorDron);
                    tvol(millorDron) = tvol(millorDron) + tiempoRecarga;
                end
                
                tvol(millorDron) = tvol(millorDron) + tiempoVueloNeces;
                BateriaRestant(millorDron) = BateriaRestant(millorDron) - consumoEnergia;
                
                canvi = true;
            end
        end
        
        % Calcular coste total actual
        [cVec, cExt, cTotal] = CostFinal(...
            Eq, preuD, vidaUtilD, CapBat0, distMax, costMant0, costE, ...
            pesos, dist, ncicles, alfa, costBat, tvol, mD, vyD, ...
            horesVolAny, costAss, K, find(all(Eq == 0,1)), ...
            costFixDron, costExternalitzacio, costTRec);
        
        if cTotal < bestCost
            bestCost = cTotal;
            EqBest = Eq;
        end
    end
    
    EqGrupo = EqBest;
end

%==========================================================================
function [EqGrupo, satisfMedia, iter, tEntregaGrupo] = NashOptimizado(...
    numP, numD, CapBat0, distMax, pesos, dist, tvol, BateriaRestant, ...
    TRecarregaMax, UmbralRecarrega, PesFactor, baseConsumD, ...
    maxEstq_Satisf, K, vyD, horaPref, paqsAOptimizar)
    
    EqGrupo = zeros(numD, numP);
    tEntregaGrupo = nan(1, numP);
    satisfMedia = 0;
    iter = 0;
    
    canvi = true;
    maxIter = 50;
    
    while canvi && iter < maxIter
        canvi = false;
        iter = iter + 1;
        
        for p = paqsAOptimizar
            if ~isnan(tEntregaGrupo(p))
                continue;  
            end
            
            millorSatisf = -inf;
            millorDron   = -1;
            millorTEntrega = inf;
            
            for d = 1:numD
                if dist(p) > distMax(d)
                    continue;  
                end
                
                distTotal = 2 * dist(p);
                tiempoVueloNeces = distTotal / vyD(d);
                consumoEnergia = baseConsumD(d) * distTotal * (1 + PesFactor * pesos(p));
                
                if BateriaRestant(d) < consumoEnergia
                    faltaEnergia = consumoEnergia - BateriaRestant(d);
                    tiempoRecarga = (faltaEnergia / CapBat0(d)) * TRecarregaMax;
                else
                    tiempoRecarga = 0;
                end
                
                if (tvol(d) + tiempoRecarga + tiempoVueloNeces) > K
                    continue;  
                end
                
                tEntregaPossible = (tvol(d) + tiempoRecarga + tiempoVueloNeces) * 60;
                tEntregaPossible = max(tEntregaPossible, horaPref(p));
                
                satisf_p = CalcularSatisfaccioExacta(tEntregaPossible, horaPref(p), K);
                
                if satisf_p > millorSatisf
                    millorSatisf = satisf_p;
                    millorDron = d;
                    millorTEntrega = tEntregaPossible;
                end
            end
            
            if millorDron ~= -1
                EqGrupo(:, p) = 0;
                EqGrupo(millorDron, p) = 1;
                tEntregaGrupo(p) = millorTEntrega;
                
                if BateriaRestant(millorDron) < (baseConsumD(millorDron) * 2 * dist(p) * (1 + PesFactor * pesos(p)))
                    % Recargar
                    BateriaRestant(millorDron) = CapBat0(millorDron);
                    tvol(millorDron) = tvol(millorDron) + tiempoRecarga;
                end
                
                tvol(millorDron) = tvol(millorDron) + tiempoVueloNeces;
                BateriaRestant(millorDron) = BateriaRestant(millorDron) - (baseConsumD(millorDron) * 2 * dist(p) * (1 + PesFactor * pesos(p)));
                
                canvi = true;
            end
        end
    end
    
    Stemp = zeros(1, numP);
    for p = paqsAOptimizar
        if ~isnan(tEntregaGrupo(p))
            Stemp(p) = CalcularSatisfaccioExacta(tEntregaGrupo(p), horaPref(p), K);
        end
    end
    satisfMedia = mean(Stemp);
end

%==========================================================================
function Sval = CalcularSatisfaccioExacta(tEntrega, tPref, tLimHores)
    difH = (tEntrega - tPref)/60;  
    epsilon = 1e-4;
    num = -(tLimHores^2);
    den = 2*log(epsilon/100);
    sigma = sqrt(num/den);

    Sval = 100*exp(-(difH^2)/(2*sigma^2));
    if abs(difH) >= tLimHores
        Sval = 0;
    end
end

%==========================================================================
function [cVec, cExt, cTotal] = CostFinal(Eq, preuD, vidaUtilD, ...
    CapBat0, distMax, costMant0, costE, pesos, dist, ...
    ncicles, alfa, costBat, tvol, mD, vyD, ...
    horesVolAny, costAss, K, paqsNoAsig, costFixDron, ...
    costExternalitzacio, costTRec)
    
    numD = size(Eq,1);
    costosTotals = zeros(1,numD);
    costExternalitzacioTotal = 0;

    for d=1:numD
        paqs = find(Eq(d,:)==1);
        if ~isempty(paqs)
            CDep = tvol(d)*(preuD(d)/max(vidaUtilD(d),1));
            CapBatReal = CapBat0(d)/1000;
            desgaste = 1 - min(alfa(d)*log(max(ncicles(d),1)), 0.5);
            eficiencia = min(1, tvol(d)/K);
            CBat = CapBatReal*desgaste*costBat(d)*eficiencia;
            CMant = costMant0(d)*(1+0.05*tvol(d))*length(paqs);
            CPes=0;
            for pp=paqs
                if dist(pp)>0 && dist(pp)<=distMax(d)
                    CPes = CPes + (mD(d)*9.81*vyD(d)*(dist(pp)/60)) * costE(d)*pesos(pp)*2;
                end
            end
            CAss = (costAss(d)/horesVolAny(d))*length(paqs)*1.5;
            costosTotals(d) = CDep + CBat + CMant + CPes + CAss + costFixDron + costTRec(d);
        else
            costosTotals(d) = costFixDron + costMant0(d)*1;
        end
    end

    if ~isempty(paqsNoAsig)
        for p=paqsNoAsig
            costExternalitzacioTotal = costExternalitzacioTotal + ...
                                       costExternalitzacio * dist(p);
        end
    end

    costTotalGeneral = sum(costosTotals) + costExternalitzacioTotal;
    cVec = sum(costosTotals);
    cExt = costExternalitzacioTotal;
    cTotal = costTotalGeneral;
end

%==========================================================================
% FUNCIONES PARA CONSTRUIR LOS GRÁFICOS (CRONOGRAMA Y CURVA DE SATISFACCIÓN)
%==========================================================================

function ConstruirCronograma(Eq, dist, vyD, numD, numP, K_horas, deltaT_min, CapBat0, baseConsumD, TRecarregaMax, UmbralRecarrega)
    % Construye el cronograma de entregas a partir de la asignación (Eq)
    % Las entregas comienzan a las 7:00 a.m.
    startTime = datetime('today','Format','dd-MMM-yyyy HH:mm') + hours(7);
    startTimeNum = datenum(startTime);
    tiempoTotal_min = K_horas * 60;
    numSlots = floor(tiempoTotal_min / deltaT_min);
    volInici = NaN(numD, numP);
    volFin = NaN(numD, numP);
    tEntrega = nan(1, numP);
    bateriaActual = CapBat0;
    tiempoVolConsumido = zeros(1, numD);
    figure('Name','Cronograma de Entregas','Color',[1 1 1]);
    hold on;
    baseColors = lines(numD);
    for d = 1:numD
        paqsDron = find(Eq(d, :) == 1);
        tiempoAcumSlots = 0;
        colorIndex = 0;
        for p = paqsDron
            if tiempoVolConsumido(d) >= K_horas
                break;
            end
            distTotal = 2 * dist(p);
            tVuelo_h = distTotal / max(vyD(d), 0.1);
            consumoEnergia = baseConsumD(d)*distTotal;
            % Recarga previa si falta batería
            if bateriaActual(d) < consumoEnergia
                energiaFalta = CapBat0(d) - bateriaActual(d);
                fracFalta = energiaFalta / CapBat0(d);
                tRec_h = fracFalta * TRecarregaMax;
                slotsRec = ceil((tRec_h*60)/deltaT_min);
                tiempoAcumSlots = tiempoAcumSlots + slotsRec;
                bateriaActual(d) = CapBat0(d);
            end
            if tiempoVolConsumido(d) + tVuelo_h > K_horas
                break;
            end
            slotsVuelo = ceil((tVuelo_h*60) / deltaT_min);
            startSlot = tiempoAcumSlots + 1;
            endSlot = tiempoAcumSlots + slotsVuelo;
            if endSlot > numSlots
                endSlot = numSlots;
            end
            volInici(d,p) = startSlot;
            volFin(d,p) = endSlot;
            bateriaActual(d) = bateriaActual(d) - consumoEnergia;
            tiempoVolConsumido(d) = tiempoVolConsumido(d) + tVuelo_h;
            tiempoAcumSlots = endSlot;
            % Recarga final si la batería es baja
            if bateriaActual(d) < UmbralRecarrega * CapBat0(d)
                energiaFalta = CapBat0(d) - bateriaActual(d);
                fracFalta = energiaFalta / CapBat0(d);
                tRec_h = fracFalta * TRecarregaMax;
                slotsRec = ceil((tRec_h*60)/deltaT_min);
                tiempoAcumSlots = tiempoAcumSlots + slotsRec;
                bateriaActual(d) = CapBat0(d);
            end
            t1_min = (startSlot - 1)*deltaT_min;
            t2_min = endSlot*deltaT_min;
            tEntrega_min = t1_min + 0.5*(t2_min - t1_min);
            tEntrega(p) = tEntrega_min;
            X1 = startTimeNum + t1_min/1440;
            X2 = startTimeNum + t2_min/1440;
            Xentrega = startTimeNum + tEntrega_min/1440;
            colorIndex = colorIndex + 1;
            if mod(colorIndex, 2)==0
                thisColor = baseColors(d, :);
            else
                thisColor = baseColors(d, :)*0.6;
            end
            plot([X1, X2],[d d],'LineWidth',4,'Color',thisColor);
            plot(Xentrega, d,'k*','MarkerSize',8,'LineWidth',1.5);
            xm = (X1 + X2)/2;
            ym = d;
            text(xm, ym, sprintf('P%d', p), 'HorizontalAlignment','center','VerticalAlignment','bottom');
        end
    end
    datetick('x','HH:MM','keeplimits','keepticks');
    ylabel('Dron');
    ylim([0.5, numD+0.5]);
    yticks(1:numD);
    grid on;
    title(sprintf('Cronograma de Entregas (intervalos de %d min).', deltaT_min));
    xlabel('Hora del día');
    hold off;
end

%==========================================================================
function GraficarSatisfaccio(Eq, tEntrega, horaPref, K)
    % Dibuja la curva de la función de satisfacción y superpone (en rojo)
    % los paquetes entregados.
    
    % Elegir la preferencia de referencia (mediana)
    tPref0 = median(horaPref);
    % Definir el rango de desviación en minutos
    x_min = linspace(-K*60/2, K*60/2, 300);
    % Calcular la curva de satisfacción
    Scurve = arrayfun(@(delta) CalcularSatisfaccioExacta(tPref0 + delta, tPref0, K/2), x_min);
    x_hours = x_min / 60;
    
    figure('Name','Funció de Satisfacció','Color',[1 1 1]);
    hold on;
    plot(x_hours, Scurve, 'k-', 'LineWidth', 2);
    
    % Superponer los paquetes entregados
    deliveredIdx = ~isnan(tEntrega);
    x_delivered = (tEntrega(deliveredIdx) - horaPref(deliveredIdx)) / 60;
    S_delivered = arrayfun(@(t, tp) CalcularSatisfaccioExacta(t, tp, K/2), tEntrega(deliveredIdx), horaPref(deliveredIdx));
    plot(x_delivered, S_delivered, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    xlabel('Desviació respecte a la preferència (h)');
    ylabel('Satisfacció (%)');
    title('Funció de Satisfacció i Paquets Entregats');
    grid on;
    legend('Corbaa de Satisfacció', 'Paquets Entregats','Location','Best');
    hold off;
end

%==========================================================================
function MostrarResultatsFinals(Eq, tEntregaFinal, bestIter, ...
    horaPref, numP, numD, preuD, CapBat0, distMax, ...
    costMant0, costE, costAss, pesos, dist, ...
    ncicles, alfa, costBat, tvol, mD, vyD, ...
    horesVolAny, K, costFixDron, costExternalitzacio)
    % Recalcula la satisfacción final por cada paquete
    S = zeros(1, numP);
    for p = 1:numP
        if ~isnan(tEntregaFinal(p))
            S(p) = CalcularSatisfaccioExacta(tEntregaFinal(p), horaPref(p), K);
        else
            S(p) = 0;
        end
    end

    % Calcular el coste total
    costTRec_dummy = zeros(1,numD); 
    [costosTotals, costExtTotal, costTotalGeneral] = CostFinal(...
        Eq, preuD, zeros(1,numD)+1000, ...  % vidaUtilD dummy
        CapBat0, distMax, costMant0, costE, ...
        pesos, dist, ncicles, alfa, costBat, tvol, mD, vyD, ...
        horesVolAny, costAss, K, find(all(Eq == 0,1)), ...
        costFixDron, costExternalitzacio, costTRec_dummy);

    % Calcular la satisfacción media de los paquetes entregados
    dronsAmbPaquets = any(Eq == 1, 2);
    paquetsAmbSatisfaccio = any(Eq == 1, 1);
    if isempty(paquetsAmbSatisfaccio)
        satisfFinal = 0;
    else
        S_ambPaquets = S(paquetsAmbSatisfaccio);
        satisfFinal = mean(S_ambPaquets);
    end

    % Mostrar resultados
    fprintf('\n----- RESULTADOS FINALES -----\n');
    fprintf('Mejor asignación encontrada en la iteración %d\n', bestIter);
    fprintf('Coste total general: %.2f euros\n', costTotalGeneral);
    fprintf('Coste externalización: %.2f euros\n', costExtTotal);
    fprintf('Satisfacción media global (paquetes entregados): %.2f\n', satisfFinal);
end
