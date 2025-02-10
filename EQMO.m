function [costTotalGeneral, satisfFinal, bestIter] = ...
    EQMO(porcentajeOpt, porcentajeGroup, ...
                                         numIteracions, maxNoImprovement,numP,numD)
    fprintf('>> DistribucioMultiobjectiu4c_Modified2\n');
    rng(145763);
    
    % numP = 100;  
    % numD = 25;
    K    = 12;   % Horas totals de la jornada
    
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
    
    % --- MODIFICACIÓ: Ajust mínim en la lògica d'agrupament ---
    % Llavors, volem que:
    %   - Si porcentajeGroup <= 0, cada paquet es col·loca en el seu grup (grupTamany = 1).
    %   - Si porcentajeGroup >= 100, tots els paquets van en un sol grup (grupTamany = numP).
    % Per tant:
    if porcentajeGroup < 0
        grupTamany = 1;          % 0%: cada paquet en grup individual (sense agrupament)
    elseif porcentajeGroup > 100
        grupTamany = numP;       % 100%: tots els paquets en un sol grup
    else
        % Interpolar de forma que a mesura que porcentajeGroup augmenta, grupTamany
        % creixi (i per tant el nombre de grups disminueixi).
        grupTamany = max(1, round(numP * ((100 - porcentajeGroup) / 100)));
    end

    % -----------------------------------------------------------
    
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
end

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
    
    costMaxSim = numP*costExternalitzacio*sum(dist);  % LLindar maxim per externalització
    
    % Inicialització
    Eq = zeros(numD, numP);
    tEntrega = nan(1, numP);
    
    bestWeightedSum = -realmax; 
    bestEq = Eq;
    best_tEntrega = tEntrega;
    bestIter = 0;
    noImprovementCounter = 0;
    
    % Divisió en grups
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
    
    % Aleatorietat
    indicesGrupos = randperm(numGrupos);
    gruposSatisf = indicesGrupos(1:numGruposSatisf);
    gruposCoste  = indicesGrupos(numGruposSatisf+1:end);
    asignacionGrupos = zeros(1, numGrupos);
    asignacionGrupos(gruposSatisf) = 1; 
    
    % Historial per graficar
    costHistory = zeros(1, numIteracions);
    satisfHistory = zeros(1, numIteracions);
    
    for it = 1:numIteracions
        fprintf('\nITERACIÓ #%d de %d\n', it, numIteracions);
        
        if it > 1 && mod(it, ceil(numIteracions / 10)) == 0
            numReasignar = max(1, round(0.1 * numGrupos));
            gruposAReasignar = randperm(numGrupos, numReasignar);
            for g = gruposAReasignar
                if asignacionGrupos(g) == 1 && numGruposSatisf > 0
                    asignacionGrupos(g) = 0;  % canviar a Coste
                    numGruposSatisf = numGruposSatisf - 1;
                    numGruposCoste = numGruposCoste + 1;
                elseif asignacionGrupos(g) == 0 && numGruposSatisf < gruposTotal
                    asignacionGrupos(g) = 1;  % canviar a Satisf
                    numGruposSatisf = numGruposSatisf + 1;
                    numGruposCoste = numGruposCoste - 1;
                end
            end
            fprintf('  -> Reasignació de %d grups.\n', numReasignar);
        end
        
        for g = 1:numGrupos
            paqsG = grupos{g};
            
            if asignacionGrupos(g) == 1
                % Satisfacció
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
        
        % Satisfacció mitjana
        SIter = zeros(1, numP);
        for p = 1:numP
            if ~isnan(tEntrega(p))
                SIter(p) = CalcularSatisfaccioExacta(tEntrega(p), horaPref(p), K);
            end
        end
        satIter = mean(SIter);
        
        costHistory(it) = cTotal;
        satisfHistory(it) = satIter;
     
        satNorm = satIter / 100;   % satisfacció va de 0 a 100
        costNorm = min(1, cTotal / costMaxSim); 
        alpha = (porcentajeOpt / 100);  % pes de satisfacció
        
        % Maximitzar la combinació => més satisfacció i menys cost
        % (1 - costNorm) farà que a cost alt => el valor baixi

        weightedSum = alpha * satNorm + (1 - alpha) * (1 - costNorm);
        
        fprintf('Fin iter %d: Coste= %.2f, Satisf= %.2f, WeightedSum= %.4f\n', ...
            it, cTotal, satIter, weightedSum);
        
        if (weightedSum > bestWeightedSum)
            bestWeightedSum = weightedSum;
            bestEq = Eq;
            best_tEntrega = tEntrega;
            bestIter = it;
            noImprovementCounter = 0;
            fprintf('Mejor asignació trobada en la iteració %d.\n', it);
        else
            noImprovementCounter = noImprovementCounter + 1;
            fprintf('Sin mejora significativa. Contador: %d\n', noImprovementCounter);
        end
        
        if noImprovementCounter >= maxNoImprovement
            fprintf('Convergència assolida en la iteració %d.\n', it);
            break;
        end
    end
    
    % Millor solució trobada
    EqFinal = bestEq;
    tEntregaFinal = best_tEntrega;
    
    assignin('base', 'costHistory', costHistory(1:bestIter));
    assignin('base', 'satisfHistory', satisfHistory(1:bestIter));
end

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
                
                % Calcular cost incremental
                distTotal = 2 * dist(p);
                tiempoVueloNeces = distTotal / vyD(d);
                consumoEnergia = baseConsumD(d) * distTotal * (1 + PesFactor * pesos(p));
                
                if BateriaRestant(d) < consumoEnergia
                    faltaEnergia = consumoEnergia - BateriaRestant(d);
                    tiempoRecarga = (faltaEnergia / CapBat0(d)) * TRecarregaMax;
                else
                    tiempoRecarga = 0;
                end
                
                % Restriccions
                if (tvol(d) + tiempoRecarga + tiempoVueloNeces) > K
                    continue;  
                end
                
                EqSim = Eq;
                EqSim(:, p) = 0;
                EqSim(d, p) = 1;
                
                % Calcular cost
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
                
                % Actualitzar estat del dron
                if BateriaRestant(millorDron) < consumoEnergia
                    BateriaRestant(millorDron) = CapBat0(millorDron);
                    tvol(millorDron) = tvol(millorDron) + tiempoRecarga;
                end
                
                tvol(millorDron) = tvol(millorDron) + tiempoVueloNeces;
                BateriaRestant(millorDron) = BateriaRestant(millorDron) - consumoEnergia;
                
                canvi = true;
            end
        end
        
        % Calcular cost total actual
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

function [cVec, cExt, cTotal] = CostFinal(Eq, preuD, vidaUtilD, ...
    CapBat0, distMax, costMant0, costE, pesos, dist, ...
    ncicles, alfa, costBat, tvol, mD, vyD, ...
    horesVolAny, costAss, K, paqsNoAsig, ...
    costFixDron, costExternalitzacio, costTRec)
    
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

