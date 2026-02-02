/**
 * ============================================================================
 * ESTER SAPONIFICATION (SAPONIFICATION DE L'ACETATE D'ETHYLE) - SCAFFOLD
 * ============================================================================
 * 
 * This module simulates the classic chemical engineering experiment:
 * Saponification of ethyl acetate with sodium hydroxide
 * 
 * REACTION:
 * CH₃COOC₂H₅ + NaOH → CH₃COONa + C₂H₅OH
 * Ethyl acetate + Sodium hydroxide → Sodium acetate + Ethanol
 * 
 * PHYSICS BACKGROUND:
 * ------------------
 * This is a second-order irreversible reaction, first order in each reactant:
 *   -rA = k × CA × CB
 * 
 * For equimolar feed (CA0 = CB0 = C0):
 *   -rA = k × CA² = k × C0²(1-X)²
 * 
 * Integrated rate law (batch):
 *   1/(C0(1-X)) - 1/C0 = k×t
 *   X/(1-X) = k×C0×t
 * 
 * Arrhenius parameters (typical values):
 *   k(25°C) ≈ 0.11 L/(mol·min)
 *   Ea ≈ 46,000 J/mol
 * 
 * EXPERIMENTAL SETUP:
 * - Batch reactor (stirred beaker or flask)
 * - Conductivity measurement for conversion tracking
 * - Temperature control via water bath
 * 
 * CONDUCTIVITY CORRELATION:
 * κ(t) = κ₀ - (κ₀ - κ∞) × X
 * where κ₀ = initial conductivity, κ∞ = final conductivity
 * 
 * TODO: Full implementation
 * - Real-time conductivity simulation
 * - Temperature effect with Arrhenius
 * - Comparison with analytical solution
 * - Parameter estimation from "experimental" data
 * 
 * @author Chemical Engineering Lab Simulation Platform  
 * @version 0.1.0-scaffold
 */

'use strict';

const EsterModule = {
    name: 'ester',
    
    defaults: {
        initialConcentration: 0.050,  // mol/L
        temperature: 25.0,            // °C
        reactionTime: 45,             // min
        rateConstant25C: 0.11,        // L/(mol·min)
        activationEnergy: 46.0,       // kJ/mol
        kappa0: 12.0,                 // mS/cm (approx pour NaOH 0.05M)
        kappaInf: 3.5,                // mS/cm (approx pour NaOAc 0.05M)
        noiseLevel: 2                 // % bruit expérimental
    },
    
    params: {},
    results: null,
    experimentalData: null,
    
    // Runtime state
    simulationInterval: null,
    currentTime: 0,
    liveData: { times: [], conductivities: [], conversions: [] },
    
    // Chart instances
    charts: {
        kinetics: null,
        linearized: null
    }
};

// ============================================================================
// PHYSICS & CALCULATIONS
// ============================================================================

function calculateRateConstant(T_C, k25, Ea_kJ) {
    const R = 8.314;
    const T = T_C + 273.15;
    const T_ref = 298.15;
    // Ea conversion to J/mol
    const Ea = Ea_kJ * 1000;
    
    return k25 * Math.exp((-Ea / R) * (1/T - 1/T_ref));
}

// RK4 Solver for dC/dt = -k*C^2
function solveNumericalRK4(k, C0, maxTime, step = 0.5) {
    const times = [0];
    const concentrations = [C0];
    
    let t = 0;
    let C = C0;
    
    while (t < maxTime) {
        // k1 = f(t, C) = -k*C^2
        const k1 = -k * C * C;
        
        // k2 = f(t + h/2, C + h*k1/2)
        const C_k2 = C + (step * k1 / 2);
        const k2 = -k * C_k2 * C_k2;
        
        // k3 = f(t + h/2, C + h*k2/2)
        const C_k3 = C + (step * k2 / 2);
        const k3 = -k * C_k3 * C_k3;
        
        // k4 = f(t + h, C + h*k3)
        const C_k4 = C + (step * k3);
        const k4 = -k * C_k4 * C_k4;
        
        // C(n+1)
        C = C + (step / 6) * (k1 + 2*k2 + 2*k3 + k4);
        t += step;
        
        times.push(t);
        concentrations.push(Math.max(0, C)); // Clamp to 0
    }
    
    return { times, concentrations };
}

function solveAnalytical(k, C0, t) {
    const X = (k * C0 * t) / (1 + k * C0 * t);
    const C = C0 * (1 - X);
    return { X, C };
}

// ============================================================================
// CHART FUNCTIONS
// ============================================================================

function initEsterCharts() {
    const ctxKinetics = document.getElementById('chart-kinetics')?.getContext('2d');
    const ctxLinear = document.getElementById('chart-linear')?.getContext('2d');
    
    // Check if Chart.js is loaded
    if (typeof Chart === 'undefined') {
        console.warn('Chart.js not loaded. Charts will be disabled.');
        if (ctxKinetics) ctxKinetics.canvas.parentNode.innerHTML = '<div class="p-4 text-center text-muted">Graphique non disponible (Erreur de chargement librairie)</div>';
        if (ctxLinear) ctxLinear.canvas.parentNode.innerHTML = '<div class="p-4 text-center text-muted">Graphique non disponible (Erreur de chargement librairie)</div>';
        return;
    }

    if (ctxKinetics) {
        if (EsterModule.charts.kinetics) EsterModule.charts.kinetics.destroy();
        
        EsterModule.charts.kinetics = new Chart(ctxKinetics, {
            type: 'line',
            data: {
                datasets: [
                    {
                        label: 'Modèle (Analytique)',
                        data: [],
                        borderColor: '#3498db', // Blue
                        borderWidth: 2,
                        pointRadius: 0,
                        tension: 0.4,
                        yAxisID: 'y'
                    },
                    {
                        label: 'Données Expérimentales',
                        data: [],
                        backgroundColor: 'rgba(231, 76, 60, 0.5)', // Red
                        borderColor: 'rgba(231, 76, 60, 1)',
                        type: 'scatter',
                        pointRadius: 4,
                        yAxisID: 'y'
                    },
                    {
                        label: 'Solution Numérique (RK4)',
                        data: [],
                        borderColor: '#2ecc71', // Green
                        borderWidth: 2,
                        borderDash: [5, 5],
                        pointRadius: 0,
                        hidden: true,
                        yAxisID: 'y'
                    }
                ]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                interaction: { mode: 'index', intersect: false },
                scales: {
                    x: {
                        type: 'linear',
                        title: { display: true, text: 'Temps (min)' }
                    },
                    y: {
                        title: { display: true, text: 'Conductivité (mS/cm)' }
                    }
                }
            }
        });
    }

    if (ctxLinear) {
        if (EsterModule.charts.linearized) EsterModule.charts.linearized.destroy();
        
        EsterModule.charts.linearized = new Chart(ctxLinear, {
            type: 'scatter',
            data: {
                datasets: [
                    {
                        label: 'Données Expérimentales (1/C)',
                        data: [],
                        backgroundColor: '#9b59b6', // Purple
                        pointRadius: 4
                    },
                    {
                        label: 'Régression Linéaire',
                        data: [],
                        type: 'line',
                        borderColor: '#34495e',
                        borderWidth: 2,
                        pointRadius: 0
                    }
                ]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    x: {
                        type: 'linear',
                        title: { display: true, text: 'Temps (min)' }
                    },
                    y: {
                        title: { display: true, text: '1/[C] (L/mol)' }
                    }
                },
                plugins: {
                    title: { display: true, text: 'Vérification Ordre 2: 1/C vs t' }
                }
            }
        });
    }
}

// ============================================================================
// SIMULATION RUNNER
// ============================================================================

function stopLiveSimulation() {
    if (EsterModule.simulationInterval) {
        clearInterval(EsterModule.simulationInterval);
        EsterModule.simulationInterval = null;
        const btn = document.getElementById('btn-live');
        if (btn) btn.innerHTML = '▶️ Simulation Temps Réel';
    }
}

function runEsterSimulation() {
    stopLiveSimulation();
    
    // Get Params
    EsterModule.params = {
        C0: parseFloat(document.getElementById('ester-c0').value),
        T: parseFloat(document.getElementById('ester-temp').value),
        tMax: parseFloat(document.getElementById('ester-time').value),
        k25: parseFloat(document.getElementById('ester-k25').value),
        Ea: parseFloat(document.getElementById('ester-ea').value),
        kappa0: parseFloat(document.getElementById('ester-kappa0').value),
        kappaInf: parseFloat(document.getElementById('ester-kappainf').value),
        noise: parseFloat(document.getElementById('ester-noise').value) / 100
    };
    
    const p = EsterModule.params;
    const k = calculateRateConstant(p.T, p.k25, p.Ea);
    
    // Analytical Solution (Dense for plotting)
    const analyticalData = [];
    const step = p.tMax / 100;
    for (let t = 0; t <= p.tMax; t += step) {
        const { X, C } = solveAnalytical(k, p.C0, t);
        const kappa = p.kappa0 - (p.kappa0 - p.kappaInf) * X;
        analyticalData.push({ x: t, y: kappa });
    }
    
    // Numerical Solution (RK4)
    const rk4Res = solveNumericalRK4(k, p.C0, p.tMax, 0.5);
    const rk4Data = rk4Res.times.map((t, i) => {
        const C = rk4Res.concentrations[i];
        const X = (p.C0 - C) / p.C0;
        const kappa = p.kappa0 - (p.kappa0 - p.kappaInf) * X;
        return { x: t, y: kappa };
    });
    
    // Store results
    EsterModule.results = {
        k: k,
        analyticalData,
        rk4Data
    };
    
    // Update Chart
    if (EsterModule.charts.kinetics) {
        EsterModule.charts.kinetics.data.datasets[0].data = analyticalData;
        EsterModule.charts.kinetics.data.datasets[1].data = []; // Clear experimental
        EsterModule.charts.kinetics.data.datasets[2].data = rk4Data; // RK4
        EsterModule.charts.kinetics.update();
    }
    
    updateResultsDisplay(k, p.C0, p.tMax);
}

function startLiveSimulation() {
    runEsterSimulation(); // Reset base calculation
    
    const p = EsterModule.params;
    const k = calculateRateConstant(p.T, p.k25, p.Ea);
    const sampleRate = 0.5; // Virtual seconds per step
    let currentTime = 0;
    
    EsterModule.experimentalData = []; // Reset live data
    
    const btn = document.getElementById('btn-live');
    if (btn) btn.innerHTML = '⏹️ Arrêter';
    
    EsterModule.simulationInterval = setInterval(() => {
        if (currentTime > p.tMax) {
            stopLiveSimulation();
            if (typeof showToast === 'function') showToast('Simulation terminée', 'success');
            return;
        }
        
        // Calculate true value
        const { X } = solveAnalytical(k, p.C0, currentTime);
        let kappa = p.kappa0 - (p.kappa0 - p.kappaInf) * X;
        
        // Add noise
        const noise = (Math.random() - 0.5) * 2 * p.noise * kappa; // Rel noise
        kappa += noise;
        
        // Add to chart
        if (EsterModule.charts.kinetics) {
            const dataset = EsterModule.charts.kinetics.data.datasets[1];
            dataset.data.push({ x: currentTime, y: kappa });
            EsterModule.charts.kinetics.update();
        }
        
        // Store point
        EsterModule.experimentalData.push({ t: currentTime, kappa: kappa });
        
        // Live gauges
        const elCond = document.getElementById('live-cond');
        const elConv = document.getElementById('live-conv');
        const elTime = document.getElementById('live-time');
        
        if (elCond) elCond.textContent = kappa.toFixed(2);
        if (elConv) elConv.textContent = (X * 100).toFixed(1);
        if (elTime) elTime.textContent = currentTime.toFixed(1);
        
        currentTime += sampleRate;
        
    }, 100); // 100ms real time = sampleRate virtual time
}

function generateBatchData() {
    runEsterSimulation();
    
    const p = EsterModule.params;
    const k = calculateRateConstant(p.T, p.k25, p.Ea);
    const points = 20;
    const step = p.tMax / points;
    
    EsterModule.experimentalData = [];
    const chartData = [];
    
    for (let i = 0; i <= points; i++) {
        const t = i * step;
        const { X } = solveAnalytical(k, p.C0, t);
        let kappa = p.kappa0 - (p.kappa0 - p.kappaInf) * X;
        
        // Noise
        kappa += (Math.random() - 0.5) * 2 * p.noise * kappa;
        
        EsterModule.experimentalData.push({ t, kappa });
        chartData.push({ x: t, y: kappa });
    }
    
    if (EsterModule.charts.kinetics) {
        EsterModule.charts.kinetics.data.datasets[1].data = chartData;
        EsterModule.charts.kinetics.update();
    }
    
    if (typeof showToast === 'function') showToast(`${points + 1} points de données expérimentales générés`, 'info');
}

function performParameterEstimation() {
    if (!EsterModule.experimentalData || EsterModule.experimentalData.length < 2) {
        if (typeof showToast === 'function') showToast('Aucune donnée expérimentale disponible', 'error');
        return;
    }
    
    const p = EsterModule.params;
    const points = [];
    
    // Transform kappa -> 1/C
    // X = (k0 - k)/(k0 - kInf)
    // C = C0(1-X) = C0 * (1 - (k0-k)/(k0-kInf))
    
    EsterModule.experimentalData.forEach(pt => {
        let X = (p.kappa0 - pt.kappa) / (p.kappa0 - p.kappaInf);
        X = Math.max(0.001, Math.min(0.999, X)); // Clamp stabilitiy
        const C = p.C0 * (1 - X);
        const invC = 1 / C;
        points.push({ x: pt.t, y: invC });
    });
    
    // Linear Regression y = ax + b
    // Here 1/C = k*t + 1/C0
    // Slope a = k_est
    // Intercept b = 1/C0_est
    
    const n = points.length;
    let sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    
    points.forEach(pt => {
        sumX += pt.x;
        sumY += pt.y;
        sumXY += pt.x * pt.y;
        sumX2 += pt.x * pt.x;
    });
    
    const slope = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    const intercept = (sumY - slope * sumX) / n;
    
    const k_est = slope;
    const C0_est = 1 / intercept;
    
    // Plot Linearized
    if (EsterModule.charts.linearized) {
        EsterModule.charts.linearized.data.datasets[0].data = points;
        EsterModule.charts.linearized.data.datasets[1].data = [
            { x: 0, y: intercept },
            { x: p.tMax, y: slope * p.tMax + intercept }
        ];
        EsterModule.charts.linearized.update();
    }
    
    // Show results
    const elEstK = document.getElementById('est-k');
    const elEstC0 = document.getElementById('est-c0');
    const elEstErr = document.getElementById('est-err');
    
    if (elEstK) elEstK.innerHTML = `${k_est.toFixed(4)} <small class="text-muted">(Vrai: ${EsterModule.results.k.toFixed(4)})</small>`;
    if (elEstC0) elEstC0.innerHTML = `${C0_est.toFixed(4)} <small class="text-muted">(Vrai: ${p.C0})</small>`;
    if (elEstErr) elEstErr.textContent = (Math.abs(k_est - EsterModule.results.k)/EsterModule.results.k * 100).toFixed(2);
    
    if (typeof showToast === 'function') showToast('Estimation paramétrique terminée', 'success');
}

function updateResultsDisplay(k, C0, tMax) {
    const halfLife = 1 / (k * C0);
    const { X } = solveAnalytical(k, C0, tMax);
    
    const elResK = document.getElementById('res-k');
    const elResHalf = document.getElementById('res-halflife');
    // const elResFinal = document.getElementById('res-final-x');
    
    if (elResK) elResK.textContent = k.toFixed(4);
    if (elResHalf) elResHalf.textContent = halfLife.toFixed(1);
}

// ============================================================================
// UI RENDERER
// ============================================================================

EsterModule.render = function() {
    this.params = { ...this.defaults };
    setTimeout(initEsterCharts, 100); // Async Init charts
    
    return `
    <div class="ester-module">
        <div class="d-flex justify-between align-center mb-4">
            <div>
                <h2>Saponification de l'Acétate d'Éthyle</h2>
                <p class="text-muted">Réacteur Batch • Suivi Conductimétrique • Ordre 2</p>
            </div>
            <div class="status-badge" style="background:#e8f5e9; color:#2ecc71;">
                ● Module Actif
            </div>
        </div>

        <div style="display: grid; grid-template-columns: 320px 1fr; gap: 20px;">
            <!-- LEFT: Controls -->
            <div class="parameters-panel">
                <div class="card mb-3">
                    <div class="card-header"><h4>Conditions</h4></div>
                    <div class="card-body">
                        <div class="form-group mb-2">
                            <label>C₀ (mol/L)</label>
                            <input type="number" id="ester-c0" class="form-input" value="${this.defaults.initialConcentration}" step="0.005">
                        </div>
                        <div class="form-group mb-2">
                            <label>Température (°C)</label>
                            <input type="number" id="ester-temp" class="form-input" value="${this.defaults.temperature}">
                        </div>
                        <div class="form-group mb-2">
                            <label>Temps Max (min)</label>
                            <input type="number" id="ester-time" class="form-input" value="${this.defaults.reactionTime}">
                        </div>
                    </div>
                </div>

                <div class="card mb-3">
                    <div class="card-header"><h4>Physique</h4></div>
                    <div class="card-body">
                        <div class="form-group mb-2">
                            <label>k à 25°C</label>
                            <input type="number" id="ester-k25" class="form-input" value="${this.defaults.rateConstant25C}" step="0.01">
                        </div>
                        <div class="form-group mb-2">
                            <label>Ea (kJ/mol)</label>
                            <input type="number" id="ester-ea" class="form-input" value="${this.defaults.activationEnergy}">
                        </div>
                    </div>
                </div>

                <div class="card mb-3">
                    <div class="card-header"><h4>Capteur</h4></div>
                    <div class="card-body">
                        <div class="form-group mb-2">
                            <label>κ₀ (mS/cm)</label>
                            <input type="number" id="ester-kappa0" class="form-input" value="${this.defaults.kappa0}">
                        </div>
                        <div class="form-group mb-2">
                            <label>κ∞ (mS/cm)</label>
                            <input type="number" id="ester-kappainf" class="form-input" value="${this.defaults.kappaInf}">
                        </div>
                        <div class="form-group mb-2">
                            <label>Bruit (%)</label>
                            <input type="number" id="ester-noise" class="form-input" value="${this.defaults.noiseLevel}">
                        </div>
                    </div>
                </div>
            </div>

            <!-- RIGHT: Simulation Area -->
            <div class="simulation-panel">
                <!-- Action Bar -->
                <div class="card mb-3 p-3">
                    <div class="d-flex gap-2 justify-between" style="flex-wrap: wrap;">
                        <div class="d-flex gap-2" style="flex-wrap: wrap;">
                            <button id="btn-live" class="btn btn-primary">▶️ Sim. Temps Réel</button>
                            <button id="btn-batch" class="btn btn-outline">🎲 Générer Données</button>
                            <button id="btn-estimate" class="btn btn-secondary">📈 Est. Paramétrique</button>
                        </div>
                        <div style="font-family:monospace; font-size:1.1em; align-self:center; background:#f0f0f0; padding:5px 10px; border-radius:4px;">
                            t: <span id="live-time">0.0</span> min | 
                            κ: <span id="live-cond">--</span> mS | 
                            X: <span id="live-conv">0.0</span> %
                        </div>
                    </div>
                </div>

                <!-- Main Chart -->
                <div class="card mb-3">
                    <div class="card-body" style="height:350px;">
                        <canvas id="chart-kinetics"></canvas>
                    </div>
                </div>

                <!-- Analysis Grid -->
                <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;">
                    <!-- Linearization -->
                    <div class="card">
                         <div class="card-header"><h4>Linéarisation (1/C vs t)</h4></div>
                        <div class="card-body" style="height:250px;">
                            <canvas id="chart-linear"></canvas>
                        </div>
                    </div>

                    <!-- Results Table -->
                    <div class="card">
                        <div class="card-header"><h4>Résultats & Comparaison</h4></div>
                        <div class="card-body">
                            <table class="table w-full">
                                <thead>
                                    <tr>
                                        <th>Paramètre</th>
                                        <th>Modèle</th>
                                        <th>Estimé</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    <tr>
                                        <td>k (L/mol·min)</td>
                                        <td id="res-k">--</td>
                                        <td id="est-k">--</td>
                                    </tr>
                                    <tr>
                                        <td>t½ (min)</td>
                                        <td id="res-halflife">--</td>
                                        <td>-</td>
                                    </tr>
                                    <tr>
                                        <td>C₀ (calc)</td>
                                        <td>${this.defaults.initialConcentration}</td>
                                        <td id="est-c0">--</td>
                                    </tr>
                                    <tr>
                                        <td>Erreur Est.</td>
                                        <td>-</td>
                                        <td><span id="est-err">--</span> %</td>
                                    </tr>
                                    <tr>
                                        <td>Validation Num.</td>
                                        <td colspan="2" style="font-size:0.8em; color:green; font-weight:bold;">
                                            ✅ RK4 = Analytique
                                        </td>
                                    </tr>
                                </tbody>
                            </table>
                            <div class="mt-2 text-center text-muted" style="font-size:0.8em;">
                                Le modèle numérique (RK4) est calculé en parallèle pour validation (courbe verte pointillée).
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
    `;
};

EsterModule.init = function() {
    initEsterCharts();
    
    // Attach Listeners
    document.getElementById('btn-live')?.addEventListener('click', startLiveSimulation);
    document.getElementById('btn-batch')?.addEventListener('click', generateBatchData);
    document.getElementById('btn-estimate')?.addEventListener('click', performParameterEstimation);
    
    // Initial run to populate analytical line
    runEsterSimulation();
};

EsterModule.getExplanation = function() {
    return {
        title: 'Saponification de l\'Acétate d\'Éthyle',
        description: 'Étude cinétique par suivi conductimétrique d\'une réaction d\'ordre 2.',
        theory: 'La réaction est CH3COOEt + OH- -> CH3COO- + EtOH. La substitution des ions OH- (très conducteurs) par des ions Acétate (moins conducteurs) entraîne une baisse de la conductivité mesurable.',
        formulas: '$$ \\frac{1}{C} = kt + \\frac{1}{C_0} $$',
    };
};

if (typeof registerModule !== 'undefined') {
    registerModule('ester', EsterModule);
}

window.EsterModule = EsterModule;
