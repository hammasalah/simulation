/**
 * ============================================================================
 * DISTILLATION COLUMN SIMULATION - SCAFFOLD
 * ============================================================================
 * 
 * This module provides distillation column design using:
 *   - McCabe-Thiele graphical method
 *   - HETP (Height Equivalent to Theoretical Plate) for packed columns
 * 
 * PHYSICS BACKGROUND:
 * ------------------
 * Distillation separates liquid mixtures based on differences in volatility.
 * 
 * KEY CONCEPTS:
 * 
 * Vapor-Liquid Equilibrium (VLE):
 *   y* = αx / (1 + (α-1)x)  [constant relative volatility]
 *   where α = volatility of light component / volatility of heavy component
 * 
 * McCabe-Thiele Assumptions:
 *   - Constant molar overflow (CMO)
 *   - Equal latent heats of vaporization
 *   - Adiabatic column
 * 
 * Operating Lines:
 *   Rectifying: y = (R/(R+1))x + xD/(R+1)
 *   Stripping: y = (L̄/V̄)x - (L̄-B)/V̄ × xB
 *   
 * where R = L/D (reflux ratio), xD = distillate composition, xB = bottoms
 * 
 * Minimum Reflux Ratio (Rmin):
 *   Found at pinch point where operating line touches equilibrium curve
 * 
 * Theoretical Stages:
 *   Step off between operating lines and equilibrium curve
 * 
 * Packed Columns:
 *   Z = HETP × N_theoretical
 *   HETP depends on packing type, vapor/liquid loading
 * 
 * TODO: Full implementation
 * - Non-ideal VLE (activity coefficients)
 * - Murphree tray efficiency
 * - Enthalpy-concentration (Ponchon-Savarit) method
 * - Pressure drop calculations
 * - Multi-component distillation (FUG method)
 * 
 * @author Chemical Engineering Lab Simulation Platform
 * @version 0.1.0-scaffold
 */

'use strict';

const DistillationModule = {
    name: 'distillation',
    
    defaults: {
        feedComposition: 0.5,        // zF (mole fraction light component)
        distillateComposition: 0.95, // xD
        bottomsComposition: 0.05,    // xB
        feedCondition: 1.0,          // q (1 = sat. liquid, 0 = sat. vapor)
        relativeVolatility: 2.5,     // α
        refluxRatio: 2.0,            // R = L/D
        columnType: 'tray',          // 'tray' or 'packed'
        hetp: 0.5,                   // m (for packed columns)
        feedFlowRate: 100,           // kmol/h
        trayEfficiency: 0.7          // Murphree efficiency
    },
    
    params: null,
    results: null,
    charts: {}
};

// ============================================================================
// VLE CALCULATIONS
// ============================================================================

/**
 * Calculate equilibrium vapor composition for constant relative volatility
 * 
 * @param {number} x - Liquid mole fraction
 * @param {number} alpha - Relative volatility
 * @returns {number} Equilibrium vapor mole fraction y*
 */
function equilibriumY(x, alpha) {
    return (alpha * x) / (1 + (alpha - 1) * x);
}

/**
 * Generate equilibrium curve data points
 * 
 * @param {number} alpha - Relative volatility
 * @returns {Array} Array of {x, y} points
 */
function generateEquilibriumCurve(alpha) {
    const points = [];
    for (let i = 0; i <= 100; i++) {
        const x = i / 100;
        const y = equilibriumY(x, alpha);
        points.push({ x, y });
    }
    return points;
}

// ============================================================================
// MCCABE-THIELE METHOD
// ============================================================================

/**
 * Calculate minimum reflux ratio at pinch point
 * 
 * For feed at bubble point (q=1):
 * Rmin = (xD - y'eq) / (y'eq - x'eq)
 * where (x', y') is intersection of q-line with equilibrium curve
 * 
 * @param {Object} params - Column parameters
 * @returns {number} Minimum reflux ratio
 */
function calculateMinReflux(params) {
    const { feedComposition, distillateComposition, feedCondition, relativeVolatility } = params;
    
    const zF = feedComposition;
    const xD = distillateComposition;
    const q = feedCondition;
    const alpha = relativeVolatility;
    
    // q-line: y = (q/(q-1))x - zF/(q-1)
    // For q=1 (saturated liquid), q-line is vertical at x = zF
    
    // At x = zF (for q=1)
    const yEq = equilibriumY(zF, alpha);
    
    // Minimum reflux: operating line passes through (zF, yEq) and (xD, xD)
    const Rmin = (xD - yEq) / (yEq - zF);
    
    return Math.max(0, Rmin);
}

/**
 * Calculate number of theoretical stages using McCabe-Thiele
 * 
 * @param {Object} params - Column parameters
 * @returns {Object} Results including number of stages
 */
function mcCabeThiele(params) {
    const {
        feedComposition,
        distillateComposition,
        bottomsComposition,
        feedCondition,
        relativeVolatility,
        refluxRatio
    } = params;
    
    const zF = feedComposition;
    const xD = distillateComposition;
    const xB = bottomsComposition;
    const q = feedCondition;
    const alpha = relativeVolatility;
    const R = refluxRatio;
    
    // Calculate minimum reflux
    const Rmin = calculateMinReflux(params);
    
    // Check if specified R > Rmin
    if (R < Rmin) {
        return {
            error: true,
            message: `Reflux ratio (${R.toFixed(2)}) must be greater than Rmin (${Rmin.toFixed(2)})`
        };
    }
    
    // Rectifying section operating line: y = (R/(R+1))x + xD/(R+1)
    const rectSlope = R / (R + 1);
    const rectIntercept = xD / (R + 1);
    
    // q-line intersection with rectifying line
    // For saturated liquid (q=1): x = zF
    const xIntersect = zF;
    const yIntersect = rectSlope * xIntersect + rectIntercept;
    
    // Stripping section operating line passes through (xB, xB) and (xIntersect, yIntersect)
    const stripSlope = (yIntersect - xB) / (xIntersect - xB);
    
    // Step off stages from top (xD) to bottom (xB)
    const stages = [];
    let x = xD;
    let y = xD;
    let n = 0;
    let feedStage = 0;
    
    while (x > xB && n < 100) {  // Max 100 stages safety
        n++;
        
        // From y, find equilibrium x (step horizontally to eq curve)
        // Solve: y = αx/(1+(α-1)x) for x
        // x = y / (α - y(α-1))
        const xEq = y / (alpha - y * (alpha - 1));
        x = Math.max(xB, xEq);
        
        stages.push({ n, x: xEq, y });
        
        // Determine which operating line to use
        if (x > xIntersect) {
            // Rectifying section
            y = rectSlope * x + rectIntercept;
        } else {
            // Stripping section
            if (feedStage === 0) feedStage = n;
            y = stripSlope * (x - xB) + xB;
        }
        
        y = Math.max(xB, Math.min(1, y));
    }
    
    return {
        theoreticalStages: n,
        feedStage,
        Rmin,
        RoverRmin: R / Rmin,
        stages,
        operatingLines: {
            rectifying: { slope: rectSlope, intercept: rectIntercept },
            stripping: { slope: stripSlope, intercept: xB },
            qLine: { x: xIntersect, y: yIntersect }
        }
    };
}

/**
 * Calculate column height for packed column
 * 
 * @param {Object} params - Column parameters
 * @returns {Object} Height calculation results
 */
function calculatePackedHeight(params) {
    const mcResult = mcCabeThiele(params);
    
    if (mcResult.error) return mcResult;
    
    const N = mcResult.theoreticalStages;
    const HETP = params.hetp;
    
    const packingHeight = N * HETP;
    
    return {
        ...mcResult,
        hetp: HETP,
        packingHeight,
        packingType: 'Structured (assumed)'
    };
}

// ============================================================================
// MAIN CALCULATION
// ============================================================================

function calculateDistillation(params) {
    // Material balance
    const F = params.feedFlowRate;
    const zF = params.feedComposition;
    const xD = params.distillateComposition;
    const xB = params.bottomsComposition;
    
    // D and B from material balance
    // F = D + B
    // F*zF = D*xD + B*xB
    const D = F * (zF - xB) / (xD - xB);
    const B = F - D;
    
    // McCabe-Thiele
    const mcResult = params.columnType === 'packed' 
        ? calculatePackedHeight(params) 
        : mcCabeThiele(params);
    
    if (mcResult.error) {
        return mcResult;
    }
    
    // Actual stages with tray efficiency
    const actualStages = params.columnType === 'tray'
        ? Math.ceil(mcResult.theoreticalStages / params.trayEfficiency)
        : mcResult.theoreticalStages;
    
    const calculationSteps = [
        {
            description: 'Bilan matière global',
            formula: `F = D + B → ${F.toFixed(1)} = D + B`,
            result: `D = ${D.toFixed(2)} kmol/h, B = ${B.toFixed(2)} kmol/h`
        },
        {
            description: 'Reflux minimum',
            formula: `Rmin calculé à partir du pincement`,
            result: `Rmin = ${mcResult.Rmin.toFixed(3)}`
        },
        {
            description: 'Rapport R/Rmin',
            formula: `R/Rmin = ${params.refluxRatio.toFixed(2)} / ${mcResult.Rmin.toFixed(3)}`,
            result: `R/Rmin = ${mcResult.RoverRmin.toFixed(2)}`
        },
        {
            description: 'Nombre de plateaux théoriques',
            formula: 'Méthode de McCabe-Thiele',
            result: `N = ${mcResult.theoreticalStages} plateaux`
        },
        {
            description: 'Étage d\'alimentation optimal',
            formula: 'Intersection des droites opératoires',
            result: `Étage n°${mcResult.feedStage} (depuis le haut)`
        }
    ];
    
    if (params.columnType === 'packed') {
        calculationSteps.push({
            description: 'Hauteur de garnissage',
            formula: `Z = N × HETP = ${mcResult.theoreticalStages} × ${params.hetp}`,
            result: `Z = ${mcResult.packingHeight.toFixed(2)} m`
        });
    } else {
        calculationSteps.push({
            description: 'Nombre de plateaux réels',
            formula: `N_réel = N_théo / η = ${mcResult.theoreticalStages} / ${params.trayEfficiency}`,
            result: `N_réel = ${actualStages} plateaux`
        });
    }
    
    return {
        inputs: { ...params },
        materialBalance: { F, D, B, zF, xD, xB },
        ...mcResult,
        actualStages,
        calculationSteps,
        equilibriumCurve: generateEquilibriumCurve(params.relativeVolatility),
        timestamp: new Date().toISOString()
    };
}

// ============================================================================
// UI RENDERING - SCAFFOLD
// ============================================================================

DistillationModule.render = function() {
    this.params = { ...this.defaults };
    
    return `
        <div class="distillation-module" id="distillation-content">
            <!-- Header -->
            <div class="d-flex justify-between align-center mb-5">
                <div>
                    <h2>Simulation de Distillation</h2>
                    <p class="text-muted">
                        Dimensionnement de colonnes à plateaux et garnies par 
                        la méthode de McCabe-Thiele.
                    </p>
                </div>
                <div class="status-badge idle" id="simulation-status">
                    <span class="status-dot"></span>
                    <span>Prêt</span>
                </div>
            </div>
            
            <!-- Development Notice -->
            <div class="card mb-5" style="border-color: var(--warning); background: rgba(243, 156, 18, 0.1);">
                <div class="card-body">
                    <h3>⚠️ Module en cours de développement</h3>
                    <p>Fonctionnalités à implémenter :</p>
                    <ul>
                        <li>VLE non-idéal (coefficients d'activité)</li>
                        <li>Méthode de Ponchon-Savarit</li>
                        <li>Distillation multi-composants (FUG)</li>
                        <li>Calcul des charges de perte et dimensionnement hydraulique</li>
                    </ul>
                </div>
            </div>
            
            <!-- Main Layout -->
            <div style="display: grid; grid-template-columns: 1fr 1fr; gap: var(--space-5);">
                <!-- Left: Inputs -->
                <div>
                    <div class="card mb-5">
                        <div class="card-header">
                            <h3 class="card-title">Type de colonne</h3>
                        </div>
                        <div class="card-body">
                            <div class="mode-selector" style="width: 100%;">
                                <button class="mode-btn active" data-column="tray">Plateaux</button>
                                <button class="mode-btn" data-column="packed">Garnie</button>
                            </div>
                            
                            <div class="diagram-container mt-4" id="distillation-diagram">
                                <!-- SVG will load here -->
                            </div>
                        </div>
                    </div>
                    
                    <div class="card">
                        <div class="card-header">
                            <h3 class="card-title">Paramètres</h3>
                        </div>
                        <div class="card-body">
                            <div class="param-grid">
                                <div class="form-group">
                                    <label class="form-label">Composition alimentation (zF)</label>
                                    <div class="form-range-container">
                                        <input type="range" class="form-range" id="feed-comp"
                                               value="50" min="10" max="90" step="1">
                                        <span class="range-value">0.50</span>
                                    </div>
                                </div>
                                
                                <div class="form-group">
                                    <label class="form-label">Composition distillat (xD)</label>
                                    <div class="form-range-container">
                                        <input type="range" class="form-range" id="distillate-comp"
                                               value="95" min="80" max="99" step="1">
                                        <span class="range-value">0.95</span>
                                    </div>
                                </div>
                                
                                <div class="form-group">
                                    <label class="form-label">Composition résidu (xB)</label>
                                    <div class="form-range-container">
                                        <input type="range" class="form-range" id="bottoms-comp"
                                               value="5" min="1" max="20" step="1">
                                        <span class="range-value">0.05</span>
                                    </div>
                                </div>
                                
                                <div class="form-group">
                                    <label class="form-label">Volatilité relative (α)</label>
                                    <div class="input-with-unit">
                                        <input type="number" class="form-input" id="alpha"
                                               value="2.5" min="1.1" max="10" step="0.1">
                                        <span class="input-unit">-</span>
                                    </div>
                                </div>
                                
                                <div class="form-group">
                                    <label class="form-label">Taux de reflux (R)</label>
                                    <div class="input-with-unit">
                                        <input type="number" class="form-input" id="reflux"
                                               value="2.0" min="0.5" max="10" step="0.1">
                                        <span class="input-unit">-</span>
                                    </div>
                                </div>
                                
                                <div class="form-group">
                                    <label class="form-label">Débit alimentation</label>
                                    <div class="input-with-unit">
                                        <input type="number" class="form-input" id="feed-flow"
                                               value="100" min="10" max="1000" step="10">
                                        <span class="input-unit">kmol/h</span>
                                    </div>
                                </div>
                            </div>
                            
                            <div class="btn-group mt-5">
                                <button class="btn btn-primary" id="run-simulation">
                                    ▶️ Calculer
                                </button>
                                <button class="btn btn-outline" id="load-example">
                                    💡 Benzène/Toluène
                                </button>
                            </div>
                        </div>
                    </div>
                </div>
                
                <!-- Right: Results -->
                <div>
                    <div class="card mb-5">
                        <div class="card-header">
                            <h3 class="card-title">Résultats</h3>
                        </div>
                        <div class="card-body">
                            <div class="results-grid">
                                <div class="result-item">
                                    <div class="result-label">Plateaux théoriques</div>
                                    <div class="result-value" id="result-stages">--</div>
                                    <div class="result-unit">-</div>
                                </div>
                                <div class="result-item">
                                    <div class="result-label">Étage alimentation</div>
                                    <div class="result-value" id="result-feed-stage">--</div>
                                    <div class="result-unit">-</div>
                                </div>
                                <div class="result-item">
                                    <div class="result-label">Rmin</div>
                                    <div class="result-value" id="result-rmin">--</div>
                                    <div class="result-unit">-</div>
                                </div>
                                <div class="result-item">
                                    <div class="result-label">R/Rmin</div>
                                    <div class="result-value" id="result-r-ratio">--</div>
                                    <div class="result-unit">-</div>
                                </div>
                                <div class="result-item">
                                    <div class="result-label">Débit distillat</div>
                                    <div class="result-value" id="result-d">--</div>
                                    <div class="result-unit">kmol/h</div>
                                </div>
                                <div class="result-item">
                                    <div class="result-label">Débit résidu</div>
                                    <div class="result-value" id="result-b">--</div>
                                    <div class="result-unit">kmol/h</div>
                                </div>
                            </div>
                            
                            <div class="mt-4">
                                <button class="btn btn-outline w-full" id="toggle-calculations">
                                    📝 Afficher les calculs
                                </button>
                            </div>
                            <div class="calc-steps mt-4 hidden" id="calculation-steps"></div>
                        </div>
                    </div>
                    
                    <div class="card">
                        <div class="card-header">
                            <h3 class="card-title">Diagramme McCabe-Thiele</h3>
                        </div>
                        <div class="card-body">
                            <div class="chart-container tall" id="mccabe-thiele-chart">
                                <p class="text-muted text-center">
                                    TODO: Diagramme x-y avec courbe d'équilibre et droites opératoires
                                </p>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    `;
};

DistillationModule.init = async function() {
    // Load diagram
    const container = document.getElementById('distillation-diagram');
    try {
        const response = await fetch('assets/distillation.svg');
        if (response.ok) {
            container.innerHTML = await response.text();
        }
    } catch (e) {
        container.innerHTML = '<p class="text-muted">Schéma non disponible</p>';
    }
    
    setupDistillationEventListeners();
    runDistillationSimulation();
};

function setupDistillationEventListeners() {
    document.getElementById('run-simulation')?.addEventListener('click', runDistillationSimulation);
    
    document.querySelectorAll('[data-column]').forEach(btn => {
        btn.addEventListener('click', (e) => {
            document.querySelectorAll('[data-column]').forEach(b => b.classList.remove('active'));
            e.target.classList.add('active');
            DistillationModule.params.columnType = e.target.dataset.column;
        });
    });
    
    // Range sliders
    ['feed-comp', 'distillate-comp', 'bottoms-comp'].forEach(id => {
        document.getElementById(id)?.addEventListener('input', (e) => {
            e.target.nextElementSibling.textContent = (e.target.value / 100).toFixed(2);
        });
    });
    
    document.getElementById('toggle-calculations')?.addEventListener('click', () => {
        document.getElementById('calculation-steps').classList.toggle('hidden');
    });
}

function runDistillationSimulation() {
    const params = {
        feedComposition: parseFloat(document.getElementById('feed-comp')?.value || 50) / 100,
        distillateComposition: parseFloat(document.getElementById('distillate-comp')?.value || 95) / 100,
        bottomsComposition: parseFloat(document.getElementById('bottoms-comp')?.value || 5) / 100,
        feedCondition: 1.0,
        relativeVolatility: parseFloat(document.getElementById('alpha')?.value || 2.5),
        refluxRatio: parseFloat(document.getElementById('reflux')?.value || 2.0),
        columnType: document.querySelector('[data-column].active')?.dataset.column || 'tray',
        hetp: 0.5,
        feedFlowRate: parseFloat(document.getElementById('feed-flow')?.value || 100),
        trayEfficiency: 0.7
    };
    
    DistillationModule.params = params;
    const results = calculateDistillation(params);
    DistillationModule.results = results;
    
    if (results.error) {
        showToast(results.message, 'error');
        return;
    }
    
    // Update display
    document.getElementById('result-stages').textContent = results.theoreticalStages;
    document.getElementById('result-feed-stage').textContent = results.feedStage;
    document.getElementById('result-rmin').textContent = results.Rmin.toFixed(3);
    document.getElementById('result-r-ratio').textContent = results.RoverRmin.toFixed(2);
    document.getElementById('result-d').textContent = results.materialBalance.D.toFixed(1);
    document.getElementById('result-b').textContent = results.materialBalance.B.toFixed(1);
    
    // Update steps
    const stepsContainer = document.getElementById('calculation-steps');
    stepsContainer.innerHTML = results.calculationSteps.map((step, i) => `
        <div class="calc-step">
            <div class="step-number">${i + 1}</div>
            <div class="step-content">
                <div class="step-description">${step.description}</div>
                <div class="step-formula">${step.formula}</div>
                <div class="step-result">${step.result}</div>
            </div>
        </div>
    `).join('');
    
    showToast('Calcul de distillation terminé', 'success');
}

DistillationModule.getExplanation = function() {
    return {
        title: 'Simulation de Distillation',
        description: `
            La distillation est une opération unitaire de séparation basée sur 
            la différence de volatilité des composants d'un mélange liquide.
        `,
        theory: `
            <p><strong>Méthode de McCabe-Thiele:</strong></p>
            <ul>
                <li>Hypothèse de flux molaires constants (CMO)</li>
                <li>Construction graphique étage par étage</li>
                <li>Droite opératoire de rectification et de stripping</li>
            </ul>
        `,
        formulas: `
            <p><strong>Courbe d'équilibre (volatilité constante):</strong></p>
            $$ y^* = \\frac{\\alpha x}{1 + (\\alpha - 1)x} $$
            
            <p><strong>Droite de rectification:</strong></p>
            $$ y = \\frac{R}{R+1}x + \\frac{x_D}{R+1} $$
        `,
        references: [
            'McCabe, W.L., Thiele, E.W. "Graphical Design of Fractionating Columns" (1925)',
            'Seader, Henley "Separation Process Principles"'
        ]
    };
};

if (typeof registerModule !== 'undefined') {
    registerModule('distillation', DistillationModule);
}

window.DistillationModule = DistillationModule;
