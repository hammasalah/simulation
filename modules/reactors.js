/**
 * ============================================================================
 * CHEMICAL REACTOR SIMULATION MODULE - ADVANCED
 * ============================================================================
 * 
 * Simulations for Batch, CSTR, PFR, and PBR reactors with advanced features:
 * - Multiple Reactions (Series/Parallel)
 * - Non-Isothermal Operation (Energy Balance)
 * - Pressure Drop in PBR (Ergun Equation)
 * - Residence Time Distribution (RTD)
 * 
 * @version 0.2.0-advanced
 */

'use strict';

console.log("Loading reactors.js module...");

var ReactorsModule = {
    name: 'reactors',
    
    defaults: {
        reactorType: 'pfr',         // batch, cstr, pfr, pbr
        reactionScheme: 'AtoB',     // AtoB, AtoBtoC, 2AtoB
        k1: 0.1,                    // Rate constant 1 (1/s or L/mol.s)
        k2: 0.05,                   // Rate constant 2 (series/parallel)
        CA0: 1.0,                   // mol/L
        v0: 10,                     // L/s
        V: 100,                     // Reactor Volume (L) or Weight (kg for PBR)
        
        // Energy
        nonIsothermal: false,
        T0: 300,                    // Inlet Temp (K)
        Ta: 300,                    // Ambient/Coolant Temp (K)
        Ua: 100,                    // Overall Heat Transfer * Area (J/s.K.m3) -> U*a
        dHrxn: -50000,              // Heat of reaction (J/mol) (Exothermic < 0)
        Cp: 4000,                   // Heat capacity (J/kg.K) approx water
        rho: 1000,                  // Density (kg/m3)
        
        // PBR Specific
        alpha: 0.01,                // Pressure drop parameter (1/kg)
        
        // RTD
        sigma: 5                    // Dispersion/Time parameter
    },
    
    params: null,
    results: null
};

// IMMEDIATE REGISTRATION
// This ensures the module is available even if subsequent code (init/solver) fails
if (typeof window !== 'undefined') {
    window.ReactorsModule = ReactorsModule;
    console.log("ReactorsModule registered globally (Safety Check)");
}
if (typeof registerModule !== 'undefined') {
    registerModule('reactors', ReactorsModule);
}

// ============================================================================
// ODE SOLVER (Runge-Kutta 4) - Local Implementation
// ============================================================================
function reactorRK4(dydx, x0, y0, xf, step) {
    let x = x0;
    let y = [...y0]; // State vector
    const history = [{x, y: [...y]}];
    
    const nSteps = Math.ceil((xf - x0) / step);
    const h = (xf - x0) / nSteps; // Adjust step exact fit
    
    for (let i = 0; i < nSteps; i++) {
        const k1 = dydx(x, y);
        const y_k1 = y.map((val, j) => val + 0.5 * h * k1[j]);
        
        const k2 = dydx(x + 0.5 * h, y_k1);
        const y_k2 = y.map((val, j) => val + 0.5 * h * k2[j]);
        
        const k3 = dydx(x + 0.5 * h, y_k2);
        const y_k3 = y.map((val, j) => val + h * k3[j]);
        
        const k4 = dydx(x + h, y_k3);
        
        // Update y
        y = y.map((val, j) => val + (h / 6) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]));
        x += h;
        
        history.push({x, y: [...y]});
    }
    return history;
}

// ============================================================================
// KINETICS & PHYSICS
// ============================================================================

/**
 * Calculates rates of change for PFR/PBR/Batch
 * Returns dy/dVar array.
 * 
 * State Vector Map:
 * pfr/pbr: [Fa, Fb, Fc, T, y]  (y = P/P0)
 * batch:   [Ca, Cb, Cc, T]
 */
function getDerivatives(type, volumeOrTime, state, params) {
    const { 
        reactionScheme, k1, k2, CA0, v0, 
        nonIsothermal, dHrxn, Ua, Ta, T0, Cp, rho,
        alpha, reactorType 
    } = params;

    // Unpack state
    // Common: Concentrations need to be calculated first
    let Fa, Fb, Fc, T, P_ratio;
    let Ca, Cb, Cc;
    
    // Arrhenius correction
    // k(T) = k(T0) * exp(E/R * (1/T0 - 1/T)) 
    // Simplified: Just scaling k linearly or using provided Arrhenius params if we had them.
    // For this demo, let's assume k input is at T0, and strictly correct using approximation
    // Let's implement real Arrhenius if E is provided, else assume small T effect or input k is at T0
    // Using simple doubling per 10K idea or explicit E/R? 
    // Let's use E/R = 5000 (approx Ea=40kJ) for demo sensitivity
    const E_R = 5000; 
    
    // Arrhenius Helper
    const arrheniusTerm = (T_current) => Math.exp((-E_R) * (1/T_current - 1/T0));

    const currentT = state[3] || T0;

    // Dynamic rates
    const kr1 = k1 * (nonIsothermal ? arrheniusTerm(currentT) : 1);
    const kr2 = k2 * (nonIsothermal ? arrheniusTerm(currentT) : 1);

    if (type === 'pfr' || type === 'pbr') {
        Fa = state[0];
        Fb = state[1];
        Fc = state[2];
        T = state[3];
        P_ratio = state[4];

        const F_total = Fa + Fb + Fc;
        const P = P_ratio; // relative pressure
        // v = v0 * (F/F0) * (P0/P) * (T/T0)
        // Gas Phase correction needed? Assume Liquid for CSTR/Batch/PFR simpliciy unless PBR gas
        
        // Let's assume Liquid Phase (constant density) for PFR/CSTR default
        // PBR usually Gas. 
        // Let's simplify: Constant Volumetric Flow v = v0 for Liquid PFR.
        // For PBR Gas: v = v0 * (P0/P) * (T/T0) (assuming F_total constant or neglected mole change)
        
        let v = v0;
        if (reactorType === 'pbr') {
            v = v0 * (1/P_ratio) * (T/T0);
        }

        Ca = Fa / v;
        Cb = Fb / v;
        Cc = Fc / v;

        // Rate Laws
        let rA = 0, rB = 0, rC = 0;
        
        if (reactionScheme === 'AtoB') {
            // A -> B
            rA = -kr1 * Ca;
            rB = kr1 * Ca;
        } else if (reactionScheme === 'AtoBtoC') {
            // A -> B -> C
            rA = -kr1 * Ca;
            rB = kr1 * Ca - kr2 * Cb;
            rC = kr2 * Cb;
        } else if (reactionScheme === '2AtoB') {
            // 2A -> B
            rA = -kr1 * Ca * Ca;
            rB = 0.5 * kr1 * Ca * Ca;
        }

        // Mass Balances (dF/dV = r)
        // For PBR dF/dW = r' (rate per kg)
        // We will treat V input as W for PBR, so equations same form just param change
        
        const dFadV = rA; 
        const dFbdV = rB;
        const dFcdV = rC;

        // Energy Balance
        // dT/dV = (Q_generated - Q_removed) / (Sum(Fi*Cpi))
        // Q_gen = (-dHr) * (-rA)  (approx for single reaction A->B)
        // If complex: Sum(r_i * dH_i)? 
        // Simplified: dT/dV = ( (Ua*(Ta - T)) + (-dHrxn * (-rA)) ) / (F_total * Cp)
        // Note: Cp is usually molar here if F is molar. Input Cp is J/kg.K.
        // Let's convert: m_dot * Cp_mass ~ F_total * MW_avg * Cp_mass?
        // Let's assume parameters are consistent. F*Cp term represents Heat Capacity Flow rate.
        // F_total (mol/s) * 40 (J/mol.K) approx? Let's use simplified m_dot*Cp
        const m_dot_Cp = v0 * rho * Cp / 1000; // J/s.K approx fixed
        
        let dTdV = 0;
        if (nonIsothermal) {
             const Q_gen = (-dHrxn) * (kr1 * Ca); // Only A->B heat considered for demo
             const Q_rem = Ua * (T - Ta); 
             dTdV = (Q_gen + Q_rem) / m_dot_Cp;
        }

        // Pressure Drop (Ergun simplified)
        // dy/dW = -alpha/2y * (T/T0)
        let dydV = 0;
        if (reactorType === 'pbr') {
            dydV = -(alpha / (2 * P_ratio)) * (T / T0);
        }

        return [dFadV, dFbdV, dFcdV, dTdV, dydV];
    
    } else if (type === 'batch') {
        // [Ca, Cb, Cc, T]
        Ca = state[0];
        Cb = state[1];
        Cc = state[2];
        T = state[3];
        
        let rA = 0, rB = 0, rC = 0;
        // Same Rate Laws ...
        if (reactionScheme === 'AtoB') {
             rA = -kr1 * Ca; rB = kr1 * Ca;
        } else if (reactionScheme === 'AtoBtoC') {
             rA = -kr1 * Ca; rB = kr1 * Ca - kr2 * Cb; rC = kr2 * Cb;
        } else if (reactionScheme === '2AtoB') {
             rA = -kr1 * Ca * Ca; rB = 0.5 * kr1 * Ca * Ca;
        }

        // Energy Batch: dT/dt = (Q_gen + Q_exchange) / (V * rho * Cp)
        // V is fixed batch vol
        // Heat Capacity = V * rho * Cp
        const HeatCap = params.V * rho * Cp / 1000; // J/K
        let dTdt = 0;
        if (nonIsothermal) {
             const Q_gen = (-dHrxn) * (-rA) * params.V; // J/s
             const Q_rem = Ua * params.V * (Ta - T); // Ua here per unit volume? Params Ua is defined as U*a (energy/vol). 
             // Normally U*Area. Let's assume input Ua is U*Area for batch or U*a*V. 
             // Let's assume param Ua is "UA" total for batch, but "Ua" (per vol) for PFR?
             // To keep code simple, let's assume param Ua value is scaled to total Watts/K.
             dTdt = (Q_gen + Ua * (Ta - T)) / HeatCap; 
        }

        return [rA, rB, rC, dTdt, 0];
    }
}

// ============================================================================
// SIMULATION ENGINE
// ============================================================================

function solveSystem(params) {
    const { reactorType, CA0, V, v0, T0 } = params;
    
    // Initial Conditions
    // PFR/PBR: Flow rates [Fa0, Fb0, Fc0, T0, y0=1]
    const FA0 = CA0 * v0;
    const initialPFR = [FA0, 0, 0, T0, 1.0];
    
    // Batch: Concentrations [Ca0, Cb0, Cc0, T0]
    const initialBatch = [CA0, 0, 0, T0, 0];

    // CSTR is algebraic, handle separately
    if (reactorType === 'cstr') {
        return solveCSTR_Algebraic(params);
    }
    
    // Run Solver
    const isFlow = (reactorType === 'pfr' || reactorType === 'pbr');
    const y0 = isFlow ? initialPFR : initialBatch;
    const endVar = isFlow ? V : 100; // Volume for flow, Time for batch (fixed 100s or calculated?)
    // For batch we simulate time. How much time? Until 99% conv or fixed limit.
    // Let's do fixed reasonable time based on k1. tau ~ 3/k1
    const tEnd = (reactorType === 'batch') ? (4.0 / params.k1) : V; 
    
    const results = reactorRK4(
        (x, y) => getDerivatives(isFlow ? reactorType : 'batch', x, y, params),
        0, y0, tEnd, tEnd/100
    );

    // Process Results
    return results.map(step => {
        const x = step.x; // Time or Vol
        const y = step.y;
        
        let Ca, Cb, Cc, T, P;
        
        if (isFlow) {
            T = y[3];
            let v_local = v0;
             // Correction for PBR density change
            if (reactorType === 'pbr') {
                v_local = v0 * (1/y[4]) * (T/params.T0);
            }
            Ca = y[0] / v_local;
            Cb = y[1] / v_local;
            Cc = y[2] / v_local;
            P = y[4];
        } else {
            Ca = y[0]; Cb = y[1]; Cc = y[2]; T = y[3]; P = 1;
        }

        return {
             metric: x, // Time or Volume
             Ca: Math.max(0, Ca),
             Cb: Math.max(0, Cb),
             Cc: Math.max(0, Cc),
             T: T,
             P: P,
             Conversion: (CA0 - Ca)/CA0
        };
    });
}

function solveCSTR_Algebraic(params) {
    // CSTR Single Point calculation
    // Algebraic V = FA0 * X / -rA
    // Must iterate if non-isothermal or complex kinetics
    // Simplified: Isothermal CSTR A->B
    const { CA0, v0, V, k1, reactionScheme } = params;
    const tau = V / v0;
    
    // 1st Order A->B: Ca = Ca0 / (1 + k*tau)
    let Ca = CA0 / (1 + k1 * tau);
    let Cb = CA0 - Ca;
    let Cc = 0;
    
    if (reactionScheme === '2AtoB') {
        // 2nd Order: Ca0 - Ca = k*tau*Ca^2 -> k*tau*Ca^2 + Ca - Ca0 = 0
        // Quad formula
        const a = k1 * tau;
        const b = 1;
        const c = -CA0;
        Ca = (-b + Math.sqrt(b*b - 4*a*c))/(2*a);
        Cb = 0.5 * (CA0 - Ca);
    } else if (reactionScheme === 'AtoBtoC') {
        // A -> B -> C
        // Ca = Ca0 / (1 + k1*tau)
        // Cb = (k1*Ca*tau) / (1 + k2*tau)
        Ca = CA0 / (1 + k1 * tau);
        Cb = (k1 * Ca * tau) / (1 + params.k2 * tau);
        Cc = CA0 - Ca - Cb;
    }
    
    return [{
        metric: V,
        Ca, Cb, Cc, 
        T: params.T0, // Adiabatic CSTR not implemented in this simplest algebraic solver
        P: 1,
        Conversion: (CA0 - Ca)/CA0
    }];
}

// ============================================================================
// RTD FUNCTIONS
// ============================================================================
function calculateRTD(params) {
    // E-curve generation
    const { reactorType, v0, V, sigma } = params;
    const tau = V / v0;
    const data = [];
    
    const tMax = 3 * tau;
    const dt = tMax / 100;
    
    for(let t=0; t <= tMax; t+=dt) {
        let E = 0;
        if (reactorType === 'cstr') {
            // Ideal CSTR: E(t) = (1/tau) * exp(-t/tau)
            E = (1/tau) * Math.exp(-t/tau);
        } else if (reactorType === 'pfr' || reactorType === 'pbr') {
            // Dispersion Model (Gaussian approx)
            // E(t) = 1/(2*sqrt(pi*D/uL)*t) ... simplified normal distribution around tau
            const sig = sigma || (tau * 0.1); // width
            E = (1 / (sig * Math.sqrt(2*Math.PI))) * Math.exp( -0.5 * Math.pow((t-tau)/sig, 2) );
        } else {
             // Batch: No flow, RTD undefined (delta at t_residence)
             E = (t > tau-dt && t < tau+dt) ? 1/dt : 0;
        }
        data.push({t, E});
    }
    return data;
}

// ============================================================================
// UI
// ============================================================================

ReactorsModule.render = function() {
    this.params = { ...this.defaults };
    
    return `
        <style>
            .d-flex { display: flex; }
            .justify-between { justify-content: space-between; }
            .align-center { align-items: center; }
            .mb-5 { margin-bottom: 1.5rem; }
            .mb-4 { margin-bottom: 1rem; }
            .mb-3 { margin-bottom: 0.75rem; }
            .mt-4 { margin-top: 1rem; }
            .gap-2 { gap: 0.5rem; }
            .p-2 { padding: 0.5rem; }
            .pt-3 { padding-top: 0.75rem; }
            .pb-0 { padding-bottom: 0; }
            .w-full { width: 100%; }
            .text-sm { font-size: 0.875rem; }
            .text-muted { color: #6c757d; }
            .font-bold { font-weight: bold; }
            
            .param-grid.two-col { display: grid; grid-template-columns: 1fr 1fr; gap: 10px; }
            .tabs { display: flex; gap: 5px; margin-bottom: 10px; }
            .tab-btn { padding: 8px 15px; border: none; background: transparent; cursor: pointer; border-bottom: 2px solid transparent; transition: all 0.2s; }
            .tab-btn.active { border-bottom-color: var(--primary, #008080); color: var(--primary, #008080); font-weight: bold; }
            .tab-btn:hover { background: rgba(0,0,0,0.05); }
            
            .mode-selector { display: flex; gap: 5px; background: #eee; padding: 4px; border-radius: 6px; }
            .mode-btn { flex: 1; padding: 6px; border: none; background: transparent; border-radius: 4px; cursor: pointer; transition: all 0.2s; }
            .mode-btn.active { background: white; box-shadow: 0 1px 3px rgba(0,0,0,0.1); font-weight: bold; color: var(--primary, #008080); }
            
            .results-grid { display: grid; grid-template-columns: repeat(2, 1fr); gap: 10px; }
            .result-item { background: #f8f9fa; padding: 10px; border-radius: 6px; text-align: center; border: 1px solid #eee; }
            .result-value { font-size: 1.5em; font-weight: bold; color: var(--primary, #008080); margin: 5px 0; }
            .result-label { font-size: 0.8em; text-transform: uppercase; color: #666; letter-spacing: 0.5px; }
            
            .btn-outline { border: 1px solid #ddd; background: white; border-radius: 4px; padding: 4px 8px; cursor: pointer; }
            .btn-outline:hover { background: #f8f9fa; border-color: #ccc; }
        </style>
        <div class="reactors-module" id="reactors-content">
            <div class="d-flex justify-between align-center mb-5">
                <div>
                    <h2>Simulation de Réacteurs Chimiques v2.0</h2>
                    <p class="text-muted">Multi-réactions, Thermique & RTD</p>
                </div>
                <div class="status-badge success">
                    <span class="status-dot"></span>
                    <span>Modèle Avancé</span>
                </div>
            </div>

            <div style="display: grid; grid-template-columns: 1fr 1.5fr; gap: 20px;">
                <!-- CONTROLS -->
                <div>
                    <!-- Type Selector -->
                    <div class="card mb-4">
                        <div class="card-body p-2">
                             <div class="mode-selector text-sm">
                                <button class="mode-btn active" onclick="setRType('cstr')">CSTR</button>
                                <button class="mode-btn" onclick="setRType('pfr')">PFR</button>
                                <button class="mode-btn" onclick="setRType('pbr')">PBR</button>
                                <button class="mode-btn" onclick="setRType('batch')">Batch</button>
                            </div>
                        </div>
                    </div>

                    <div class="card">
                        <div class="card-header pb-0">
                            <div class="tabs" style="border-bottom: 1px solid #eee;">
                                <button class="tab-btn active" onclick="switchRTab('kinetics')" id="tab-k">Cinétique</button>
                                <button class="tab-btn" onclick="switchRTab('energy')" id="tab-e">Énergie</button>
                                <button class="tab-btn" onclick="switchRTab('phys')" id="tab-p">Physique</button>
                            </div>
                        </div>
                        <div class="card-body pt-3">
                            
                            <!-- KINETICS TAB -->
                            <div id="panel-kinetics">
                                <div class="form-group">
                                    <label>Schéma Réactionnel</label>
                                    <select class="form-select" id="rxn-scheme">
                                        <option value="AtoB">Simple A → B</option>
                                        <option value="AtoBtoC">Série A → B → C</option>
                                        <option value="2AtoB">Ordre 2 (2A → B)</option>
                                    </select>
                                </div>
                                <div class="param-grid two-col">
                                    <div class="form-group">
                                        <label>k1 (s⁻¹)</label>
                                        <input type="number" class="form-input" id="k1" value="0.1" step="0.01">
                                    </div>
                                    <div class="form-group" id="grp-k2" style="display:none;">
                                        <label>k2 (s⁻¹)</label>
                                        <input type="number" class="form-input" id="k2" value="0.05" step="0.01">
                                    </div>
                                    <div class="form-group">
                                        <label>CA0 (mol/L)</label>
                                        <input type="number" class="form-input" id="ca0" value="1.0">
                                    </div>
                                </div>
                            </div>

                            <!-- ENERGY TAB -->
                            <div id="panel-energy" style="display:none;">
                                <div class="form-group mb-3">
                                    <label class="form-check">
                                        <input type="checkbox" id="non-iso">
                                        <strong>Mode Non-Isotherme</strong>
                                    </label>
                                </div>
                                <div class="param-grid two-col">
                                    <div class="form-group">
                                        <label>T entrée (K)</label>
                                        <input type="number" class="form-input" id="t0" value="300">
                                    </div>
                                    <div class="form-group">
                                        <label>T fluide (K)</label>
                                        <input type="number" class="form-input" id="ta" value="300">
                                    </div>
                                    <div class="form-group">
                                        <label>ΔHrxn (J/mol)</label>
                                        <input type="number" class="form-input" id="dhrxn" value="-50000">
                                    </div>
                                    <div class="form-group">
                                        <label>UA (W/K)</label>
                                        <input type="number" class="form-input" id="ua" value="500">
                                    </div>
                                </div>
                            </div>

                            <!-- PHYSICS TAB -->
                            <div id="panel-phys" style="display:none;">
                                <div class="form-group">
                                    <label>Volume/Poids (L ou kg)</label>
                                    <input type="number" class="form-input" id="vol" value="100">
                                </div>
                                <div class="form-group">
                                    <label>Débit v0 (L/s)</label>
                                    <input type="number" class="form-input" id="v0" value="10">
                                </div>
                                <div id="pbr-params" style="display:none;" class="mt-3 border-t pt-2">
                                    <label class="text-sm font-bold text-muted">PBR Spécifique</label>
                                    <div class="form-group">
                                        <label>Perte Charge α (kg⁻¹)</label>
                                        <input type="number" class="form-input" id="alpha" value="0.01" step="0.001">
                                    </div>
                                </div>
                            </div>

                            <button class="btn btn-primary w-full mt-4" onclick="runRxnSim()">Lancer Simulation</button>
                        </div>
                    </div>
                </div>

                <!-- OUTPUTS -->
                <div>
                    <!-- KPI Cards -->
                    <div class="results-grid mb-3">
                        <div class="result-item">
                            <div class="result-label">Conversion</div>
                            <div class="result-value" id="res-x">--</div>
                            <div class="result-unit">%</div>
                        </div>
                        <div class="result-item">
                            <div class="result-label">Selectivité B</div>
                            <div class="result-value" id="res-sel">--</div>
                            <div class="result-unit">%</div>
                        </div>
                        <div class="result-item">
                            <div class="result-label">Sortie T</div>
                            <div class="result-value" id="res-t">--</div>
                            <div class="result-unit">K</div>
                        </div>
                        <div class="result-item">
                            <div class="result-label">Chute Pression</div>
                            <div class="result-value" id="res-dp">--</div>
                            <div class="result-unit">%</div>
                        </div>
                    </div>

                    <!-- Charts TABS -->
                    <div class="d-flex gap-2 mb-2">
                        <button class="btn-sm btn-outline" onclick="showChart('conc')">Concentration</button>
                        <button class="btn-sm btn-outline" onclick="showChart('temp')">Température</button>
                        <button class="btn-sm btn-outline" onclick="showChart('rtd')">DTS (RTD)</button>
                    </div>
                    
                    <div class="card" style="height: 300px;">
                        <div class="card-body" id="chart-area" style="height: 100%; position: relative;">
                            <canvas id="rxnCanvas"></canvas>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    `;
};

// ============================================================================
// CLIENT SIDE LOGIC
// ============================================================================

ReactorsModule.init = function() {
    console.log("ReactorsModule.init() called");
    
    try {
        window.currentRType = 'cstr';
        
        // Helper exports
        window.setRType = (type) => {
            window.currentRType = type;
            document.querySelectorAll('.mode-btn').forEach(b => {
                b.classList.toggle('active', b.innerText.toLowerCase() === type || (b.innerText==='Batch' && type==='batch'));
            });
            const pbrEl = document.getElementById('pbr-params');
            if (pbrEl) pbrEl.style.display = type === 'pbr' ? 'block' : 'none';
            
            // Update labels based on type
            const volLabel = document.querySelector('label[for="vol"]');
            if(volLabel) volLabel.innerText = type === 'pbr' ? 'Poids Catalyseur (kg)' : 'Volume Réacteur (L)';
            
            if (typeof runRxnSim === 'function') runRxnSim();
        };

        window.switchRTab = (tab) => {
            ['kinetics', 'energy', 'phys'].forEach(t => {
                const panel = document.getElementById('panel-' + t);
                const tabBtn = document.getElementById('tab-' + t[0]);
                if (panel) panel.style.display = t === tab ? 'block' : 'none';
                if (tabBtn) tabBtn.classList.toggle('active', t === tab);
            });
        };

        const rxnScheme = document.getElementById('rxn-scheme');
        if (rxnScheme) {
            rxnScheme.addEventListener('change', (e) => {
                const grpK2 = document.getElementById('grp-k2');
                if (grpK2) grpK2.style.display = e.target.value === 'AtoBtoC' ? 'block' : 'none';
            });
        }
    
    // QuickChart Helper
    window.showChart = (type) => {
        if(!ReactorsModule.lastResults) return;
        const res = ReactorsModule.lastResults;
        const ctx = document.getElementById('rxnCanvas').getContext('2d');
        const W = ctx.canvas.width = ctx.canvas.parentElement.clientWidth;
        const H = ctx.canvas.height = ctx.canvas.parentElement.clientHeight;
        
        ctx.clearRect(0,0,W,H);
        
        // Simple Plotter
        const padding = 40;
        const datas = []; 
        let xLabel = "";
        
        if (type === 'conc') {
            xLabel = "Distance / Temps";
            datas.push({l: 'Ca', c: '#e74c3c', pts: res.map(p => ({x: p.metric, y: p.Ca}))});
            datas.push({l: 'Cb', c: '#3498db', pts: res.map(p => ({x: p.metric, y: p.Cb}))});
            if (res[0].Cc !== undefined) datas.push({l: 'Cc', c: '#2ecc71', pts: res.map(p => ({x: p.metric, y: p.Cc}))});
        } else if (type === 'temp') {
            xLabel = "Distance / Temps";
            datas.push({l: 'Temp (K)', c: '#e67e22', pts: res.map(p => ({x: p.metric, y: p.T}))});
        } else if (type === 'rtd') {
            xLabel = "Temps (s)";
            const rtd = calculateRTD(ReactorsModule.params);
            datas.push({l: 'E(t)', c: '#9b59b6', pts: rtd.map(p => ({x: p.t, y: p.E}))});
        }

        // Draw Axes
        const xMax = datas[0].pts[datas[0].pts.length-1].x;
        const yMax = Math.max(...datas.flatMap(d => d.pts.map(p => p.y))) * 1.1;
        
        const sx = x => padding + (x/xMax)*(W-2*padding);
        const sy = y => H - padding - (y/yMax)*(H-2*padding);
        
        ctx.beginPath();
        ctx.strokeStyle = '#666';
        ctx.moveTo(padding, 0); ctx.lineTo(padding, H-padding); ctx.lineTo(W, H-padding);
        ctx.stroke();

        // Draw Lines
        datas.forEach(d => {
            ctx.beginPath();
            ctx.strokeStyle = d.c;
            ctx.lineWidth = 2;
            d.pts.forEach((p, i) => {
                if(i===0) ctx.moveTo(sx(p.x), sy(p.y));
                else ctx.lineTo(sx(p.x), sy(p.y));
            });
            ctx.stroke();
            // Legend
            ctx.fillStyle = d.c;
            ctx.fillText(d.l, sx(d.pts[d.pts.length-1].x)-30, sy(d.pts[d.pts.length-1].y)-10);
        });
    };

    window.runRxnSim = () => {
        try {
            const rxnSchemeEl = document.getElementById('rxn-scheme');
            const k1El = document.getElementById('k1');
            const k2El = document.getElementById('k2');
            const ca0El = document.getElementById('ca0');
            const t0El = document.getElementById('t0');
            const taEl = document.getElementById('ta');
            const dhrxnEl = document.getElementById('dhrxn');
            const uaEl = document.getElementById('ua');
            const volEl = document.getElementById('vol');
            const v0El = document.getElementById('v0');
            const alphaEl = document.getElementById('alpha');
            const nonIsoEl = document.getElementById('non-iso');
            
            // Guard against missing elements
            if (!rxnSchemeEl || !k1El || !ca0El || !volEl || !v0El) {
                console.warn("runRxnSim: Required DOM elements not found yet");
                return;
            }
            
            const params = {
                reactorType: window.currentRType || 'cstr',
                reactionScheme: rxnSchemeEl.value,
                k1: parseFloat(k1El.value) || 0.1,
                k2: parseFloat(k2El?.value) || 0.05,
                CA0: parseFloat(ca0El.value) || 1.0,
                T0: parseFloat(t0El?.value) || 300,
                Ta: parseFloat(taEl?.value) || 300,
                dHrxn: parseFloat(dhrxnEl?.value) || -50000,
                Ua: parseFloat(uaEl?.value) || 100,
                V: parseFloat(volEl.value) || 100,
                v0: parseFloat(v0El.value) || 10,
                alpha: parseFloat(alphaEl?.value) || 0.01,
                nonIsothermal: nonIsoEl?.checked || false,
                Cp: 4000, rho: 1000
            };
            
            ReactorsModule.params = params;
            const results = solveSystem(params);
            ReactorsModule.lastResults = results;
            
            // Stats
            const final = results[results.length-1];
            const resX = document.getElementById('res-x');
            const resSel = document.getElementById('res-sel');
            const resT = document.getElementById('res-t');
            const resDp = document.getElementById('res-dp');
            
            if (resX) resX.innerText = (final.Conversion*100).toFixed(1);
            const sel = final.Cb / (final.Ca + final.Cb + final.Cc + 0.0001) * 100;
            if (resSel) resSel.innerText = sel.toFixed(1);
            if (resT) resT.innerText = final.T.toFixed(1);
            if (resDp) resDp.innerText = ((1 - final.P)*100).toFixed(1);
            
            window.showChart('conc');
            if(window.showToast) showToast('Simulation Terminée', 'success');
        } catch (err) {
            console.error("runRxnSim error:", err);
        }
    };

    // Initial Run
    setTimeout(runRxnSim, 200);
    
    // Default PBR hide
    const pbrParams = document.getElementById('pbr-params');
    if(pbrParams) pbrParams.style.display = 'none';
    
    console.log("ReactorsModule.init() completed");
    
    } catch (initError) {
        console.error("ReactorsModule.init() error:", initError);
    }
};

// Register the module
if (typeof registerModule === 'function') {
    registerModule('reactors', ReactorsModule);
} else if (window.ModuleRegistry) {
    window.ModuleRegistry['reactors'] = ReactorsModule;
}
window.ReactorsModule = ReactorsModule;

