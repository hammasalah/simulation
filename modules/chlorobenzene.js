/**
 * ============================================================================
 * CHLOROBENZENE SYNTHESIS SIMULATION
 * ============================================================================
 * 
 industrial synthesis of chlorobenzene via electrophilic aromatic substitution.
 * Supports:
 * - Series CSTRs (1 to 5)
 * - Kinetic Optimization
 * - Detailed Heat Balance
 * - Downstream Separation Logic
 */

'use strict';

const ChlorobenzeneModule = {
    name: 'chlorobenzene',
    
    defaults: {
        benzeneFlow: 100,      // F_A0 (mol/min)
        chlorineFlow: 80,      // F_B0 (mol/min)
        temperature: 55,       // T (°C)
        pressure: 1,           // P (bar)
        reactorVolume: 2000,   // V (L) total
        numberOfReactors: 3,   // N
        
        // Kinetic Parameters
        k1_ref: 0.52,          // L/mol.min @ 55°C (Benzene -> Mono)
        k2_ref: 0.08,          // L/mol.min @ 55°C (Mono -> Di)
        Ea1: 45000,            // J/mol
        Ea2: 52000             // J/mol
    },

    // Physical Properties
    props: {
        benzene: { Mw: 78.11, rho: 876, Cp: 136, name: 'Benzène', bp: 80.1 },
        chlorine: { Mw: 70.9, rho: 1400, Cp: 65, name: 'Chlore', bp: -34 }, 
        monochlor: { Mw: 112.56, rho: 1106, Cp: 145, name: 'Chlorobenzène', bp: 131.6 },
        dichlor: { Mw: 147.0, rho: 1300, Cp: 170, name: 'Dichlorobenzène', bp: 180.5 },
        hcl: { Mw: 36.46, rho: 1000, Cp: 29, name: 'HCl', bp: -85 },
        
        dH_rxn1: -131800, // J/mol (Benzo + Cl2 -> Mono + HCl)
        dH_rxn2: -125000  // J/mol (Mono + Cl2 -> Di + HCl)
    }
};

// ============================================================================
// SIMULATION ENGINE
// ============================================================================

/**
 * Arrhenius Law
 */
function getK(T_K, k_ref, Ea) {
    const R = 8.314;
    const T_ref = 273.15 + 55; // Ref temp 55°C
    return k_ref * Math.exp((-Ea / R) * (1/T_K - 1/T_ref));
}

/**
 * Solve a Single CSTR via Iteration (Newton-ish on Chlorine Concentration)
 * Inputs: Fin (mol/min), V (L), k1, k2
 */
function solveSingleCSTR(Fin, V, k1, k2, v0) {
    // Inlet Concentrations
    const Cai = Fin.benzene / v0;
    const Cbi = Fin.chlorine / v0;
    const Cpi = Fin.monochlor / v0;
    const Cqi = Fin.dichlor / v0;
    const Chi = Fin.hcl / v0;

    const tau = V / v0;

    // Reactions:
    // 1) A + B -> P + H  (r1 = k1*Ca*Cb)
    // 2) P + B -> Q + H  (r2 = k2*Cp*Cb)
    
    // Algebraic relations from mass balance:
    // Ca = Cai / (1 + k1*tau*Cb)
    // Cp = (Cpi + k1*tau*Ca*Cb) / (1 + k2*tau*Cb)
    // Cb = Cbi - tau*(r1 + r2)  <-- Implicit eq to solve
    
    // Geometric search for Cb in [0, Cbi]
    let Cb_guess = Cbi * 0.5;
    let err = 1;
    let iter = 0;
    
    let Ca_ss = 0, Cp_ss = 0;
    
    while(err > 1e-6 && iter < 100) {
        // Calculate dependent concentrations based on current Cb_guess
        const denom1 = (1 + k1 * tau * Cb_guess);
        const Ca = Cai / denom1;
        
        const denom2 = (1 + k2 * tau * Cb_guess);
        const term1 = k1 * tau * Ca * Cb_guess;
        const Cp = (Cpi + term1) / denom2;
        
        // Calculate Rates
        const r1 = k1 * Ca * Cb_guess;
        const r2 = k2 * Cp * Cb_guess;
        
        // Check Chlorine Balance
        // F_B_out = F_B_in - V(r1 + r2)
        // Cb_calc = Cbi - tau*(r1 + r2)
        const Cb_calc = Cbi - tau * (r1 + r2);
        
        // Update Cb (Damped substitution)
        const diff = Cb_calc - Cb_guess;
        Cb_guess = Cb_guess + 0.5 * diff;
        
        // Constraints
        if (Cb_guess < 0) Cb_guess = 1e-9;
        if (Cb_guess > Cbi) Cb_guess = Cbi;
        
        err = Math.abs(diff);
        iter++;
        
        Ca_ss = Ca;
        Cp_ss = Cp;
    }
    
    // Final concentrations of others
    const r1_final = k1 * Ca_ss * Cb_guess;
    const r2_final = k2 * Cp_ss * Cb_guess;
    
    const Cq_ss = Cqi + tau * r2_final;
    const Ch_ss = Chi + tau * (r1_final + r2_final);

    return {
        concentrations: {
            benzene: Ca_ss, chlorine: Cb_guess, monochlor: Cp_ss, dichlor: Cq_ss, hcl: Ch_ss
        },
        flows: {
            benzene: Ca_ss * v0, chlorine: Cb_guess * v0, monochlor: Cp_ss * v0, dichlor: Cq_ss * v0, hcl: Ch_ss * v0
        },
        rates: { r1: r1_final, r2: r2_final }
    };
}


/**
 * Separation System Simulation
 * Simplified multicomponent split logic
 */
function solveSeparationSystem(reactorOutlet) {
    const flows = reactorOutlet.flows;
    
    // Unit 1: Flash / Decanter (Remove volatiles: HCl, Cl2)
    // Simplifying: 100% HCl removal, 95% Cl2 removal.
    const stream1_Gas = {
        hcl: flows.hcl * 0.99,
        chlorine: flows.chlorine * 0.95,
        benzene: flows.benzene * 0.02 // Loss
    };
    
    const stream2_Liquid = {
        hcl: flows.hcl * 0.01,
        chlorine: flows.chlorine * 0.05,
        benzene: flows.benzene * 0.98,
        monochlor: flows.monochlor,
        dichlor: flows.dichlor
    };
    
    // Unit 2: Benzene Column (Recycle)
    // Separates Benzene (Light) from Chlorobenzenes (Heavy)
    const stream3_Recycle = {
        benzene: stream2_Liquid.benzene * 0.99,
        chlorine: stream2_Liquid.chlorine, 
        monochlor: stream2_Liquid.monochlor * 0.01
    };
    
    const stream4_Bottoms = {
        benzene: stream2_Liquid.benzene * 0.01,
        monochlor: stream2_Liquid.monochlor * 0.99,
        dichlor: stream2_Liquid.dichlor
    };
    
    // Unit 3: Product Column
    // Separates Mono (131°C) from Di (180°C)
    const stream5_Product = {
        benzene: stream4_Bottoms.benzene,
        monochlor: stream4_Bottoms.monochlor * 0.995,
        dichlor: stream4_Bottoms.dichlor * 0.02
    };
    
    const stream6_Heavies = {
        monochlor: stream4_Bottoms.monochlor * 0.005,
        dichlor: stream4_Bottoms.dichlor * 0.98
    };
    
    return {
        gas: stream1_Gas,
        recycle: stream3_Recycle,
        product: stream5_Product,
        waste: stream6_Heavies
    };
}

/**
 * Main Calculation Engine
 */
function calculateChlorobenzeneProcess(params) {
    // 1. Setup Flows and Constants
    const T_K = params.temperature + 273.15;
    const k1 = getK(T_K, params.k1_ref, params.Ea1);
    const k2 = getK(T_K, params.k2_ref, params.Ea2);
    
    const N = params.numberOfReactors;
    const V_i = params.reactorVolume / N; 
    
    // Initial volumetric flow estimate
    const massFlowBenzene = params.benzeneFlow * ChlorobenzeneModule.props.benzene.Mw;
    const massFlowChlorine = params.chlorineFlow * ChlorobenzeneModule.props.chlorine.Mw;
    const rho_mix = 880; // g/L
    const v0 = (massFlowBenzene + massFlowChlorine) / rho_mix; // L/min
    
    let currentFlows = {
        benzene: params.benzeneFlow,
        chlorine: params.chlorineFlow,
        monochlor: 0,
        dichlor: 0,
        hcl: 0
    };
    
    const reactorProfiles = [];
    
    // 2. Solve Series Reactors
    for(let i=0; i<N; i++) {
        const sol = solveSingleCSTR(currentFlows, V_i, k1, k2, v0);
        
        // Heat duty calculation: Q = - V * Rate * dH
        const q_gen_kW = (V_i * (sol.rates.r1 * ChlorobenzeneModule.props.dH_rxn1 + sol.rates.r2 * ChlorobenzeneModule.props.dH_rxn2)) / (60 * 1000); 
        
        reactorProfiles.push({
            id: i+1,
            flows: sol.flows,
            concs: sol.concentrations,
            rates: sol.rates,
            duty: -q_gen_kW // kW cooling required
        });
        
        currentFlows = sol.flows;
    }
    
    // 3. Separation System
    const separation = solveSeparationSystem({ flows: currentFlows });
    
    // 4. Performance Metrics
    const conversionBenzene = (params.benzeneFlow - currentFlows.benzene) / params.benzeneFlow * 100;
    
    const molMono = currentFlows.monochlor;
    const molDi = currentFlows.dichlor;
    const selectivity = molMono > 0 ? molMono / (molMono + molDi) * 100 : 0;
    
    const purity = separation.product.monochlor / (separation.product.monochlor + separation.product.benzene + separation.product.dichlor) * 100;
    
    return {
        reactors: reactorProfiles,
        separation: separation,
        metrics: {
            conversionBenzene,
            selectivity,
            purity
        },
        inputs: params
    };
}

// ============================================================================
// UI RENDERING
// ============================================================================

ChlorobenzeneModule.render = function() {
    return `
    <div class="chlorobenzene-module">
        <!-- Header Section -->
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: var(--space-6);">
            <div>
                <h2 style="margin-bottom: var(--space-1);">Synthèse Industrielle du Chlorobenzène</h2>
                <div style="color: var(--text-secondary); font-size: var(--font-size-sm);">
                    Chloration du Benzène • Réacteurs CSTR en Série • Séparation Intégrée
                </div>
            </div>
            <div style="padding: 0.5rem 1rem; background-color: rgba(39, 174, 96, 0.1); color: var(--success); border-radius: var(--radius-full); font-size: var(--font-size-xs); font-weight: 600; border: 1px solid rgba(39, 174, 96, 0.2);">
                ● Simulation Active
            </div>
        </div>

        <div style="display: grid; grid-template-columns: 1fr 2fr; gap: var(--space-6);">
            
            <!-- Left Panel: Controls -->
            <div style="display: flex; flex-direction: column; gap: var(--space-4);">
                
                <!-- Operating Parameters -->
                <div class="card">
                    <div class="card-header">
                        <h3>Paramètres Opératoires</h3>
                    </div>
                    <div class="card-body">
                        <div class="form-group">
                            <label class="form-label" style="display: flex; justify-content: space-between;">
                                <span>Débit Benzène</span>
                                <span style="font-family: var(--font-mono); color: var(--primary);" id="bz-val-disp">100</span>
                            </label>
                            <input type="range" id="bz-flow" min="10" max="250" value="100" class="form-input" style="padding: 0;">
                            <div class="form-text">mol/min</div>
                        </div>
                        
                        <div class="form-group">
                            <label class="form-label" style="display: flex; justify-content: space-between;">
                                <span>Débit Chlore</span>
                                <span style="font-family: var(--font-mono); color: var(--danger);" id="cl-val-disp">70</span>
                            </label>
                            <input type="range" id="cl-flow" min="10" max="200" value="70" class="form-input" style="padding: 0;">
                            <div class="form-text" id="ratio-display">Ratio Bz/Cl2: 1.4</div>
                        </div>
                        
                        <div class="form-group">
                            <label class="form-label">Configuration CSTR (Série)</label>
                            <div style="display: flex; align-items: center; gap: var(--space-3);">
                                <input type="range" id="n-reactors" min="1" max="5" value="3" class="form-input" style="padding: 0; flex: 1;">
                                <span style="background: var(--bg-tertiary); padding: 0.2rem 0.6rem; border-radius: var(--radius-sm); font-family: var(--font-mono); font-weight: bold;" id="n-reactors-val">3</span>
                            </div>
                        </div>

                        <div class="form-group">
                            <label class="form-label" style="display: flex; justify-content: space-between;">
                                <span>Volume Total (L)</span>
                                <span style="font-family: var(--font-mono);" id="vol-disp">2000</span>
                            </label>
                            <input type="range" id="vol-total" min="500" max="5000" step="100" value="2000" class="form-input" style="padding: 0;">
                        </div>

                        <div class="form-group">
                            <label class="form-label" style="display: flex; justify-content: space-between;">
                                <span>Température (°C)</span>
                                <span style="font-family: var(--font-mono); color: var(--warning);" id="temp-val">55°C</span>
                            </label>
                            <input type="range" id="temp-ctrl" min="20" max="80" value="55" class="form-input" style="padding: 0;">
                        </div>

                        <div class="btn-group" style="margin-top: var(--space-4);">
                            <button id="solve-btn-cb" class="btn btn-primary" style="flex: 1;">
                                ⚡ Actualiser
                            </button>
                            <button id="optimize-btn-cb" class="btn btn-outline" style="flex: 1;">
                                🎯 Optimiser
                            </button>
                        </div>
                    </div>
                </div>

                <!-- Heat Report -->
                <div class="card">
                    <div class="card-header" style="border-left: 4px solid var(--danger);">
                        <h3 style="color: var(--danger); font-size: var(--font-size-base); margin: 0;">Bilan Thermique (Refroidissement)</h3>
                    </div>
                    <div class="card-body">
                        <div id="heat-report" style="font-family: var(--font-mono); font-size: var(--font-size-xs); background: var(--bg-tertiary); padding: var(--space-3); border-radius: var(--radius-md); max-height: 150px; overflow-y: auto;">
                           Calcul en cours...
                        </div>
                    </div>
                </div>
            </div>

            <!-- Right Panel: Visualization -->
            <div style="display: flex; flex-direction: column; gap: var(--space-4);">
                
                <!-- KPI Cards -->
                <div style="display: grid; grid-template-columns: repeat(4, 1fr); gap: var(--space-3);">
                    <div class="card" style="text-align: center; padding: var(--space-3);">
                        <div style="font-size: var(--font-size-xs); text-transform: uppercase; color: var(--text-muted); letter-spacing: 0.05em;">Conversion Bz</div>
                        <div style="font-size: var(--font-size-xl); font-weight: 700; margin-top: var(--space-1);" id="res-conv">--%</div>
                    </div>
                    <div class="card" style="text-align: center; padding: var(--space-3);">
                        <div style="font-size: var(--font-size-xs); text-transform: uppercase; color: var(--text-muted); letter-spacing: 0.05em;">Sélectivité</div>
                        <div style="font-size: var(--font-size-xl); font-weight: 700; color: var(--success); margin-top: var(--space-1);" id="res-sel">--%</div>
                    </div>
                    <div class="card" style="text-align: center; padding: var(--space-3);">
                        <div style="font-size: var(--font-size-xs); text-transform: uppercase; color: var(--text-muted); letter-spacing: 0.05em;">Production MCB</div>
                        <div style="font-size: var(--font-size-xl); font-weight: 700; color: var(--info); margin-top: var(--space-1);" id="res-prod">-- kg/h</div>
                    </div>
                    <div class="card" style="text-align: center; padding: var(--space-3);">
                        <div style="font-size: var(--font-size-xs); text-transform: uppercase; color: var(--text-muted); letter-spacing: 0.05em;">Pureté</div>
                        <div style="font-size: var(--font-size-xl); font-weight: 700; color: #9b59b6; margin-top: var(--space-1);" id="res-purity">--%</div>
                    </div>
                </div>

                <!-- Profile Chart -->
                <div class="card" style="flex: 1; min-height: 0;">
                    <div class="card-header">
                        <h3>Profil de Concentration (Réacteurs en Série)</h3>
                    </div>
                    <div class="card-body" style="height: 350px; position: relative;">
                        <canvas id="profileChart"></canvas>
                    </div>
                </div>
                
                <!-- Separation Scheme Metrics -->
                <div class="card">
                    <div class="card-header">
                        <h3 style="font-size: var(--font-size-base);">Performance Système de Séparation</h3>
                    </div>
                    <div class="card-body">
                         <div style="display: flex; justify-content: space-between; align-items: center; background: var(--bg-tertiary); padding: var(--space-3); border-radius: var(--radius-md); font-size: var(--font-size-sm);">
                            <div style="text-align: center;">
                                <div style="color: var(--text-muted); font-size: var(--font-size-xs);">Sortie Réacteurs</div>
                                <div style="font-family: var(--font-mono); font-weight: 600;" id="flow-reactor-out">-- MCB</div>
                            </div>
                            <div style="height: 1px; background: var(--border-color); flex: 1; margin: 0 var(--space-4); position: relative;">
                                <span style="position: absolute; top: -8px; left: 50%; transform: translateX(-50%); font-size: 10px; color: var(--text-muted); background: var(--bg-tertiary); padding: 0 4px;">Flux</span>
                            </div>
                            <div style="text-align: center;">
                                <div style="color: var(--text-muted); font-size: var(--font-size-xs);">Recycle Benzène</div>
                                <div style="font-family: var(--font-mono); font-weight: 600; color: var(--info);" id="recycle-info">--</div>
                            </div>
                            <div style="height: 1px; background: var(--border-color); flex: 1; margin: 0 var(--space-4);"></div>
                            <div style="text-align: center;">
                                <div style="color: var(--text-muted); font-size: var(--font-size-xs);">Produit Final</div>
                                <div style="font-family: var(--font-mono); font-weight: 600; color: var(--success);" id="flow-product">--</div>
                            </div>
                        </div>
                    </div>
                </div>

            </div>
        </div>
    </div>
    `;
};

// ============================================================================
// SIMULATION LOGIC & BINDINGS
// ============================================================================

ChlorobenzeneModule.init = async function() {
    this.bindEvents();
    this.updateSim();
};

ChlorobenzeneModule.bindEvents = function() {
    const inputs = ['bz-flow', 'cl-flow', 'n-reactors', 'vol-total', 'temp-ctrl'];
    inputs.forEach(id => {
        const el = document.getElementById(id);
        if(el) {
            el.addEventListener('input', (e) => {
                // Update value displays immediately
                if(id === 'bz-flow') document.getElementById('bz-val-disp').innerText = e.target.value;
                if(id === 'cl-flow') document.getElementById('cl-val-disp').innerText = e.target.value;
                if(id === 'vol-total') document.getElementById('vol-disp').innerText = e.target.value;
                if(id === 'n-reactors') document.getElementById('n-reactors-val').innerText = e.target.value;
                if(id === 'temp-ctrl') document.getElementById('temp-val').innerText = e.target.value + '°C';
                
                this.updateSim();
            });
        }
    });
    
    document.getElementById('solve-btn-cb')?.addEventListener('click', () => this.updateSim());
    document.getElementById('optimize-btn-cb')?.addEventListener('click', () => this.runOptimization());
};

ChlorobenzeneModule.getParams = function() {
    return {
        benzeneFlow: parseFloat(document.getElementById('bz-flow').value),
        chlorineFlow: parseFloat(document.getElementById('cl-flow').value),
        numberOfReactors: parseInt(document.getElementById('n-reactors').value),
        reactorVolume: parseFloat(document.getElementById('vol-total').value),
        temperature: parseFloat(document.getElementById('temp-ctrl').value),
        k1_ref: this.defaults.k1_ref,
        k2_ref: this.defaults.k2_ref,
        Ea1: this.defaults.Ea1,
        Ea2: this.defaults.Ea2
    };
};

ChlorobenzeneModule.updateSim = function() {
    const params = this.getParams();
    
    // Update ratio
    const ratio = params.benzeneFlow / params.chlorineFlow;
    const ratioEl = document.getElementById('ratio-display');
    if(ratioEl) ratioEl.innerText = `Ratio Bz/Cl2: ${ratio.toFixed(2)}`;
    if(ratioEl) ratioEl.className = ratio < 1 ? "text-xs text-red-400 mt-1 font-bold" : "text-xs text-gray-500 mt-1";

    // Run Calculation
    const res = calculateChlorobenzeneProcess(params);
    this.results = res;

    // Render Metrics
    document.getElementById('res-conv').innerText = res.metrics.conversionBenzene.toFixed(1) + '%';
    document.getElementById('res-sel').innerText = res.metrics.selectivity.toFixed(1) + '%';
    
    // Prod kg/h
    const prodKgH = res.separation.product.monochlor * 112.56 * 0.06; // mol/min * g/mol * 60min/h * 1kg/1000g
    document.getElementById('res-prod').innerText = prodKgH.toFixed(1) + ' kg/h';
    document.getElementById('res-purity').innerText = res.metrics.purity.toFixed(1) + '%';
    
    // Separation Info
    const lastReactor = res.reactors[res.reactors.length-1];
    document.getElementById('flow-reactor-out').innerText = lastReactor.flows.monochlor.toFixed(1) + ' mol/m';
    document.getElementById('flow-product').innerText = res.separation.product.monochlor.toFixed(1) + ' mol/m';
    document.getElementById('recycle-info').innerText = res.separation.recycle.benzene.toFixed(1) + ' mol/m';

    // Heat Report
    const totalDuty = res.reactors.reduce((a,b)=>a+b.duty, 0);
    const heatHtml = res.reactors.map(r => 
        `<div class="flex justify-between"><span>R${r.id} (${(r.flows.monochlor).toFixed(1)} MCB):</span> <span class="text-red-400">${r.duty.toFixed(1)} kW</span></div>`
    ).join('');
    document.getElementById('heat-report').innerHTML = heatHtml + 
        `<div class="border-t border-gray-600 mt-2 pt-1 flex justify-between font-bold text-white"><span>TOTAL:</span> <span>${totalDuty.toFixed(1)} kW</span></div>`;

    // Charts
    this.renderCharts(res);
};

ChlorobenzeneModule.runOptimization = function() {
    // Strategy: High Bz/Cl2 ratio, Moderate Temp
    document.getElementById('bz-flow').value = 150;
    document.getElementById('cl-flow').value = 50; 
    document.getElementById('temp-ctrl').value = 45;
    document.getElementById('n-reactors').value = 4;
    
    // Trigger update
    document.getElementById('bz-val-disp').innerText = 150;
    document.getElementById('cl-val-disp').innerText = 50;
    document.getElementById('temp-val').innerText = '45°C';
    document.getElementById('n-reactors-val').innerText = 4;
    
    this.updateSim();
};

let chartInstance = null;
ChlorobenzeneModule.renderCharts = function(res) {
    const ctx = document.getElementById('profileChart');
    if(!ctx) return;
    
    // Check if Chart.js is available
    if (typeof Chart === 'undefined') {
        const parent = ctx.parentElement || ctx.parentNode;
        if (parent && !parent.querySelector('.chart-error-msg')) {
             parent.innerHTML = '<div class="chart-error-msg p-4 text-center text-muted">Graphique non disponible (Erreur de chargement librairie)</div>';
        }
        return;
    }

    const labels = [0, ...res.reactors.map(r => `R${r.id}`)];
    const dataBz = [res.inputs.benzeneFlow, ...res.reactors.map(r => r.flows.benzene)];
    const dataCl = [res.inputs.chlorineFlow, ...res.reactors.map(r => r.flows.chlorine)];
    const dataMono = [0, ...res.reactors.map(r => r.flows.monochlor)];
    const dataDi = [0, ...res.reactors.map(r => r.flows.dichlor)];

    if(chartInstance) chartInstance.destroy();
    
    chartInstance = new Chart(ctx, {
        type: 'line',
        data: {
            labels: labels,
            datasets: [
                { label: 'Benzène', data: dataBz, borderColor: '#3b82f6', tension: 0.1 },
                { label: 'Chlore', data: dataCl, borderColor: '#ef4444', tension: 0.1 },
                { label: 'Mono-Cl (Produit)', data: dataMono, borderColor: '#10b981', tension: 0.1, borderWidth: 3 },
                { label: 'Di-Cl', data: dataDi, borderColor: '#8b5cf6', tension: 0.1, borderDash: [5,5] }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            interaction: { mode: 'index', intersect: false },
            scales: {
                y: { 
                    beginAtZero: true,
                    grid: { color: 'rgba(255,255,255,0.05)' },
                    ticks: { color: '#9ca3af' }
                },
                x: {
                    grid: { color: 'rgba(255,255,255,0.05)' },
                    ticks: { color: '#9ca3af' }
                }
            },
            plugins: {
                legend: { position: 'bottom', labels: { color: '#d1d5db', boxWidth: 10 } }
            }
        }
    });
};

ChlorobenzeneModule.getExplanation = function() {
    return {
        title: 'Synthèse du Chlorobenzène',
        description: 'Simulation avancée incluant bilans thermiques et séparation.',
        theory: 'CSTR en série avec cinétique consécutive et recycle.',
        formulas: '$$ r_1 = k_1 C_A C_B $$',
        references: ['Levenspiel', 'Fogler']
    };
};

if (typeof registerModule !== 'undefined') {
    registerModule('chlorobenzene', ChlorobenzeneModule);
}

window.ChlorobenzeneModule = ChlorobenzeneModule;
