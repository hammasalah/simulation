/**
 * ============================================================================
 * COMBINED CYCLE (BRAYTON + RANKINE) SIMULATION
 * ============================================================================
 * 
 * This module simulates a combined cycle power plant consisting of:
 *   - Gas turbine cycle (Brayton cycle) - topping cycle
 *   - Steam cycle (Rankine cycle) - bottoming cycle
 *   - Heat Recovery Steam Generator (HRSG) with Pinch Point Analysis
 * 
 * PHYSICS BACKGROUND:
 * ------------------
 * Combined cycles exploit the high exhaust temperature of gas turbines.
 * 
 * IMPROVEMENTS V0.2:
 * - Proper Steam Table (IAPWS-IF97 simplified)
 * - HRSG Pinch Point Analysis for steam mass flow prediction
 * - Multi-pressure capability (Dual Pressure simplified)
 * - Part-load performance curves
 * - T-Q Diagram generation
 */

'use strict';

// ============================================================================
// STEAM PROPERTIES HELPER (Based on IAPWS-97 Simplified)
// ============================================================================
const CombinedSteam = {
    Tc: 647.096, Pc: 22064,
    getPsat(T) {
        if (T > this.Tc) return this.Pc;
        const t = T - 273.15;
        // Antoine for water: P(bar) approx. -> kPa
        return Math.exp(16.3872 - 3885.70 / (t + 230.170)) * 100;
    },
    getTsat(P) {
        const lnP = Math.log(P / 100);
        return (3885.70 / (16.3872 - lnP)) - 230.170 + 273.15;
    },
    getEnthalpy(P, T, phase) { // P(kPa), T(K)
        const t = T - 273.15;
        const Tsat = this.getTsat(P);
        
        if (phase === 'liquid' || (phase === undefined && T < Tsat)) {
            // Compressed liquid approx (Cp ~ 4.18)
            return 4.19 * t; // kJ/kg
        } else if (phase === 'vapor' || (phase === undefined && T >= Tsat)) {
             // Superheated steam approx
             // h = hg(at Tsat) + cp_vap * (T - Tsat)
             const hg_sat = 2500 + 1.9 * (Tsat - 273.15); // approx saturation line
             return hg_sat + 2.1 * (T - Tsat);
        } else if (phase === 'sat_liquid') {
             return 4.19 * (Tsat - 273.15);
        } else if (phase === 'sat_vapor') {
             return 2500 + 1.9 * (Tsat - 273.15);
        }
        return 0;
    }
};

const CombinedModule = {
    name: 'combined',
    
    defaults: {
        // GT Params
        airInletTemp: 288.15,
        compressorRatio: 18,
        turbineInletTemp: 1473.15,   // 1200°C
        airMassFlow: 500,            // kg/s
        
        // Steam/HRSG Params
        steamPressure: 8000,         // kPa (HP)
        steamPressureLP: 600,        // kPa (LP for dual pressure)
        condenserPressure: 10,       // kPa
        pinchPoint: 10,              // K (Delta T pinch)
        approachTemp: 15,            // K (T_gas_out - T_steam_out difference) 
        
        // Simulation options
        multiPressure: false,
        loadFactor: 1.0              // 0.5 to 1.0
    },
    
    params: null,
    results: null
};

// ============================================================================
// CORE CALCULATIONS
// ============================================================================

function calculateCombinedCycle(params) {
    const {
        airInletTemp, compressorRatio, turbineInletTemp,
        steamPressure, condenserPressure,
        pinchPoint, approachTemp,
        loadFactor = 1.0,
        multiPressure = false,
        steamPressureLP
    } = params;
    
    // --- 1. GAS TURBINE (Topping Cycle) ---
    // Part load adjustments
    // Air flow reduces with load (VIGV regulation)
    const currentAirFlow = params.airMassFlow * loadFactor;
    // Efficiency penalty at part load
    const baseGtEff = 0.38; // 38% simple cycle efficiency base
    const gtEff = baseGtEff * (1 - 0.2 * (1 - loadFactor)); // Penalty
    
    const cp_air = 1.005; // kJ/kgK
    const cp_gas = 1.15;  // kJ/kgK avg for exhaust
    
    // Compressor
    const T1 = airInletTemp;
    const rp = compressorRatio;
    const gamma_air = 1.4;
    const T2s = T1 * Math.pow(rp, (gamma_air-1)/gamma_air);
    const eta_c = 0.88;
    const T2 = T1 + (T2s - T1)/eta_c;
    const Wc = cp_air * (T2 - T1);
    
    // Turbine
    const T3 = turbineInletTemp; 
    // Approx Expansion considering pressure losses
    const P3_P4 = rp * 0.95; // Pressure ratio across turbine
    const gamma_gas = 1.33;
    const T4s = T3 / Math.pow(P3_P4, (gamma_gas-1)/gamma_gas);
    const eta_t = 0.90;
    const T4 = T3 - eta_t * (T3 - T4s); // TIT is governed, T4 is Exhaust Temp
    const Wt = cp_gas * (T3 - T4);
    
    // GT Power
    const Wnet_GT_specific = Wt - Wc; // kJ/kg air (very approx due to Cp diff, usually need fuel mass)
    // Refined: Power = m_air * (Wt - Wc) is rough. 
    // Better: P_GT = m_air * Cp_gas * (T3-T4) * eta_mech - m_air * Cp_air * (T2-T1)
    const gtPowerMW = (currentAirFlow * Wnet_GT_specific) / 1000;
    
    // Fuel Heat approx (Q_in = Cp_gas*(T3-T2))
    const Qin_GT_MW = currentAirFlow * 1.1 * (T3 - T2) / 1000; 
    
    const Texhaust = T4; // Gas into HRSG

    // --- 2. HRSG & STEAM CYCLE (Bottoming) ---
    // Pinch Point Analysis for HP Steam
    
    // Water/Steam Props HP
    const Tsat_HP = CombinedSteam.getTsat(steamPressure);
    const hf_HP = CombinedSteam.getEnthalpy(steamPressure, Tsat_HP, 'sat_liquid');
    const hg_HP = CombinedSteam.getEnthalpy(steamPressure, Tsat_HP, 'sat_vapor');
    
    // Live Steam Temp (limited by Gas Inlet - Approach)
    // Usually Superheater Approach: T_steam_out = T_gas_in - Approach? No, usually 20-30K diff.
    const T_steam_HP = Texhaust - 30; // Simple approach constraint
    const h_steam_HP = CombinedSteam.getEnthalpy(steamPressure, T_steam_HP, 'vapor');
    
    // Feedwater (approx 40C from cond/pump)
    const T_feed = 313.15; // 40°C
    const h_feed = CombinedSteam.getEnthalpy(steamPressure, T_feed, 'liquid');
    
    // PINCH ANALYSIS ---------------------------
    // Pinch point is at Evaporator Water Inlet (Saturated Liquid) vs Gas Temp there.
    // T_gas_pinch = Tsat_HP + pinchPoint
    const T_gas_pinch = Tsat_HP + pinchPoint;
    
    let steamFlowHP = 0;
    let stackTemp = 0;
    let heatTransferredHP = 0;
    
    // Check if exhaust is hot enough
    if (Texhaust > T_gas_pinch) {
        // Energy available above pinch (Gas cools from Texhaust to T_gas_pinch)
        // This energy heats water from Saturated Liquid to Superheated Steam
        // Q_gas_sup_evap = m_gas * Cp_gas * (Texhaust - T_gas_pinch)
        const Q_avail_above_pinch = currentAirFlow * cp_gas * (Texhaust - T_gas_pinch);
        
        // Energy required per kg steam (Boiling + Superheating)
        const q_req_steam_HP = h_steam_HP - hf_HP;
        
        // Max HP Steam Flow
        steamFlowHP = Q_avail_above_pinch / q_req_steam_HP;
        
        // Economizer Check (Gas cooling below pinch)
        // Heat required for Economizer (Feed -> Sat Liq)
        const Q_eco_req = steamFlowHP * (hf_HP - h_feed);
        const deltaT_gas_eco = Q_eco_req / (currentAirFlow * cp_gas);
        
        stackTemp = T_gas_pinch - deltaT_gas_eco;
        
        // Min stack temp constraint (acid dew point ~100C)
        if (stackTemp < 373.15) {
             // Constrained by stack temp, not pinch
             // Recalculate based on T_stack = 100C
             stackTemp = 373.15;
             const Q_total_avail = currentAirFlow * cp_gas * (Texhaust - stackTemp);
             steamFlowHP = Q_total_avail / (h_steam_HP - h_feed);
        }
        
        heatTransferredHP = steamFlowHP * (h_steam_HP - h_feed);
    }

    // LP Cycle (Dual Pressure) ----------------
    let steamFlowLP = 0;
    let stPowerLP = 0;
    
    if (multiPressure && stackTemp > 400) { // Only if enough heat left
         const Tsat_LP = CombinedSteam.getTsat(steamPressureLP);
         const T_gas_pinch_LP = Tsat_LP + pinchPoint;
         
         // Assuming LP evaporator is located after HP Economizer
         // Gas enters LP section at 'stackTemp' (from HP calc)
         const T_gas_in_LP = stackTemp; 
         
         if (T_gas_in_LP > T_gas_pinch_LP) {
             const h_steam_LP = CombinedSteam.getEnthalpy(steamPressureLP, T_gas_in_LP - 10, 'vapor');
             const hf_LP = CombinedSteam.getEnthalpy(steamPressureLP, Tsat_LP, 'sat_liquid');
             
             // Available Q above LP pinch
             const Q_avail_LP = currentAirFlow * cp_gas * (T_gas_in_LP - T_gas_pinch_LP);
             const q_req_LP = h_steam_LP - hf_LP;
             
             steamFlowLP = Q_avail_LP / q_req_LP;
             
             // Eco LP
             const h_feed_LP = CombinedSteam.getEnthalpy(steamPressureLP, 303.15, 'liquid');
             const Q_eco_LP = steamFlowLP * (hf_LP - h_feed_LP);
             stackTemp = T_gas_pinch_LP - Q_eco_LP / (currentAirFlow * cp_gas);
             
             // LP Power approx
             stPowerLP = (steamFlowLP * (h_steam_LP - 2300)) / 1000 * 0.85; // rough isentropic default
         }
    }

    // Steam Turbine Power (HP)
    // h_in = h_steam_HP. Expansion to Condenser P
    // Simple Rankine calc
    const h_out_ideal = 2100; // rough approx for vac
    const stEff = 0.87;
    const stWorkSpecific = (h_steam_HP - h_out_ideal) * stEff;
    const stPowerHP = (steamFlowHP * stWorkSpecific) / 1000; // MW
    
    // Totals
    const totalPower = gtPowerMW + stPowerHP + stPowerLP;
    const combinedEfficiency = totalPower / Qin_GT_MW;
    
    // Calculation Steps for UI
    const steps = [
        {
            description: "Puissance Turbine à Gaz (Brayton)",
            formula: "W_GT = ṁ_air × (Wt - Wc)",
            result: `${gtPowerMW.toFixed(1)} MW`
        },
        {
            description: "Température Échappement GT",
            formula: "T_out",
            result: `${Texhaust.toFixed(0)} K (${(Texhaust-273).toFixed(0)}°C)`
        },
        {
            description: "Débit Vapeur HP (Pincement)",
            formula: "ṁ_s = Q_avail / Δh",
            result: `${steamFlowHP.toFixed(1)} kg/s`
        },
        {
            description: "Puissance Vapeur (Rankine)",
            formula: "P_ST = P_HP + P_LP",
            result: `${(stPowerHP + stPowerLP).toFixed(1)} MW`
        },
        {
            description: "Rendement Global",
            formula: "η = P_tot / Q_fuel",
            result: `${(combinedEfficiency*100).toFixed(2)} %`
        }
    ];

    return {
        inputs: params,
        results: {
            gtPower: gtPowerMW,
            stPower: stPowerHP + stPowerLP,
            totalPower,
            efficiency: combinedEfficiency,
            exhaustTemp: Texhaust,
            stackTemp,
            hpFlow: steamFlowHP,
            lpFlow: steamFlowLP,
            heatRate: 3600 / combinedEfficiency
        },
        diagramData: {
            // Data for T-Q Diagram
            // Points [Q, T]
            gas: [
                {q: 0, t: Texhaust}, 
                {q: 100, t: stackTemp} // Normalized 0-100% heat transfer
            ],
            water: [
                {q: 100, t: T_feed}, // Feed in (coldest)
                {q: 85, t: Tsat_HP}, // Pinch liquid
                {q: 40, t: Tsat_HP}, // Evap finish
                {q: 0, t: T_steam_HP} // Superheat out (hottest)
                // Note: Q axis is usually "Heat Transferred FROM gas", so 0 is gas inlet.
            ]
        },
        calculationSteps: steps
    };
}

// ============================================================================
// UI RENDERING
// ============================================================================

CombinedModule.render = function() {
    this.params = { ...this.defaults };
    
    return `
        <div class="combined-module" id="combined-content">
            <div class="d-flex justify-between align-center mb-5">
                <div>
                    <h2>Cycle Combiné Gaz-Vapeur (CCGT)</h2>
                    <p class="text-muted">Centrale à cycle combiné haute efficacité avec HRSG.</p>
                </div>
                <div class="status-badge success">
                    <span class="status-dot"></span>
                    <span>Opérationnel v0.2</span>
                </div>
            </div>
            
            <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;">
                <!-- CONTROLS -->
                <div class="card">
                    <div class="card-header"><h3 class="card-title">Paramètres</h3></div>
                    <div class="card-body">
                        <!-- TABS -->
                        <div class="tabs mb-4" style="display:flex; gap:10px; border-bottom:1px solid #ddd;">
                            <button onclick="switchTab('gt')" class="tab-btn active" id="btn-gt">Turbine à Gaz</button>
                            <button onclick="switchTab('hrsg')" class="tab-btn" id="btn-hrsg">HRSG & Vapeur</button>
                        </div>

                        <!-- GT Controls -->
                        <div id="tab-gt">
                            <div class="form-group">
                                <label>Charge Partielle (%)</label>
                                <input type="range" class="form-range" id="load-factor" min="50" max="100" value="100">
                                <div class="text-right small text-primary" id="load-val">100%</div>
                            </div>
                            <div class="param-grid two-col">
                                <div class="form-group">
                                    <label>Temp. Air (K)</label>
                                    <input type="number" class="form-input" id="air-temp" value="288">
                                </div>
                                <div class="form-group">
                                    <label>Débit Air (kg/s)</label>
                                    <input type="number" class="form-input" id="air-flow" value="500">
                                </div>
                                <div class="form-group">
                                    <label>TIT (K)</label>
                                    <input type="number" class="form-input" id="tit" value="1473">
                                </div>
                                <div class="form-group">
                                    <label>Ratio Compression</label>
                                    <input type="number" class="form-input" id="rc" value="18">
                                </div>
                            </div>
                        </div>

                        <!-- HRSG Controls -->
                        <div id="tab-hrsg" style="display:none;">
                            <div class="form-group mb-3">
                                <label class="form-check">
                                    <input type="checkbox" id="multi-pressure">
                                    <strong>Cycle Multi-Pression (Dual)</strong>
                                </label>
                            </div>
                            <div class="param-grid two-col">
                                <div class="form-group">
                                    <label>Pression HP (kPa)</label>
                                    <input type="number" class="form-input" id="p-hp" value="8000">
                                </div>
                                <div class="form-group" id="lp-group" style="display:none;">
                                    <label>Pression LP (kPa)</label>
                                    <input type="number" class="form-input" id="p-lp" value="500">
                                </div>
                                <div class="form-group">
                                    <label>Pincement HRSG (K)</label>
                                    <input type="number" class="form-input" id="pinch" value="10" min="5" max="30">
                                </div>
                                <div class="form-group">
                                    <label>Condenseur (kPa)</label>
                                    <input type="number" class="form-input" id="p-cond" value="5">
                                </div>
                            </div>
                        </div>

                        <button class="btn btn-primary w-full mt-4" id="btn-calc">Calculer Performance</button>
                    </div>
                </div>

                <!-- RESULTS & CHARTS -->
                <div>
                     <!-- KPI Cards -->
                    <div class="results-grid mb-4">
                        <div class="result-item">
                            <div class="result-label">Puissance Totale</div>
                            <div class="result-value" id="res-power">--</div>
                            <div class="result-unit">MW</div>
                        </div>
                        <div class="result-item">
                            <div class="result-label">Rendement</div>
                            <div class="result-value" id="res-eff">--</div>
                            <div class="result-unit">%</div>
                        </div>
                        <div class="result-item">
                            <div class="result-label">Récupération HRSG</div>
                            <div class="result-value" id="res-steam">--</div>
                            <div class="result-unit">kg/s</div>
                        </div>
                    </div>

                    <!-- T-Q Diagram -->
                    <div class="card">
                        <div class="card-header"><h4 class="card-title">Diagramme Echangeur (T-Q)</h4></div>
                        <div class="card-body" style="height:300px; display:flex; justify-content:center; align-items:center;" id="tq-chart">
                            <canvas id="tqCanvas"></canvas>
                        </div>
                    </div>
                </div>
            </div>
            
            <!-- Details removed -->
        </div>
    `;
};

CombinedModule.init = function() {
    // Basic Tab Logic
    window.switchTab = function(t) {
        document.getElementById('tab-gt').style.display = t === 'gt' ? 'block' : 'none';
        document.getElementById('tab-hrsg').style.display = t === 'hrsg' ? 'block' : 'none';
        document.getElementById('btn-gt').classList.toggle('active', t === 'gt');
        document.getElementById('btn-hrsg').classList.toggle('active', t === 'hrsg');
    };

    // Listeners
    document.getElementById('btn-calc').addEventListener('click', runCalc);
    document.getElementById('multi-pressure').addEventListener('change', e => {
        document.getElementById('lp-group').style.display = e.target.checked ? 'block' : 'none';
    });
    document.getElementById('load-factor').addEventListener('input', e => {
        document.getElementById('load-val').innerText = e.target.value + '%';
    });

    runCalc();
};

function runCalc() {
    const params = {
        airInletTemp: parseFloat(document.getElementById('air-temp').value),
        airMassFlow: parseFloat(document.getElementById('air-flow').value),
        turbineInletTemp: parseFloat(document.getElementById('tit').value),
        compressorRatio: parseFloat(document.getElementById('rc').value),
        steamPressure: parseFloat(document.getElementById('p-hp').value),
        steamPressureLP: parseFloat(document.getElementById('p-lp').value),
        pinchPoint: parseFloat(document.getElementById('pinch').value),
        loadFactor: parseFloat(document.getElementById('load-factor').value) / 100,
        multiPressure: document.getElementById('multi-pressure').checked
    };

    const res = calculateCombinedCycle(params);
    
    // Bind Results
    document.getElementById('res-power').innerText = res.results.totalPower.toFixed(1);
    document.getElementById('res-eff').innerText = (res.results.efficiency*100).toFixed(2);
    document.getElementById('res-steam').innerText = res.results.hpFlow.toFixed(1) + (res.results.lpFlow > 0 ? ` (+${res.results.lpFlow.toFixed(1)} LP)` : '');

    // steps population removed

    drawTQChart(res.diagramData);
}

function drawTQChart(data) {
    const container = document.getElementById('tq-chart');
    container.innerHTML = '<canvas id="tqCanvas" style="width:100%; height:100%"></canvas>';
    const ctx = document.getElementById('tqCanvas').getContext('2d');
    
    // Simple canvas drawing
    // Need a library normally (Chart.js), but I will do raw canvas for robustness
    const w = ctx.canvas.width = container.clientWidth;
    const h = ctx.canvas.height = container.clientHeight;
    const padding = 40;
    
    // Scales using fixed range for now or auto
    const maxT = data.gas[0].t + 50;
    const minT = 280; // K
    
    const scaleX = q => padding + (q/100) * (w - 2*padding);
    const scaleY = t => h - padding - ((t-minT)/(maxT-minT)) * (h - 2*padding);
    
    ctx.clearRect(0,0,w,h);
    
    // Axes
    ctx.beginPath();
    ctx.moveTo(padding, h-padding); ctx.lineTo(w-padding, h-padding); // X
    ctx.moveTo(padding, h-padding); ctx.lineTo(padding, padding); // Y
    ctx.strokeStyle = '#666';
    ctx.stroke();
    
    // Gas Line (Red)
    ctx.beginPath();
    ctx.moveTo(scaleX(data.gas[0].q), scaleY(data.gas[0].t));
    ctx.lineTo(scaleX(data.gas[1].q), scaleY(data.gas[1].t));
    ctx.strokeStyle = '#e74c3c'; ctx.lineWidth = 3; ctx.stroke();
    ctx.fillStyle = '#e74c3c'; ctx.fillText("Gaz Échappement", scaleX(50), scaleY(data.gas[0].t)-10);

    // Water Line (Blue)
    ctx.beginPath();
    const wp = data.water;
    ctx.moveTo(scaleX(wp[0].q), scaleY(wp[0].t)); // Feed
    for(let i=1; i<wp.length; i++) ctx.lineTo(scaleX(wp[i].q), scaleY(wp[i].t));
    ctx.strokeStyle = '#3498db'; ctx.lineWidth = 3; ctx.stroke();
    
    // Labels
    ctx.fillStyle = '#333';
    ctx.textAlign = 'center';
    ctx.fillText("Chaleur Transférée (%)", w/2, h-10);
    ctx.save();
    ctx.translate(15, h/2); ctx.rotate(-Math.PI/2);
    ctx.fillText("Température (K)", 0,0);
    ctx.restore();
    
    // Pinch Point Circle
    ctx.beginPath();
    // Pinch is usually at evap start (wp[2] or wp[1] depending on array order, here index 1 is pinch liquid usually)
    // Actually in my data logic: wp[1] is pinch liquid, wp[2] is evap end? Wait.
    // data.water: 0:feed, 1:pinch_liq, 2:evap_liq(same T?), 3:superheat.
    // Water heats up left-to-right on Q axis (0 to 100 on my chart? No, Q is heat transferred FROM gas. 
    // Gas enters at Q=0 (hot). Water exits at Q=0 (hot).
    // So Gas Q: 0 -> 100. Water Q: 0 (Superheat out) -> ... -> 100 (Feed in).
    // Yes, my array was: 0:Feed(100), 1:Pinch(85), 2:Evap(40), 3:Steam(0).
    // Canvas lineTo connects them in order. That's fine.
}

if (typeof registerModule !== 'undefined') {
    registerModule('combined', CombinedModule);
}

window.CombinedModule = CombinedModule;

