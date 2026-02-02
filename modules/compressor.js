/**
 * ============================================================================
 * COMPRESSOR SIMULATION MODULE - FULLY IMPLEMENTED
 * ============================================================================
 * 
 * This module provides a complete simulation of a single-stage centrifugal or
 * piston compressor with both ideal (isentropic) and real (efficiency-corrected)
 * calculations.
 * 
 * PHYSICS BACKGROUND:
 * ------------------
 * A compressor is a device that increases the pressure of a gas by doing work on it.
 * The work required depends on the thermodynamic path taken during compression.
 * 
 * For an ideal gas undergoing isentropic (reversible adiabatic) compression:
 *   T₂ₛ/T₁ = (P₂/P₁)^((γ-1)/γ)
 * 
 * The isentropic efficiency relates actual to ideal work:
 *   ηₛ = Wₛ/Wₐctual = (T₂ₛ - T₁)/(T₂ - T₁)
 * 
 * Power required:
 *   Ẇ = ṁ × cₚ × (T₂ - T₁)
 * 
 * REFERENCES:
 * -----------
 * - Çengel, Y.A. & Boles, M.A. "Thermodynamics: An Engineering Approach"
 * - Perry's Chemical Engineers' Handbook, 8th Edition
 * - Smith, Van Ness, Abbott "Introduction to Chemical Engineering Thermodynamics"
 * 
 * @author Chemical Engineering Lab Simulation Platform
 * @version 1.0.0
 */

'use strict';

// ============================================================================
// MODULE CONFIGURATION AND CONSTANTS
// ============================================================================

const CompressorModule = {
    name: 'compressor',
    
    // Default parameters for simulation
    defaults: {
        inletTemperature: 298.15,    // K (25°C)
        inletPressure: 101325,       // Pa (1 atm)
        pressureRatio: 3.0,          // dimensionless
        massFlow: 1.0,               // kg/s
        isentropicEfficiency: 0.85,  // 85%
        gasType: 'air',
        compressorType: 'centrifugal'
    },
    
    // Current parameter values
    params: null,
    
    // Simulation results
    results: null,
    
    // Transient simulation data
    transientData: null,
    
    // Charts reference
    charts: {},
    
    // UI state
    showCalculations: false
};

// ============================================================================
// CORE CALCULATION FUNCTIONS
// ============================================================================

/**
 * Perform steady-state compressor calculations
 * 
 * @param {Object} params - Simulation parameters
 * @returns {Object} Calculation results
 */
function calculateCompressor(params) {
    console.log("calculateCompressor called with:", params);
    try {
        const {
            inletTemperature,    // T₁ in K
            inletPressure,       // P₁ in Pa
            pressureRatio,       // r = P₂/P₁
            massFlow,            // ṁ in kg/s
            isentropicEfficiency, // ηₛ (0-1)
            gasType
        } = params;
        
        // Debug check parameters
        if (isNaN(inletTemperature) || isNaN(inletPressure) || isNaN(pressureRatio)) {
            console.error("Invalid input parameters (NaN detected)");
            throw new Error("Paramètres invalides (valeurs non numériques)");
        }

        // Validation
        if (!gasType) throw new Error("Type de gaz non défini");
        
        // Get gas properties
        // Ensure MathUtils is available
        if (typeof MathUtils === 'undefined') {
            console.error("MathUtils is undefined in global scope");
            throw new Error("Bibliothèque MathUtils non chargée");
        }
        
        console.log("MathUtils available. GAS_PROPERTIES keys:", Object.keys(MathUtils.GAS_PROPERTIES || {}));

        const gas = MathUtils.GAS_PROPERTIES[gasType] || MathUtils.GAS_PROPERTIES.air;
        if (!gas) throw new Error("Propriétés du gaz introuvables: " + gasType);
        
        console.log("Gas properties found:", gas);

        const { R, cp, cv, gamma } = gas;
        
        // Step 1: Calculate outlet pressure
        const outletPressure = inletPressure * pressureRatio;  // P₂ = P₁ × r
        
        // Step 2: Calculate isentropic outlet temperature
        // T₂ₛ = T₁ × (P₂/P₁)^((γ-1)/γ)
        const isentropicExponent = (gamma - 1) / gamma;
        const T2_isentropic = inletTemperature * Math.pow(pressureRatio, isentropicExponent);
        
        // Step 3: Calculate actual outlet temperature using isentropic efficiency
        // ηₛ = (T₂ₛ - T₁)/(T₂ - T₁)  →  T₂ = T₁ + (T₂ₛ - T₁)/ηₛ
        const T2_actual = inletTemperature + (T2_isentropic - inletTemperature) / isentropicEfficiency;
        
        // Step 4: Calculate work/power
        // Isentropic work: Wₛ = ṁ × cₚ × (T₂ₛ - T₁)
        const isentropicWork = massFlow * cp * (T2_isentropic - inletTemperature);
        
        // Actual work: Wₐ = ṁ × cₚ × (T₂ - T₁)
        const actualWork = massFlow * cp * (T2_actual - inletTemperature);
        
        console.log("Interim calculations:", { outletPressure, T2_isentropic, T2_actual, actualWork });

        // Step 5: Calculate polytropic exponent
        // From T₂/T₁ = (P₂/P₁)^((n-1)/n)
        let polytropicN = gamma; // Default
        try {
            polytropicN = MathUtils.polytropicExponent(
                inletTemperature, T2_actual, inletPressure, outletPressure
            );
        } catch(e) { console.warn("Error calculating polytropic n", e); }
        
        // Step 6: Calculate densities and specific volumes
        const rho1 = MathUtils.idealGasDensity(inletPressure, R, inletTemperature);
        const rho2 = MathUtils.idealGasDensity(outletPressure, R, T2_actual);
        const v1 = 1 / rho1;  // Specific volume inlet
        const v2 = 1 / rho2;  // Specific volume outlet
        
        // Step 7: Calculate entropy change
        const entropyChange = MathUtils.entropyChange(
            cp, R, inletTemperature, T2_actual, inletPressure, outletPressure
        );
        
        // Step 8: Calculate volumetric flow rates
        const volumeFlowIn = massFlow / rho1;   // m³/s at inlet
        const volumeFlowOut = massFlow / rho2;  // m³/s at outlet
        
        // Store calculation steps for educational display
        const calculationSteps = [
            {
                description: 'Calcul de la pression de sortie',
                formula: `P₂ = P₁ × r = ${(inletPressure/1000).toFixed(2)} kPa × ${pressureRatio.toFixed(2)}`,
                result: `P₂ = ${(outletPressure/1000).toFixed(2)} kPa`
            },
            {
                description: 'Température de sortie isentropique',
                formula: `T₂ₛ = T₁ × (P₂/P₁)^((γ-1)/γ) = ${inletTemperature.toFixed(1)} K × ${pressureRatio.toFixed(2)}^((${gamma.toFixed(3)}-1)/${gamma.toFixed(3)})`,
                result: `T₂ₛ = ${T2_isentropic.toFixed(2)} K (${(T2_isentropic - 273.15).toFixed(1)}°C)`
            },
            {
                description: 'Température de sortie réelle (avec rendement)',
                formula: `T₂ = T₁ + (T₂ₛ - T₁)/ηₛ = ${inletTemperature.toFixed(1)} + (${T2_isentropic.toFixed(1)} - ${inletTemperature.toFixed(1)})/${isentropicEfficiency.toFixed(2)}`,
                result: `T₂ = ${T2_actual.toFixed(2)} K (${(T2_actual - 273.15).toFixed(1)}°C)`
            },
            {
                description: 'Puissance isentropique',
                formula: `Ẇₛ = ṁ × cₚ × (T₂ₛ - T₁) = ${massFlow.toFixed(2)} × ${cp.toFixed(0)} × (${T2_isentropic.toFixed(1)} - ${inletTemperature.toFixed(1)})`,
                result: `Ẇₛ = ${(isentropicWork/1000).toFixed(2)} kW`
            },
            {
                description: 'Puissance réelle requise',
                formula: `Ẇ = ṁ × cₚ × (T₂ - T₁) = ${massFlow.toFixed(2)} × ${cp.toFixed(0)} × (${T2_actual.toFixed(1)} - ${inletTemperature.toFixed(1)})`,
                result: `Ẇ = ${(actualWork/1000).toFixed(2)} kW`
            },
            {
                description: 'Coefficient polytropique',
                formula: `n = 1 / (1 - ln(T₂/T₁)/ln(P₂/P₁))`,
                result: `n = ${polytropicN.toFixed(4)}`
            }
        ];
        
        const finalResults = {
            // Input echo
            inputs: { ...params, gas },
            
            // Primary outputs
            outletTemperature: T2_actual,           // K
            outletTemperatureIsentropic: T2_isentropic, // K
            outletPressure,                         // Pa
            actualPower: actualWork,                // W
            isentropicPower: isentropicWork,        // W
            
            // Additional outputs
            polytropicExponent: polytropicN,
            entropyChange,                          // J/(kg·K)
            inletDensity: rho1,                     // kg/m³
            outletDensity: rho2,                    // kg/m³
            inletSpecificVolume: v1,                // m³/kg
            outletSpecificVolume: v2,               // m³/kg
            volumeFlowIn,                           // m³/s
            volumeFlowOut,                          // m³/s
            
            // Temperature rise
            temperatureRise: T2_actual - inletTemperature,
            temperatureRiseIsentropic: T2_isentropic - inletTemperature,
            
            // Calculation steps for educational display
            calculationSteps,
            
            // Timestamp
            timestamp: new Date().toISOString()
        };
        console.log("Calculations successful, returning:", finalResults);
        return finalResults;

    } catch (error) {
        console.error("Critical error in calculateCompressor:", error);
         // Return a safe 'error' object instead of throwing, so it doesn't crash the UI entirely
        return {
            error: true,
            message: error.message,
            inputs: { ...params, gas: { R: 287, cp: 1005, gamma: 1.4 } }, // Fallback inputs
            outletTemperature: 0,
            outletTemperatureIsentropic: 0,
            outletPressure: 0,
            actualPower: 0,
            isentropicPower: 0,
            polytropicExponent: 1.4,
            entropyChange: 0,
            calculationSteps: []
        };
    }
}

/**
 * Generate T-s (Temperature-Entropy) diagram data points
 * 
 * @param {Object} results - Calculation results
 * @returns {Object} {isentropic: [], actual: []} curve data
 */
function generateTsDiagram(results) {
    const { inputs, outletTemperature, outletTemperatureIsentropic } = results;
    const { inletTemperature, inletPressure, pressureRatio } = inputs;
    const gas = inputs.gas;
    
    const points = 50;
    const isentropicCurve = [];
    const actualCurve = [];
    
    // Reference entropy at inlet (arbitrary zero)
    const s1 = 0;
    
    // Isentropic process: Δs = 0 (vertical line in T-s)
    // For display, we show a vertical line from T1 to T2s
    isentropicCurve.push({ s: s1, T: inletTemperature });
    isentropicCurve.push({ s: s1, T: outletTemperatureIsentropic });
    
    // Actual process: calculate entropy at each point
    // Δs = cp*ln(T2/T1) - R*ln(P2/P1)
    const P2 = inletPressure * pressureRatio;
    const s2 = gas.cp * Math.log(outletTemperature / inletTemperature) 
             - gas.R * Math.log(P2 / inletPressure);
    
    // Generate curved path for actual process
    for (let i = 0; i <= points; i++) {
        const fraction = i / points;
        const T = inletTemperature + fraction * (outletTemperature - inletTemperature);
        // Approximate intermediate pressure (linear in log space)
        const P = inletPressure * Math.pow(pressureRatio, fraction);
        const s = gas.cp * Math.log(T / inletTemperature) - gas.R * Math.log(P / inletPressure);
        actualCurve.push({ s, T });
    }
    
    return { isentropicCurve, actualCurve };
}

/**
 * Generate P-v (Pressure-Volume) diagram data points
 * 
 * @param {Object} results - Calculation results  
 * @returns {Object} Curve data for P-v diagram
 */
function generatePvDiagram(results) {
    const { inputs, outletPressure, outletSpecificVolume, polytropicExponent } = results;
    const { inletPressure } = inputs;
    const v1 = results.inletSpecificVolume;
    const v2 = outletSpecificVolume;
    
    const points = 50;
    const isentropicCurve = [];
    const polytropicCurve = [];
    
    const gamma = inputs.gas.gamma;
    const n = polytropicExponent;
    
    // P*v^γ = constant for isentropic
    const constIsentropic = inletPressure * Math.pow(v1, gamma);
    
    // P*v^n = constant for polytropic
    const constPolytropic = inletPressure * Math.pow(v1, n);
    
    for (let i = 0; i <= points; i++) {
        const fraction = i / points;
        const v = v1 - fraction * (v1 - v2);
        
        // Isentropic: P = const / v^γ
        const P_isen = constIsentropic / Math.pow(v, gamma);
        isentropicCurve.push({ v, P: P_isen });
        
        // Polytropic: P = const / v^n
        const P_poly = constPolytropic / Math.pow(v, n);
        polytropicCurve.push({ v, P: P_poly });
    }
    
    return { isentropicCurve, polytropicCurve };
}

/**
 * Generate power vs pressure ratio curve
 * 
 * @param {Object} baseParams - Base parameters
 * @returns {Array} Array of {pressureRatio, power} points
 */
function generatePowerCurve(baseParams) {
    const curve = [];
    
    for (let r = 1.0; r <= 10.0; r += 0.25) {
        const params = { ...baseParams, pressureRatio: r };
        const result = calculateCompressor(params);
        curve.push({
            pressureRatio: r,
            power: result.actualPower / 1000  // kW
        });
    }
    
    return curve;
}

// ============================================================================
// TRANSIENT SIMULATION
// ============================================================================

/**
 * Run transient startup simulation using RK4 integration
 * 
 * Models mass flow ramp-up from 0 to nominal value with first-order dynamics.
 * 
 * @param {Object} params - Simulation parameters
 * @param {number} rampTime - Time to reach nominal flow (s)
 * @param {number} totalTime - Total simulation time (s)
 * @param {number} dt - Time step (s)
 * @returns {Object} Transient simulation results
 */
function runTransientSimulation(params, rampTime = 5, totalTime = 20, dt = 0.05) {
    const { massFlow } = params;
    
    // State: [massFlow, outletTemperature, power]
    // Using first-order lag for mass flow ramp
    const tau = rampTime / 3;  // Time constant (63.2% of final value at τ)
    
    // Derivative function for mass flow ramp
    const derivative = (t, m) => {
        const targetFlow = t < rampTime 
            ? massFlow * (1 - Math.exp(-t / tau))
            : massFlow;
        return (targetFlow - m) / (tau / 2);
    };
    
    // Run RK4 integration
    const { times, values } = MathUtils.rk4(derivative, 0.01, 0, totalTime, dt);
    
    // Calculate compressor outputs at each time step
    const transientResults = times.map((t, i) => {
        const currentFlow = Math.max(0.01, values[i]);  // Minimum flow to avoid division by zero
        const stepParams = { ...params, massFlow: currentFlow };
        const result = calculateCompressor(stepParams);
        
        return {
            time: t,
            massFlow: currentFlow,
            outletTemperature: result.outletTemperature,
            outletPressure: result.outletPressure,
            power: result.actualPower / 1000  // kW
        };
    });
    
    return {
        times,
        data: transientResults,
        params,
        rampTime,
        totalTime
    };
}

// ============================================================================
// UI RENDERING
// ============================================================================

/**
 * Render the compressor module UI
 * @returns {string} HTML content
 */
CompressorModule.render = function() {
    // Initialize parameters
    this.params = { ...this.defaults };
    
    return `
        <div class="simulation-module" id="compressor-content">
            <!-- Header -->
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 2rem;">
                <div>
                    <h2 class="text-2xl font-bold text-gray-800" data-i18n="compressor.title">Simulation d'un compresseur</h2>
                    <p class="text-gray-600 mt-2" data-i18n="compressor.description">
                        Modélisation d'un compresseur centrifuge ou à piston avec calcul des conditions 
                        de sortie et de la puissance requise.
                    </p>
                </div>
                <div class="status-badge idle" id="simulation-status">
                    <span class="status-dot"></span>
                    <span>Prêt</span>
                </div>
            </div>
            
            <!-- Main Content Layout -->
            <div style="display: flex; gap: 20px; align-items: flex-start; flex-wrap: wrap;">
                
                <!-- Left Column: Controls & Diagram -->
                <div style="flex: 1; min-width: 300px; display: flex; flex-direction: column; gap: 20px;">
                    
                    <!-- Diagram Card -->
                    <div class="card">
                        <div class="card-header" style="display: flex; justify-content: space-between; align-items: center;">
                            <h3 class="card-title" style="display: flex; align-items: center; gap: 8px;">
                                <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <rect x="3" y="3" width="18" height="18" rx="2"/>
                                    <circle cx="12" cy="12" r="4"/>
                                </svg>
                                Schéma du compresseur
                            </h3>
                            <!-- Visualize button removed -->
                        </div>
                        <div class="card-body p-0">
                            <div class="diagram-container relative w-full" id="compressor-diagram" style="min-height: 250px; display: flex; align-items: center; justify-content: center; background: #f8f9fa;">
                                <!-- SVG will be loaded here -->
                                <span class="loading-spinner"></span>
                            </div>
                        </div>
                    </div>
                
                    <!-- Parameters Card -->
                    <div class="card">
                        <div class="card-header">
                            <h3 class="card-title" style="display: flex; align-items: center; gap: 8px;">
                                <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <path d="M12 3v18M3 12h18"/>
                                </svg>
                                Paramètres d'entrée
                            </h3>
                        </div>
                        <div class="card-body">
                            <div style="display: grid; grid-template-columns: repeat(2, 1fr); gap: 15px;">
                                <!-- Gas Type -->
                                <div class="form-group">
                                    <label class="form-label" for="gas-type" data-i18n="compressor.gasType">Type de gaz</label>
                                    <select class="form-select" id="gas-type">
                                        <option value="air" selected>Air</option>
                                        <option value="nitrogen">Azote (N₂)</option>
                                        <option value="oxygen">Oxygène (O₂)</option>
                                        <option value="carbon_dioxide">CO₂</option>
                                        <option value="methane">Méthane (CH₄)</option>
                                        <option value="steam">Vapeur d'eau</option>
                                        <option value="hydrogen">Hydrogène (H₂)</option>
                                        <option value="helium">Hélium (He)</option>
                                        <option value="ideal">Gaz parfait</option>
                                    </select>
                                </div>
                                
                                <!-- Compressor Type -->
                                <div class="form-group">
                                    <label class="form-label" for="compressor-type">Type de compresseur</label>
                                    <select class="form-select" id="compressor-type">
                                        <option value="centrifugal" selected>Centrifuge</option>
                                        <option value="piston">À piston</option>
                                    </select>
                                </div>
                                
                                <!-- Inlet Temperature -->
                                <div class="form-group">
                                    <label class="form-label" for="inlet-temp" data-i18n="compressor.inletTemp">Température d'entrée</label>
                                    <div class="input-with-unit">
                                        <input type="number" class="form-input" id="inlet-temp" 
                                               value="298" min="200" max="500" step="1">
                                        <span class="input-unit">K</span>
                                    </div>
                                    <span class="text-xs text-gray-500" id="inlet-temp-celsius">= 25.0°C</span>
                                </div>
                                
                                <!-- Inlet Pressure -->
                                <div class="form-group">
                                    <label class="form-label" for="inlet-pressure" data-i18n="compressor.inletPressure">Pression d'entrée</label>
                                    <div class="input-with-unit">
                                        <input type="number" class="form-input" id="inlet-pressure"
                                               value="101.3" min="10" max="1000" step="0.1">
                                        <span class="input-unit">kPa</span>
                                    </div>
                                </div>
                                
                                <!-- Pressure Ratio -->
                                <div class="form-group" style="grid-column: span 2;">
                                    <label class="form-label" for="pressure-ratio"><span data-i18n="compressor.pressureRatio">Taux de compression</span> (r): <span id="pressure-ratio-value" class="font-mono">3.0</span></label>
                                    <input type="range" class="form-range w-full" id="pressure-ratio"
                                           value="3" min="1.1" max="10" step="0.1">
                                </div>
                                
                                <!-- Mass Flow -->
                                <div class="form-group">
                                    <label class="form-label" for="mass-flow" data-i18n="common.massFlow">Débit massique</label>
                                    <div class="input-with-unit">
                                        <input type="number" class="form-input" id="mass-flow"
                                               value="1.0" min="0.1" max="100" step="0.1">
                                        <span class="input-unit">kg/s</span>
                                    </div>
                                </div>
                                
                                <!-- Efficiency -->
                                <div class="form-group">
                                    <label class="form-label" for="efficiency"><span data-i18n="compressor.isentropicEff">Rendement isentropique</span> (ηₛ): <span id="efficiency-value" class="font-mono">85%</span></label>
                                    <input type="range" class="form-range w-full" id="efficiency"
                                           value="85" min="50" max="100" step="1">
                                </div>
                            </div>
                            
                            <!-- Action Buttons -->
                            <div class="mt-4 flex flex-wrap gap-2">
                                <button class="btn btn-primary flex-1" id="run-simulation">
                                    ▶️ <span data-i18n="common.run">Calculer</span>
                                </button>
                                <button class="btn btn-secondary flex-1" id="reset-params">
                                    🔄 <span data-i18n="common.reset">Réinitialiser</span>
                                </button>
                            </div>
                            <div class="mt-2 flex flex-wrap gap-2">
                                <button class="btn btn-outline flex-1" id="run-transient">
                                    📈 <span data-i18n="compressor.runTransient">Simulation transitoire</span>
                                </button>
                                <button class="btn btn-outline flex-1" id="load-example">
                                    💡 <span data-i18n="common.example">Exemple</span>
                                </button>
                            </div>
                        </div>
                    </div>
                </div>
                
                <!-- Right Column: Results & Charts -->
                <div style="flex: 1; min-width: 300px; display: flex; flex-direction: column; gap: 20px;">
                    
                    <!-- Results Card -->
                    <div class="card">
                        <div class="card-header">
                            <h3 class="card-title" style="display: flex; align-items: center; gap: 8px;">
                                <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <path d="M22 12h-4l-3 9L9 3l-3 9H2"/>
                                </svg>
                                Résultats
                            </h3>
                        </div>
                        <div class="card-body">
                            <div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(140px, 1fr)); gap: 1rem;">
                                <div class="p-3 bg-gray-50 rounded text-center border">
                                    <div class="text-sm text-gray-500 mb-1">T₂ (<span data-i18n="common.outlet">sortie</span>)</div>
                                    <div class="text-xl font-bold text-blue-600" id="result-t2">--</div>
                                    <div class="text-xs text-gray-500">K</div>
                                </div>
                                <div class="p-3 bg-gray-50 rounded text-center border">
                                    <div class="text-sm text-gray-500 mb-1">P₂ (<span data-i18n="common.outlet">sortie</span>)</div>
                                    <div class="text-xl font-bold text-blue-600" id="result-p2">--</div>
                                    <div class="text-xs text-gray-500">kPa</div>
                                </div>
                                <div class="p-3 bg-gray-50 rounded text-center border">
                                    <div class="text-sm text-gray-500 mb-1" data-i18n="common.power">Puissance</div>
                                    <div class="text-xl font-bold text-green-600" id="result-power">--</div>
                                    <div class="text-xs text-gray-500">kW</div>
                                </div>
                                <div class="p-3 bg-gray-50 rounded text-center border">
                                    <div class="text-sm text-gray-500 mb-1">T₂ₛ (isentropique)</div>
                                    <div class="text-xl font-bold text-gray-700" id="result-t2s">--</div>
                                    <div class="text-xs text-gray-500">K</div>
                                </div>
                                <div class="p-3 bg-gray-50 rounded text-center border">
                                    <div class="text-sm text-gray-500 mb-1">ΔT</div>
                                    <div class="text-xl font-bold text-gray-700" id="result-dt">--</div>
                                    <div class="text-xs text-gray-500">K</div>
                                </div>
                                <div class="p-3 bg-gray-50 rounded text-center border">
                                    <div class="text-sm text-gray-500 mb-1">Δs</div>
                                    <div class="text-xl font-bold text-gray-700" id="result-ds">--</div>
                                    <div class="text-xs text-gray-500">J/(kg·K)</div>
                                </div>
                            </div>
                            
                            <div class="mt-4">
                                <!-- Detailed calculations button removed -->
                            </div>
                            
                            <div class="mt-4 hidden p-4 bg-gray-50 rounded border" id="calculation-steps">
                            </div>
                        </div>
                    </div>
                    
                    <!-- Charts Card -->
                    <div class="card">
                        <div class="card-header" style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 10px;">
                            <h3 class="card-title">Graphiques</h3>
                            <div class="btn-group" role="group">
                                <button class="btn btn-sm btn-outline active mode-btn" data-chart="power">Puissance</button>
                                <button class="btn btn-sm btn-outline mode-btn" data-chart="ts">T-s</button>
                                <button class="btn btn-sm btn-outline mode-btn" data-chart="pv">P-v</button>
                            </div>
                        </div>
                        <div class="card-body">
                            <div class="chart-wrapper" style="position: relative; height: 300px; width: 100%;">
                                <div id="power-chart-container" class="chart-container" style="height: 100%; width: 100%;">
                                    <div class="flex items-center justify-center h-full text-gray-500">
                                        Lancez une simulation pour voir le graphique
                                    </div>
                                </div>
                                <div id="ts-chart-container" class="chart-container hidden" style="height: 100%; width: 100%;"></div>
                                <div id="pv-chart-container" class="chart-container hidden" style="height: 100%; width: 100%;"></div>
                            </div>
                        </div>
                    </div>
                    
                </div>
            </div>
            
            <!-- Transient Section -->
            <div class="card mt-5 hidden" id="transient-section">
                <div class="card-header" style="display: flex; justify-content: space-between; align-items: center;">
                    <h3 class="card-title">📈 Simulation Transitoire</h3>
                    <button class="btn btn-icon" id="close-transient">&times;</button>
                </div>
                <div class="card-body">
                    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;">
                        <div style="height: 300px; position: relative;">
                            <div id="transient-flow-chart" style="height: 100%; width: 100%;"></div>
                        </div>
                        <div style="height: 300px; position: relative;">
                            <div id="transient-temp-chart" style="height: 100%; width: 100%;"></div>
                        </div>
                    </div>
                    <div style="height: 300px; position: relative; margin-top: 20px;">
                        <div id="transient-power-chart" style="height: 100%; width: 100%;"></div>
                    </div>
                </div>
            </div>
            
        </div>
    `;
};

/**
 * Initialize the module after rendering
 */
CompressorModule.init = async function() {
    // Load SVG diagram
    await loadCompressorDiagram();
    
    // Setup event listeners
    setupEventListeners();
    
    // Run initial calculation
    runSimulation();
};

/**
 * Load and display the compressor SVG diagram
 */
async function loadCompressorDiagram() {
    const container = document.getElementById('compressor-diagram');
    try {
        const response = await fetch('assets/compressor.svg');
        if (response.ok) {
            container.innerHTML = await response.text();
        } else {
            container.innerHTML = '<p class="text-muted">Schéma non disponible</p>';
        }
    } catch (error) {
        console.warn('Could not load compressor SVG:', error);
        container.innerHTML = '<p class="text-muted">Schéma non disponible</p>';
    }
}

/**
 * Setup all event listeners for the module
 */
function setupEventListeners() {
    // Input change handlers
    document.getElementById('inlet-temp')?.addEventListener('input', (e) => {
        const kelvin = parseFloat(e.target.value);
        document.getElementById('inlet-temp-celsius').textContent = 
            `= ${(kelvin - 273.15).toFixed(1)}°C`;
        CompressorModule.params.inletTemperature = kelvin;
    });
    
    document.getElementById('inlet-pressure')?.addEventListener('input', (e) => {
        CompressorModule.params.inletPressure = parseFloat(e.target.value) * 1000; // kPa to Pa
    });
    
    document.getElementById('pressure-ratio')?.addEventListener('input', (e) => {
        const value = parseFloat(e.target.value);
        document.getElementById('pressure-ratio-value').textContent = value.toFixed(1);
        CompressorModule.params.pressureRatio = value;
    });
    
    document.getElementById('mass-flow')?.addEventListener('input', (e) => {
        CompressorModule.params.massFlow = parseFloat(e.target.value);
    });
    
    document.getElementById('efficiency')?.addEventListener('input', (e) => {
        const value = parseInt(e.target.value);
        document.getElementById('efficiency-value').textContent = `${value}%`;
        CompressorModule.params.isentropicEfficiency = value / 100;
    });
    
    document.getElementById('gas-type')?.addEventListener('change', (e) => {
        CompressorModule.params.gasType = e.target.value;
    });
    
    document.getElementById('compressor-type')?.addEventListener('change', (e) => {
        CompressorModule.params.compressorType = e.target.value;
    });
    
    // Run simulation button
    document.getElementById('run-simulation')?.addEventListener('click', runSimulation);
    
    // Run transient simulation
    document.getElementById('run-transient')?.addEventListener('click', runTransient);
    
    // Load example
    document.getElementById('load-example')?.addEventListener('click', loadExample);
    
    // Reset parameters
    document.getElementById('reset-params')?.addEventListener('click', resetParameters);
    
    // Toggle calculations listener removed
    
    // Close transient section
    document.getElementById('close-transient')?.addEventListener('click', () => {
        document.getElementById('transient-section').classList.add('hidden');
    });

    // Visualize button removed
    /*
    document.getElementById('visualize-btn')?.addEventListener('click', () => {
        // Open the visualization modal
        const modal = document.createElement('div');
        modal.className = 'modal-overlay open';
        modal.id = 'visualize-modal';
        modal.innerHTML = `
            <div class="modal" style="max-width: 800px;">
                <div class="modal-header">
                    <h3 class="modal-title">Visualisation en temps réel</h3>
                    <button class="modal-close" onclick="document.getElementById('visualize-modal').remove()">&times;</button>
                </div>
                <div class="modal-body">
                    <div style="position: relative; width: 100%; height: 400px; background: #f8f9fa; border-radius: 8px; overflow: hidden;">
                        <canvas id="viz-canvas" style="width: 100%; height: 100%;"></canvas>
                        <div id="viz-stats" style="position: absolute; top: 10px; right: 10px; background: rgba(255,255,255,0.9); padding: 10px; border-radius: 4px; font-family: monospace; font-size: 12px; border: 1px solid #ddd;">
                            Ready
                        </div>
                    </div>
                    <div style="margin-top: 20px; text-align: center;">
                        <button class="btn btn-primary" id="viz-start-btn">Démarrer l'animation</button>
                    </div>
                </div>
            </div>
        `;
        document.body.appendChild(modal);

        // Animation logic
        const canvas = document.getElementById('viz-canvas');
        const ctx = canvas.getContext('2d');
        const stats = document.getElementById('viz-stats');
        
        // Set canvas resolution
        const dpr = window.devicePixelRatio || 1;
        canvas.width = canvas.offsetWidth * dpr;
        canvas.height = canvas.offsetHeight * dpr;
        ctx.scale(dpr, dpr);
        
        let animationId;
        let time = 0;
        let isAnimating = false;
        
        // Compressor parameters
        const cx = canvas.offsetWidth / 2;
        const cy = canvas.offsetHeight / 2;
        
        // Particles
        const particles = [];
        const numParticles = 50;
        
        function initParticles() {
            particles.length = 0;
            for(let i=0; i<numParticles; i++) {
                particles.push({
                    x: Math.random() * canvas.offsetWidth,
                    y: Math.random() * canvas.offsetHeight,
                    vx: (Math.random() - 0.5) * 2,
                    vy: (Math.random() - 0.5) * 2,
                    stage: 0 // 0: inlet, 1: compression, 2: outlet
                });
            }
        }
        
        function drawCompressor(rotation) {
            const w = canvas.offsetWidth;
            const h = canvas.offsetHeight;
            
            ctx.clearRect(0, 0, w, h);
            
            // Draw housing
            ctx.beginPath();
            ctx.strokeStyle = '#333';
            ctx.lineWidth = 4;
            ctx.moveTo(cx - 100, cy + 60);
            ctx.lineTo(cx - 60, cy);
            ctx.lineTo(cx + 60, cy);
            ctx.lineTo(cx + 100, cy + 60);
            ctx.stroke();
            
            // Draw rotor/shaft
            ctx.fillStyle = '#666';
            ctx.fillRect(cx - 120, cy + 50, 240, 10);
            
            // Draw blades
            ctx.save();
            ctx.translate(cx, cy + 20);
            
            // Draw static vanes (stator)
            ctx.strokeStyle = '#999';
            ctx.lineWidth = 2;
            for(let i=0; i<8; i++) {
                ctx.beginPath();
                ctx.moveTo(0, 0);
                const a = (i * Math.PI * 2) / 8;
                ctx.lineTo(Math.cos(a) * 40, Math.sin(a) * 40);
                ctx.stroke();
            }
            
            // Draw rotating blades
            ctx.rotate(rotation);
            ctx.strokeStyle = '#008080';
            ctx.lineWidth = 3;
            for(let i=0; i<6; i++) {
                ctx.beginPath();
                ctx.moveTo(0, 0);
                const a = (i * Math.PI * 2) / 6;
                ctx.lineTo(Math.cos(a) * 50, Math.sin(a) * 50);
                ctx.stroke();
            }
            ctx.restore();
            
            // Draw gas flow (simplified)
            // Left (Inlet) - Blue
            ctx.fillStyle = 'rgba(52, 152, 219, 0.2)';
            ctx.fillRect(0, 0, cx - 100, h);
            
            // Right (Outlet) - Red
            ctx.fillStyle = 'rgba(231, 76, 60, 0.2)';
            ctx.fillRect(cx + 100, 0, w - (cx + 100), h);
        }
        
        function update() {
            if (!isAnimating) return;
            
            time += 0.1;
            drawCompressor(time);
            
            // Update particles
            // Inlet -> Compressor -> Outlet
            const w = canvas.offsetWidth;
            const h = canvas.offsetHeight;
            
            particles.forEach(p => {
                // Determine stage and color
                let color = '#3498db'; // Blue (cold)
                
                // Move towards right
                p.x += 2 + Math.random(); 
                
                if (p.x > cx - 60 && p.x < cx + 60) {
                    // Inside compressor
                    color = '#9b59b6'; // Purple (transition)
                    p.y += Math.sin(time * 5 + p.x) * 2; // Turbulence
                } else if (p.x >= cx + 60) {
                    // Outlet
                    color = '#e74c3c'; // Red (hot)
                    p.x += 3; // Accelerate
                }
                
                // Respawn
                if (p.x > w) {
                    p.x = 0;
                    p.y = Math.random() * h;
                }
                
                // Draw particle
                ctx.beginPath();
                ctx.fillStyle = color;
                ctx.arc(p.x, p.y, 4, 0, Math.PI * 2);
                ctx.fill();
            });
            
            // Update stats
            const temp = 25 + Math.min(time * 2, 100); // Fake temp rise
            stats.innerHTML = `
                Temps: ${(time/10).toFixed(1)} s<br>
                RPM: ${(2000 + Math.sin(time)*100).toFixed(0)}<br>
                Temp. Sortie: ${(25 + (time > 10 ? 150 : time*15)).toFixed(1)} °C
            `;
            
            animationId = requestAnimationFrame(update);
        }
        
        // Initial draw
        initParticles();
        drawCompressor(0);
        
        document.getElementById('viz-start-btn').addEventListener('click', () => {
            if (isAnimating) {
                isAnimating = false;
                document.getElementById('viz-start-btn').textContent = "Reprendre l'animation";
                cancelAnimationFrame(animationId);
            } else {
                isAnimating = true;
                document.getElementById('viz-start-btn').textContent = "Pause";
                update();
            }
        });
        */
       // Deprecated visualization code removed
    // });
    
    // Chart type switching
    document.querySelectorAll('[data-chart]').forEach(btn => {
        btn.addEventListener('click', (e) => {
            const chartType = e.target.dataset.chart;
            
            // Update button states
            document.querySelectorAll('[data-chart]').forEach(b => b.classList.remove('active'));
            e.target.classList.add('active');
            
            // Show/hide chart containers
            document.getElementById('power-chart-container').classList.toggle('hidden', chartType !== 'power');
            document.getElementById('ts-chart-container').classList.toggle('hidden', chartType !== 'ts');
            document.getElementById('pv-chart-container').classList.toggle('hidden', chartType !== 'pv');
        });
    });
}

/**
 * Run the main simulation
 */
function runSimulation() {
    console.log("Starting runSimulation...");
    try {
        const statusBadge = document.getElementById('simulation-status');
        if (statusBadge) {
            statusBadge.className = 'status-badge running';
            statusBadge.innerHTML = '<span class="status-dot"></span><span>Calcul...</span>';
        }
        
        // Get current parameters from UI
        updateParamsFromUI();
        console.log("Params:", CompressorModule.params);
        
        // Perform calculation
        const results = calculateCompressor(CompressorModule.params);
        console.log("Results from calculateCompressor:", results);
        
        CompressorModule.results = results;
        
        // Explicitly handle error reported by calculateCompressor
        if (results.error) {
            showToast("Calcul échoué: " + results.message, 'error');
        } else if (isNaN(results.outletTemperature)) {
            showToast("Erreur: Résultats invalides (NaN)", 'error');
        }

        // Update results display
        updateResultsDisplay(results);
        
        // Update diagram values
        updateDiagramValues(results);
        
        // Generate and display charts
        updateCharts(results);
        
        // Update calculation steps
        // updateCalculationSteps(results.calculationSteps); // Feature removed
        
        // Update status
        if (statusBadge) {
            statusBadge.className = 'status-badge idle';
            statusBadge.innerHTML = '<span class="status-dot"></span><span>Terminé</span>';
        }
        
        if (!results.error && typeof showToast === 'function') {
            showToast(t('toast.simulationComplete'), 'success');
        }
    } catch (e) {
        console.error("Error in runSimulation:", e);
        if (typeof showToast === 'function') {
            showToast('Erreur critique: ' + e.message, 'error');
        }
    }
}

/**
 * Run transient simulation
 */
function runTransient() {
    const statusBadge = document.getElementById('simulation-status');
    statusBadge.className = 'status-badge running';
    statusBadge.innerHTML = '<span class="status-dot"></span><span>Simulation...</span>';
    
    updateParamsFromUI();
    
    // Show transient section
    document.getElementById('transient-section').classList.remove('hidden');
    
    // Run transient simulation
    const transientResults = runTransientSimulation(CompressorModule.params);
    CompressorModule.transientData = transientResults;
    
    // Plot transient results
    plotTransientResults(transientResults);
    
    statusBadge.className = 'status-badge idle';
    statusBadge.innerHTML = '<span class="status-dot"></span><span>Terminé</span>';
    
    showToast('Simulation transitoire terminée', 'success');
}

/**
 * Load example parameters
 */
function loadExample() {
    // Air compressor 1.2 MPa outlet example
    const example = {
        inletTemperature: 298.15,      // 25°C
        inletPressure: 101325,         // 1 atm
        pressureRatio: 3.0,            // Approx 303.9 kPa → ~0.3 MPa (for 1.2 MPa, use ratio ~12)
        massFlow: 1.5,
        isentropicEfficiency: 0.82,
        gasType: 'air',
        compressorType: 'centrifugal'
    };
    
    // For 1.2 MPa outlet from 101.3 kPa:
    example.pressureRatio = 1200 / 101.3;  // ≈ 11.85
    
    CompressorModule.params = example;
    
    // Update UI
    document.getElementById('inlet-temp').value = example.inletTemperature.toFixed(0);
    document.getElementById('inlet-temp-celsius').textContent = 
        `= ${(example.inletTemperature - 273.15).toFixed(1)}°C`;
    document.getElementById('inlet-pressure').value = (example.inletPressure / 1000).toFixed(1);
    document.getElementById('pressure-ratio').value = example.pressureRatio;
    document.getElementById('pressure-ratio-value').textContent = example.pressureRatio.toFixed(1);
    document.getElementById('mass-flow').value = example.massFlow;
    document.getElementById('efficiency').value = example.isentropicEfficiency * 100;
    document.getElementById('efficiency-value').textContent = `${(example.isentropicEfficiency * 100).toFixed(0)}%`;
    document.getElementById('gas-type').value = example.gasType;
    document.getElementById('compressor-type').value = example.compressorType;
    
    // Run simulation with example
    runSimulation();
    
    showToast('Exemple chargé: Compresseur d\'air 1.2 MPa', 'info');
}

/**
 * Reset parameters to defaults
 */
function resetParameters() {
    CompressorModule.params = { ...CompressorModule.defaults };
    const p = CompressorModule.params;
    
    document.getElementById('inlet-temp').value = p.inletTemperature.toFixed(0);
    document.getElementById('inlet-temp-celsius').textContent = 
        `= ${(p.inletTemperature - 273.15).toFixed(1)}°C`;
    document.getElementById('inlet-pressure').value = (p.inletPressure / 1000).toFixed(1);
    document.getElementById('pressure-ratio').value = p.pressureRatio;
    document.getElementById('pressure-ratio-value').textContent = p.pressureRatio.toFixed(1);
    document.getElementById('mass-flow').value = p.massFlow;
    document.getElementById('efficiency').value = p.isentropicEfficiency * 100;
    document.getElementById('efficiency-value').textContent = `${(p.isentropicEfficiency * 100).toFixed(0)}%`;
    document.getElementById('gas-type').value = p.gasType;
    document.getElementById('compressor-type').value = p.compressorType;
    
    runSimulation();
    showToast('Paramètres réinitialisés', 'info');
}

// Function toggleCalculations removed

/**
 * Update parameters from UI values
 */
function updateParamsFromUI() {
    CompressorModule.params = {
        inletTemperature: parseFloat(document.getElementById('inlet-temp').value),
        inletPressure: parseFloat(document.getElementById('inlet-pressure').value) * 1000,
        pressureRatio: parseFloat(document.getElementById('pressure-ratio').value),
        massFlow: parseFloat(document.getElementById('mass-flow').value),
        isentropicEfficiency: parseFloat(document.getElementById('efficiency').value) / 100,
        gasType: document.getElementById('gas-type').value,
        compressorType: document.getElementById('compressor-type').value
    };
}

/**
 * Update the results display with calculated values
 */
function updateResultsDisplay(results) {
    // Debug helper to see WHY values are invalid
    const safeFixed = (val, dec) => {
        if (typeof val === 'undefined') return '-- (undef)';
        if (val === null) return '-- (null)';
        if (isNaN(val)) return '-- (NaN)';
        return val.toFixed(dec);
    };
    
    document.getElementById('result-t2').textContent = safeFixed(results.outletTemperature, 1);
    document.getElementById('result-p2').textContent = safeFixed(results.outletPressure / 1000, 1);
    document.getElementById('result-power').textContent = safeFixed(results.actualPower / 1000, 2);
    document.getElementById('result-t2s').textContent = safeFixed(results.outletTemperatureIsentropic, 1);
    document.getElementById('result-dt').textContent = safeFixed(results.temperatureRise, 1);
    document.getElementById('result-ds').textContent = safeFixed(results.entropyChange, 2);
}

/**
 * Update SVG diagram with calculated values
 */
function updateDiagramValues(results) {
    const diagram = document.getElementById('compressor-diagram');
    if (!diagram) return;
    
    // Safety helper
    const safeSet = (selector, val, decimals = 1) => {
        const el = diagram.querySelector(selector);
        if (el && val !== undefined && val !== null && !isNaN(val)) {
             el.textContent = val.toFixed(decimals);
        }
    };
    
    // Update inlet values
    safeSet('.inlet-pressure', results.inputs.inletPressure / 1000, 1);
    safeSet('.inlet-temp', results.inputs.inletTemperature, 0);
    safeSet('.mass-flow', results.inputs.massFlow, 2);
    
    // Update outlet values
    safeSet('.outlet-pressure', results.outletPressure / 1000, 1);
    safeSet('.outlet-temp', results.outletTemperature, 0);
    safeSet('.power', results.actualPower / 1000, 1);
    
    // Update efficiency
    safeSet('.efficiency-value', results.inputs.isentropicEfficiency * 100, 0);
}

/**
 * Update calculation steps display
 */
function updateCalculationSteps(steps) {
    const container = document.getElementById('calculation-steps');
    
    container.innerHTML = steps.map((step, index) => `
        <div class="calc-step">
            <div class="step-number">${index + 1}</div>
            <div class="step-content">
                <div class="step-description">${step.description}</div>
                <div class="step-formula">${step.formula}</div>
                <div class="step-result">${step.result}</div>
            </div>
        </div>
    `).join('');
}

/**
 * Update all charts with new results
 */
function updateCharts(results) {
    // Power vs Pressure Ratio curve
    const powerCurve = generatePowerCurve(CompressorModule.params);
    plotPowerChart(powerCurve, results);
    
    // T-s diagram
    const tsData = generateTsDiagram(results);
    plotTsDiagram(tsData, results);
    
    // P-v diagram
    const pvData = generatePvDiagram(results);
    plotPvDiagram(pvData, results);
}

/**
 * Plot power vs pressure ratio chart
 */
function plotPowerChart(curve, results) {
    const container = document.getElementById('power-chart-container');
    
    // Destroy existing chart if any
    if (CompressorModule.charts.power) {
        CompressorModule.charts.power.destroy?.();
    }
    
    const currentRatio = results.inputs.pressureRatio;
    const currentPower = results.actualPower / 1000;
    
    const config = {
        type: 'line',
        data: {
            labels: curve.map(p => p.pressureRatio.toFixed(1)),
            datasets: [
                {
                    label: 'Puissance (kW)',
                    data: curve.map(p => p.power),
                    borderColor: '#008080',
                    backgroundColor: 'rgba(0, 128, 128, 0.1)',
                    fill: true,
                    tension: 0.3
                },
                {
                    label: 'Point actuel',
                    data: curve.map(p => Math.abs(p.pressureRatio - currentRatio) < 0.2 ? currentPower : null),
                    borderColor: '#e74c3c',
                    backgroundColor: '#e74c3c',
                    pointRadius: 8,
                    pointHoverRadius: 10,
                    showLine: false
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: {
                    display: true,
                    text: 'Puissance vs Taux de Compression'
                },
                legend: {
                    display: true,
                    position: 'top'
                }
            },
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'Taux de compression (r = P₂/P₁)'
                    }
                },
                y: {
                    title: {
                        display: true,
                        text: 'Puissance (kW)'
                    },
                    beginAtZero: true
                }
            }
        }
    };
    
    CompressorModule.charts.power = createChart('power-chart-container', config);
}

/**
 * Plot T-s diagram
 */
function plotTsDiagram(data, results) {
    const container = document.getElementById('ts-chart-container');
    
    if (CompressorModule.charts.ts) {
        CompressorModule.charts.ts.destroy?.();
    }
    
    const config = {
        type: 'scatter',
        data: {
            datasets: [
                {
                    label: 'Processus isentropique',
                    data: data.isentropicCurve.map(p => ({ x: p.s, y: p.T })),
                    borderColor: '#27ae60',
                    backgroundColor: '#27ae60',
                    showLine: true,
                    borderDash: [5, 5],
                    pointRadius: 0
                },
                {
                    label: 'Processus réel',
                    data: data.actualCurve.map(p => ({ x: p.s, y: p.T })),
                    borderColor: '#e74c3c',
                    backgroundColor: 'rgba(231, 76, 60, 0.1)',
                    showLine: true,
                    fill: true,
                    pointRadius: 0
                },
                {
                    label: 'États',
                    data: [
                        { x: 0, y: results.inputs.inletTemperature },
                        { x: 0, y: results.outletTemperatureIsentropic },
                        { x: data.actualCurve[data.actualCurve.length - 1].s, y: results.outletTemperature }
                    ],
                    borderColor: '#333',
                    backgroundColor: ['#008080', '#27ae60', '#e74c3c'],
                    pointRadius: 8,
                    showLine: false
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: {
                    display: true,
                    text: 'Diagramme T-s'
                }
            },
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'Δs (J/kg·K)'
                    }
                },
                y: {
                    title: {
                        display: true,
                        text: 'Température (K)'
                    }
                }
            }
        }
    };
    
    CompressorModule.charts.ts = createChart('ts-chart-container', config);
}

/**
 * Plot P-v diagram
 */
function plotPvDiagram(data, results) {
    const container = document.getElementById('pv-chart-container');
    
    if (CompressorModule.charts.pv) {
        CompressorModule.charts.pv.destroy?.();
    }
    
    const config = {
        type: 'scatter',
        data: {
            datasets: [
                {
                    label: 'Compression isentropique',
                    data: data.isentropicCurve.map(p => ({ x: p.v * 1000, y: p.P / 1000 })),
                    borderColor: '#27ae60',
                    showLine: true,
                    borderDash: [5, 5],
                    pointRadius: 0
                },
                {
                    label: 'Compression polytropique',
                    data: data.polytropicCurve.map(p => ({ x: p.v * 1000, y: p.P / 1000 })),
                    borderColor: '#e74c3c',
                    showLine: true,
                    pointRadius: 0
                },
                {
                    label: 'États',
                    data: [
                        { x: results.inletSpecificVolume * 1000, y: results.inputs.inletPressure / 1000 },
                        { x: results.outletSpecificVolume * 1000, y: results.outletPressure / 1000 }
                    ],
                    borderColor: '#333',
                    backgroundColor: ['#008080', '#e74c3c'],
                    pointRadius: 8,
                    showLine: false
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: {
                    display: true,
                    text: 'Diagramme P-v'
                }
            },
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'Volume spécifique (L/kg)'
                    }
                },
                y: {
                    title: {
                        display: true,
                        text: 'Pression (kPa)'
                    }
                }
            }
        }
    };
    
    CompressorModule.charts.pv = createChart('pv-chart-container', config);
}

/**
 * Plot transient simulation results
 */
function plotTransientResults(transientResults) {
    const { data, times } = transientResults;
    
    // Mass flow chart
    if (CompressorModule.charts.transientFlow) {
        CompressorModule.charts.transientFlow.destroy?.();
    }
    
    CompressorModule.charts.transientFlow = createChart('transient-flow-chart', {
        type: 'line',
        data: {
            labels: data.map(d => d.time.toFixed(1)),
            datasets: [{
                label: 'Débit massique (kg/s)',
                data: data.map(d => d.massFlow),
                borderColor: '#008080',
                backgroundColor: 'rgba(0, 128, 128, 0.1)',
                fill: true,
                tension: 0.3
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: { display: true, text: 'Évolution du débit massique' }
            },
            scales: {
                x: { title: { display: true, text: 'Temps (s)' } },
                y: { title: { display: true, text: 'Débit (kg/s)' }, beginAtZero: true }
            }
        }
    });
    
    // Temperature chart
    if (CompressorModule.charts.transientTemp) {
        CompressorModule.charts.transientTemp.destroy?.();
    }
    
    CompressorModule.charts.transientTemp = createChart('transient-temp-chart', {
        type: 'line',
        data: {
            labels: data.map(d => d.time.toFixed(1)),
            datasets: [{
                label: 'Température de sortie (K)',
                data: data.map(d => d.outletTemperature),
                borderColor: '#e74c3c',
                backgroundColor: 'rgba(231, 76, 60, 0.1)',
                fill: true,
                tension: 0.3
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: { display: true, text: 'Évolution de la température de sortie' }
            },
            scales: {
                x: { title: { display: true, text: 'Temps (s)' } },
                y: { title: { display: true, text: 'Température (K)' } }
            }
        }
    });
    
    // Power chart
    if (CompressorModule.charts.transientPower) {
        CompressorModule.charts.transientPower.destroy?.();
    }
    
    CompressorModule.charts.transientPower = createChart('transient-power-chart', {
        type: 'line',
        data: {
            labels: data.map(d => d.time.toFixed(1)),
            datasets: [{
                label: 'Puissance (kW)',
                data: data.map(d => d.power),
                borderColor: '#27ae60',
                backgroundColor: 'rgba(39, 174, 96, 0.1)',
                fill: true,
                tension: 0.3
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: { display: true, text: 'Évolution de la puissance' }
            },
            scales: {
                x: { title: { display: true, text: 'Temps (s)' } },
                y: { title: { display: true, text: 'Puissance (kW)' }, beginAtZero: true }
            }
        }
    });
}

// ============================================================================
// EXPLANATION PANEL CONTENT
// ============================================================================

CompressorModule.getExplanation = function(lang = 'fr') {
    if (lang === 'fr') {
        return {
            title: 'Simulation d\'un compresseur',
            description: `
                Un compresseur est une machine qui augmente la pression d'un gaz en effectuant 
                un travail mécanique. Cette simulation modélise un compresseur mono-étagé 
                (centrifuge ou à piston) avec prise en compte du rendement isentropique.
            `,
            theory: `
                <p><strong>Compression isentropique:</strong> Processus réversible et adiabatique 
                (pas d'échange de chaleur avec l'extérieur). C'est le cas idéal.</p>
                
                <p><strong>Compression réelle:</strong> Les irréversibilités (frottements, 
                turbulence) augmentent l'entropie et la température de sortie.</p>
                
                <p><strong>Rendement isentropique:</strong> Rapport entre le travail isentropique 
                (idéal) et le travail réel. Typiquement 70-90% pour les compresseurs industriels.</p>
            `,
            formulas: `
                <p><strong>Température isentropique:</strong></p>
                $$ T_{2s} = T_1 \\left(\\frac{P_2}{P_1}\\right)^{\\frac{\\gamma-1}{\\gamma}} $$
                
                <p><strong>Température réelle:</strong></p>
                $$ T_2 = T_1 + \\frac{T_{2s} - T_1}{\\eta_s} $$
                
                <p><strong>Puissance:</strong></p>
                $$ \\dot{W} = \\dot{m} \\cdot c_p \\cdot (T_2 - T_1) $$
                
                <p><strong>Relation polytropique:</strong></p>
                $$ P \\cdot v^n = \\text{constante} $$
            `,
            references: [
                'Çengel, Y.A. & Boles, M.A. "Thermodynamics: An Engineering Approach"',
                'Perry\'s Chemical Engineers\' Handbook, 8th Edition',
                'Smith, Van Ness, Abbott "Introduction to Chemical Engineering Thermodynamics"'
            ]
        };
    }
    
    // English
    return {
        title: 'Compressor Simulation',
        description: `
            A compressor is a machine that increases gas pressure by performing mechanical work. 
            This simulation models a single-stage compressor (centrifugal or piston) with 
            isentropic efficiency consideration.
        `,
        theory: `
            <p><strong>Isentropic compression:</strong> Reversible and adiabatic process 
            (no heat exchange with surroundings). This is the ideal case.</p>
            
            <p><strong>Real compression:</strong> Irreversibilities (friction, turbulence) 
            increase entropy and outlet temperature.</p>
            
            <p><strong>Isentropic efficiency:</strong> Ratio of isentropic (ideal) work 
            to actual work. Typically 70-90% for industrial compressors.</p>
        `,
        formulas: `
            <p><strong>Isentropic temperature:</strong></p>
            $$ T_{2s} = T_1 \\left(\\frac{P_2}{P_1}\\right)^{\\frac{\\gamma-1}{\\gamma}} $$
            
            <p><strong>Actual temperature:</strong></p>
            $$ T_2 = T_1 + \\frac{T_{2s} - T_1}{\\eta_s} $$
            
            <p><strong>Power:</strong></p>
            $$ \\dot{W} = \\dot{m} \\cdot c_p \\cdot (T_2 - T_1) $$
        `,
        references: [
            'Çengel, Y.A. & Boles, M.A. "Thermodynamics: An Engineering Approach"',
            'Perry\'s Chemical Engineers\' Handbook, 8th Edition'
        ]
    };
};

// ============================================================================
// EXPORT FUNCTIONALITY
// ============================================================================

CompressorModule.getExportData = function() {
    if (!this.results) return null;
    
    const r = this.results;
    const p = r.inputs;
    
    return {
        csv: {
            headers: ['Parameter', 'Value', 'Unit'],
            rows: [
                ['Inlet Temperature', p.inletTemperature.toFixed(2), 'K'],
                ['Inlet Pressure', (p.inletPressure / 1000).toFixed(2), 'kPa'],
                ['Pressure Ratio', p.pressureRatio.toFixed(2), '-'],
                ['Mass Flow', p.massFlow.toFixed(2), 'kg/s'],
                ['Isentropic Efficiency', (p.isentropicEfficiency * 100).toFixed(1), '%'],
                ['Gas Type', p.gasType, '-'],
                ['Outlet Temperature (Actual)', r.outletTemperature.toFixed(2), 'K'],
                ['Outlet Temperature (Isentropic)', r.outletTemperatureIsentropic.toFixed(2), 'K'],
                ['Outlet Pressure', (r.outletPressure / 1000).toFixed(2), 'kPa'],
                ['Power Required', (r.actualPower / 1000).toFixed(2), 'kW'],
                ['Isentropic Power', (r.isentropicPower / 1000).toFixed(2), 'kW'],
                ['Polytropic Exponent', r.polytropicExponent.toFixed(4), '-'],
                ['Entropy Change', r.entropyChange.toFixed(2), 'J/(kg·K)']
            ]
        },
        report: {
            title: 'Simulation de Compresseur',
            module: 'Compresseur mono-étagé',
            inputs: [
                { name: 'Température d\'entrée', value: p.inletTemperature.toFixed(2), unit: 'K' },
                { name: 'Pression d\'entrée', value: (p.inletPressure / 1000).toFixed(2), unit: 'kPa' },
                { name: 'Taux de compression', value: p.pressureRatio.toFixed(2), unit: '-' },
                { name: 'Débit massique', value: p.massFlow.toFixed(2), unit: 'kg/s' },
                { name: 'Rendement isentropique', value: (p.isentropicEfficiency * 100).toFixed(1), unit: '%' },
                { name: 'Type de gaz', value: p.gasType, unit: '-' }
            ],
            outputs: [
                { name: 'Température de sortie', value: r.outletTemperature.toFixed(2), unit: 'K' },
                { name: 'Pression de sortie', value: (r.outletPressure / 1000).toFixed(2), unit: 'kPa' },
                { name: 'Puissance requise', value: (r.actualPower / 1000).toFixed(2), unit: 'kW' },
                { name: 'Coefficient polytropique', value: r.polytropicExponent.toFixed(4), unit: '-' }
            ],
            calculations: r.calculationSteps
        }
    };
};

CompressorModule.getDefaultParameters = function() {
    return { ...this.defaults };
};

CompressorModule.onModeChange = function(mode) {
    // Handle mode changes (student/assignment/teacher)
    console.log('Compressor module: mode changed to', mode);
    
    // In teacher mode, we might show additional controls or solutions
    // TODO: Implement mode-specific UI changes
};

// ============================================================================
// UNIT TESTS (Run in console for validation)
// ============================================================================

/**
 * Run unit tests for the compressor module
 * Call from console: CompressorModule.runTests()
 */
CompressorModule.runTests = function() {
    console.log('🧪 Running Compressor Module Tests...\n');
    
    // Test 1: Basic calculation with known values
    const test1Params = {
        inletTemperature: 300,      // K
        inletPressure: 100000,      // Pa
        pressureRatio: 3,
        massFlow: 1,
        isentropicEfficiency: 1.0,  // Ideal case
        gasType: 'air'
    };
    
    const test1Result = calculateCompressor(test1Params);
    const gamma = 1.4;
    const expectedT2s = 300 * Math.pow(3, (gamma - 1) / gamma);
    
    console.log('Test 1: Isentropic compression (η = 1)');
    console.log(`  Expected T2s: ${expectedT2s.toFixed(2)} K`);
    console.log(`  Calculated T2s: ${test1Result.outletTemperatureIsentropic.toFixed(2)} K`);
    console.log(`  Match: ${Math.abs(test1Result.outletTemperatureIsentropic - expectedT2s) < 0.1 ? '✅' : '❌'}\n`);
    
    // Test 2: With efficiency
    const test2Params = { ...test1Params, isentropicEfficiency: 0.85 };
    const test2Result = calculateCompressor(test2Params);
    const expectedT2 = 300 + (expectedT2s - 300) / 0.85;
    
    console.log('Test 2: Real compression (η = 0.85)');
    console.log(`  Expected T2: ${expectedT2.toFixed(2)} K`);
    console.log(`  Calculated T2: ${test2Result.outletTemperature.toFixed(2)} K`);
    console.log(`  Match: ${Math.abs(test2Result.outletTemperature - expectedT2) < 0.1 ? '✅' : '❌'}\n`);
    
    // Test 3: Power calculation
    const cp = 1005;  // J/kg·K for air
    const expectedPower = 1 * cp * (test2Result.outletTemperature - 300);
    
    console.log('Test 3: Power calculation');
    console.log(`  Expected Power: ${(expectedPower / 1000).toFixed(2)} kW`);
    console.log(`  Calculated Power: ${(test2Result.actualPower / 1000).toFixed(2)} kW`);
    console.log(`  Match: ${Math.abs(test2Result.actualPower - expectedPower) < 100 ? '✅' : '❌'}\n`);
    
    // Test 4: Entropy change for isentropic process
    const test4Params = { ...test1Params, isentropicEfficiency: 1.0 };
    const test4Result = calculateCompressor(test4Params);
    
    console.log('Test 4: Entropy change (isentropic)');
    console.log(`  Expected Δs: ~0 J/(kg·K)`);
    console.log(`  Calculated Δs: ${test4Result.entropyChange.toFixed(4)} J/(kg·K)`);
    console.log(`  Match: ${Math.abs(test4Result.entropyChange) < 1 ? '✅' : '❌'}\n`);
    
    console.log('🏁 Tests complete!');
};

// ============================================================================
// REGISTER MODULE
// ============================================================================

// Register with the main application
if (typeof registerModule !== 'undefined') {
    registerModule('compressor', CompressorModule);
}

// Also make available globally for direct access
window.CompressorModule = CompressorModule;
