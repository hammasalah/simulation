/**
 * ============================================================================
 * MATH UTILITIES - Numerical Integration and Thermodynamic Helpers
 * ============================================================================
 * 
 * This module provides numerical methods and thermodynamic calculations
 * for the Chemical Engineering Lab Simulation Platform.
 * 
 * PHYSICS BACKGROUND:
 * ------------------
 * - Euler Method: First-order numerical integration, simple but less accurate
 *   y(n+1) = y(n) + h * f(t(n), y(n))
 * 
 * - Runge-Kutta 4th Order (RK4): Fourth-order method, much more accurate
 *   k1 = f(t(n), y(n))
 *   k2 = f(t(n) + h/2, y(n) + h*k1/2)
 *   k3 = f(t(n) + h/2, y(n) + h*k2/2)
 *   k4 = f(t(n) + h, y(n) + h*k3)
 *   y(n+1) = y(n) + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
 * 
 * REFERENCES:
 * - Numerical Recipes in C, Press et al.
 * - Perry's Chemical Engineers' Handbook
 * - Thermodynamics: An Engineering Approach, Çengel & Boles
 * 
 * @author Chemical Engineering Lab Simulation Platform
 * @version 1.0.0
 */

'use strict';

// ============================================================================
// NUMERICAL INTEGRATORS
// ============================================================================

/**
 * Euler method for numerical integration of ODEs
 * 
 * Simple first-order method. Fast but less accurate.
 * Recommended only for quick estimates or when computational speed is critical.
 * 
 * NUMERICAL STABILITY:
 * - Conditionally stable for stiff equations
 * - Time step h should satisfy h < 2/|λ| where λ is the largest eigenvalue
 * - For most compressor transients, h = 0.01s is adequate
 * 
 * @param {Function} derivative - Function f(t, y) returning dy/dt
 * @param {number} y0 - Initial value
 * @param {number} t0 - Initial time
 * @param {number} tEnd - End time
 * @param {number} h - Time step
 * @returns {Object} {times: Array, values: Array} - Solution trajectory
 * 
 * @example
 * // Solve dy/dt = -2y with y(0) = 1
 * const result = MathUtils.euler((t, y) => -2 * y, 1, 0, 5, 0.1);
 */
const euler = (derivative, y0, t0, tEnd, h) => {
    const times = [t0];
    const values = [y0];
    
    let t = t0;
    let y = y0;
    
    while (t < tEnd) {
        // Euler step: y(n+1) = y(n) + h * f(t(n), y(n))
        const dydt = derivative(t, y);
        y = y + h * dydt;
        t = t + h;
        
        times.push(t);
        values.push(y);
    }
    
    return { times, values };
};

/**
 * Euler method for systems of ODEs (vector state)
 * 
 * @param {Function} derivatives - Function f(t, Y) returning array of derivatives
 * @param {Array} Y0 - Initial state vector
 * @param {number} t0 - Initial time
 * @param {number} tEnd - End time
 * @param {number} h - Time step
 * @returns {Object} {times: Array, states: Array of Arrays}
 */
const eulerSystem = (derivatives, Y0, t0, tEnd, h) => {
    const times = [t0];
    const states = [Y0.slice()];
    
    let t = t0;
    let Y = Y0.slice();
    const n = Y.length;
    
    while (t < tEnd) {
        const dYdt = derivatives(t, Y);
        
        for (let i = 0; i < n; i++) {
            Y[i] = Y[i] + h * dYdt[i];
        }
        t = t + h;
        
        times.push(t);
        states.push(Y.slice());
    }
    
    return { times, states };
};

/**
 * Runge-Kutta 4th Order Method (RK4)
 * 
 * Fourth-order accurate method - the workhorse of numerical integration.
 * Recommended for production calculations in this simulation platform.
 * 
 * NUMERICAL STABILITY:
 * - More stable than Euler for same step size
 * - For compressor transients, h = 0.05s typically gives < 0.1% error
 * - Adaptive step size can be implemented if needed (see TODO below)
 * 
 * COMPUTATIONAL COST:
 * - 4 function evaluations per step (vs 1 for Euler)
 * - Usually allows 4x larger step size for same accuracy
 * 
 * @param {Function} derivative - Function f(t, y) returning dy/dt
 * @param {number} y0 - Initial value
 * @param {number} t0 - Initial time
 * @param {number} tEnd - End time
 * @param {number} h - Time step
 * @returns {Object} {times: Array, values: Array}
 * 
 * @example
 * // Solve dy/dt = -2y with y(0) = 1 (exact: y = e^(-2t))
 * const result = MathUtils.rk4((t, y) => -2 * y, 1, 0, 5, 0.1);
 */
const rk4 = (derivative, y0, t0, tEnd, h) => {
    const times = [t0];
    const values = [y0];
    
    let t = t0;
    let y = y0;
    
    while (t < tEnd) {
        // RK4 coefficients
        const k1 = derivative(t, y);
        const k2 = derivative(t + h/2, y + h*k1/2);
        const k3 = derivative(t + h/2, y + h*k2/2);
        const k4 = derivative(t + h, y + h*k3);
        
        // Weighted average
        y = y + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        t = t + h;
        
        times.push(t);
        values.push(y);
    }
    
    return { times, values };
};

/**
 * RK4 for systems of ODEs (vector state)
 * 
 * Used for multi-variable problems like coupled mass-energy balances.
 * 
 * @param {Function} derivatives - Function f(t, Y) returning array of derivatives
 * @param {Array} Y0 - Initial state vector
 * @param {number} t0 - Initial time
 * @param {number} tEnd - End time
 * @param {number} h - Time step
 * @returns {Object} {times: Array, states: Array of Arrays}
 */
const rk4System = (derivatives, Y0, t0, tEnd, h) => {
    const times = [t0];
    const states = [Y0.slice()];
    
    let t = t0;
    let Y = Y0.slice();
    const n = Y.length;
    
    while (t < tEnd) {
        // k1
        const k1 = derivatives(t, Y);
        
        // k2
        const Y_k2 = Y.map((y, i) => y + h * k1[i] / 2);
        const k2 = derivatives(t + h/2, Y_k2);
        
        // k3
        const Y_k3 = Y.map((y, i) => y + h * k2[i] / 2);
        const k3 = derivatives(t + h/2, Y_k3);
        
        // k4
        const Y_k4 = Y.map((y, i) => y + h * k3[i]);
        const k4 = derivatives(t + h, Y_k4);
        
        // Update state
        for (let i = 0; i < n; i++) {
            Y[i] = Y[i] + (h/6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
        }
        t = t + h;
        
        times.push(t);
        states.push(Y.slice());
    }
    
    return { times, states };
};

// TODO: Implement adaptive step-size RK4 (RK45 Dormand-Prince)
// WHY NEEDED: For stiff systems or when accuracy requirements vary during simulation
// HOW TO IMPLEMENT: 
// 1. Compute both 4th and 5th order solutions
// 2. Estimate error from difference
// 3. Adjust step size based on error tolerance
// Reference: Dormand, J. R.; Prince, P. J. (1980)

// ============================================================================
// THERMODYNAMIC CONSTANTS AND PROPERTIES
// ============================================================================

/**
 * Universal gas constant in various units
 */
const GAS_CONSTANTS = {
    SI: 8.314,        // J/(mol·K)
    kJ: 8.314e-3,     // kJ/(mol·K)
    cal: 1.987,       // cal/(mol·K)
    BTU: 1.986,       // BTU/(lbmol·R)
    atm_L: 0.08206    // (atm·L)/(mol·K)
};

/**
 * Common gas properties for compressor simulations
 * 
 * Properties include:
 * - R: Specific gas constant (J/kg·K)
 * - cp: Specific heat at constant pressure (J/kg·K)
 * - cv: Specific heat at constant volume (J/kg·K)
 * - gamma: Ratio of specific heats (cp/cv)
 * - M: Molar mass (kg/mol)
 * 
 * Note: Values are for standard conditions (25°C, 1 atm)
 * For more accurate simulations, use temperature-dependent correlations
 */
const GAS_PROPERTIES = {
    air: {
        name: 'Air',
        nameFr: 'Air',
        R: 287,           // J/(kg·K)
        cp: 1005,         // J/(kg·K)
        cv: 718,          // J/(kg·K)
        gamma: 1.4,       // dimensionless
        M: 0.0289         // kg/mol
    },
    nitrogen: {
        name: 'Nitrogen',
        nameFr: 'Azote',
        R: 296.8,
        cp: 1040,
        cv: 743,
        gamma: 1.4,
        M: 0.028
    },
    oxygen: {
        name: 'Oxygen',
        nameFr: 'Oxygène',
        R: 259.8,
        cp: 918,
        cv: 658,
        gamma: 1.395,
        M: 0.032
    },
    carbon_dioxide: {
        name: 'Carbon Dioxide',
        nameFr: 'Dioxyde de carbone',
        R: 188.9,
        cp: 844,
        cv: 655,
        gamma: 1.289,
        M: 0.044
    },
    methane: {
        name: 'Methane',
        nameFr: 'Méthane',
        R: 518.3,
        cp: 2226,
        cv: 1708,
        gamma: 1.303,
        M: 0.016
    },
    steam: {
        name: 'Steam',
        nameFr: 'Vapeur d\'eau',
        R: 461.5,
        cp: 2080,
        cv: 1619,
        gamma: 1.33,
        M: 0.018
    },
    hydrogen: {
        name: 'Hydrogen',
        nameFr: 'Hydrogène',
        R: 4124,
        cp: 14300,
        cv: 10180,
        gamma: 1.405,
        M: 0.002
    },
    helium: {
        name: 'Helium',
        nameFr: 'Hélium',
        R: 2077,
        cp: 5193,
        cv: 3116,
        gamma: 1.667,
        M: 0.004
    },
    ideal: {
        name: 'Ideal Gas',
        nameFr: 'Gaz parfait',
        R: 287,           // Default to air-like
        cp: 1000,
        cv: 713,
        gamma: 1.4,
        M: 0.029
    }
};

// ============================================================================
// THERMODYNAMIC CALCULATION FUNCTIONS
// ============================================================================

/**
 * Calculate isentropic outlet temperature for compression/expansion
 * 
 * Based on the isentropic relation for ideal gas:
 * T2/T1 = (P2/P1)^((γ-1)/γ)
 * 
 * @param {number} T1 - Inlet temperature (K)
 * @param {number} pressureRatio - P2/P1
 * @param {number} gamma - Ratio of specific heats
 * @returns {number} Isentropic outlet temperature (K)
 */
const isentropicTemperature = (T1, pressureRatio, gamma) => {
    const exponent = (gamma - 1) / gamma;
    return T1 * Math.pow(pressureRatio, exponent);
};

/**
 * Calculate actual outlet temperature considering isentropic efficiency
 * 
 * For compression:
 * η_s = (T2s - T1) / (T2 - T1)
 * Therefore: T2 = T1 + (T2s - T1) / η_s
 * 
 * @param {number} T1 - Inlet temperature (K)
 * @param {number} T2s - Isentropic outlet temperature (K)
 * @param {number} efficiency - Isentropic efficiency (0-1)
 * @param {boolean} isCompression - true for compression, false for expansion
 * @returns {number} Actual outlet temperature (K)
 */
const actualTemperature = (T1, T2s, efficiency, isCompression = true) => {
    if (isCompression) {
        // Compressor: actual work > isentropic work
        return T1 + (T2s - T1) / efficiency;
    } else {
        // Turbine: actual work < isentropic work
        return T1 - efficiency * (T1 - T2s);
    }
};

/**
 * Calculate compressor/turbine work
 * 
 * From first law (steady flow energy equation):
 * Ẇ = ṁ * cp * (T2 - T1)
 * 
 * Sign convention: 
 * - Positive for work input (compression)
 * - Negative for work output (expansion)
 * 
 * @param {number} massFlow - Mass flow rate (kg/s)
 * @param {number} cp - Specific heat at constant pressure (J/kg·K)
 * @param {number} T1 - Inlet temperature (K)
 * @param {number} T2 - Outlet temperature (K)
 * @returns {number} Power (W)
 */
const shaftWork = (massFlow, cp, T1, T2) => {
    return massFlow * cp * (T2 - T1);
};

/**
 * Calculate polytropic exponent from temperatures and pressures
 * 
 * From polytropic relation: T2/T1 = (P2/P1)^((n-1)/n)
 * Solving for n: n = 1 / (1 - ln(T2/T1) / ln(P2/P1))
 * 
 * @param {number} T1 - Inlet temperature (K)
 * @param {number} T2 - Outlet temperature (K)
 * @param {number} P1 - Inlet pressure (Pa)
 * @param {number} P2 - Outlet pressure (Pa)
 * @returns {number} Polytropic exponent n
 */
const polytropicExponent = (T1, T2, P1, P2) => {
    const tempRatio = Math.log(T2 / T1);
    const pressRatio = Math.log(P2 / P1);
    
    if (Math.abs(pressRatio) < 1e-10) return 1; // Avoid division by zero
    
    return 1 / (1 - tempRatio / pressRatio);
};

/**
 * Calculate polytropic efficiency from isentropic efficiency
 * 
 * η_p = (γ-1)/γ * ln(r_p) / ln(1 + (r_p^((γ-1)/γ) - 1)/η_s)
 * 
 * @param {number} isentropicEfficiency - Isentropic efficiency (0-1)
 * @param {number} pressureRatio - Pressure ratio
 * @param {number} gamma - Ratio of specific heats
 * @returns {number} Polytropic efficiency (0-1)
 */
const polytropicEfficiency = (isentropicEfficiency, pressureRatio, gamma) => {
    const exp = (gamma - 1) / gamma;
    const term = Math.pow(pressureRatio, exp);
    const numerator = exp * Math.log(pressureRatio);
    const denominator = Math.log(1 + (term - 1) / isentropicEfficiency);
    
    return numerator / denominator;
};

/**
 * Calculate gas density from ideal gas law
 * 
 * ρ = P / (R * T)
 * 
 * @param {number} pressure - Pressure (Pa)
 * @param {number} R - Specific gas constant (J/kg·K)
 * @param {number} temperature - Temperature (K)
 * @returns {number} Density (kg/m³)
 */
const idealGasDensity = (pressure, R, temperature) => {
    return pressure / (R * temperature);
};

/**
 * Calculate specific volume from ideal gas law
 * 
 * v = R * T / P
 * 
 * @param {number} R - Specific gas constant (J/kg·K)
 * @param {number} temperature - Temperature (K)
 * @param {number} pressure - Pressure (Pa)
 * @returns {number} Specific volume (m³/kg)
 */
const specificVolume = (R, temperature, pressure) => {
    return R * temperature / pressure;
};

/**
 * Calculate entropy change for ideal gas
 * 
 * Δs = cp * ln(T2/T1) - R * ln(P2/P1)
 * 
 * @param {number} cp - Specific heat at constant pressure (J/kg·K)
 * @param {number} R - Specific gas constant (J/kg·K)
 * @param {number} T1 - Initial temperature (K)
 * @param {number} T2 - Final temperature (K)
 * @param {number} P1 - Initial pressure (Pa)
 * @param {number} P2 - Final pressure (Pa)
 * @returns {number} Entropy change (J/kg·K)
 */
const entropyChange = (cp, R, T1, T2, P1, P2) => {
    return cp * Math.log(T2 / T1) - R * Math.log(P2 / P1);
};

// ============================================================================
// STEAM TABLE APPROXIMATIONS
// ============================================================================
// TODO: Implement more accurate steam property calculations
// WHY NEEDED: For Rankine cycle simulations requiring accurate enthalpy/entropy
// HOW TO IMPLEMENT:
// 1. Use IAPWS-IF97 correlations for industrial-grade accuracy
// 2. Or implement polynomial fits to steam tables
// 3. Reference: Wagner, W.; Kretzschmar, H.-J. (2008) IAPWS-IF97

/**
 * Approximate saturation temperature from pressure (water/steam)
 * 
 * Simple Antoine equation approximation:
 * log10(P) = A - B/(C + T)
 * 
 * Valid for P: 0.01 - 100 bar
 * 
 * @param {number} pressure - Pressure (Pa)
 * @returns {number} Saturation temperature (K)
 */
const saturationTemperature = (pressure) => {
    // Convert to bar for correlation
    const P_bar = pressure / 1e5;
    
    // Antoine coefficients for water (approximate)
    const A = 5.11564;
    const B = 1687.537;
    const C = 230.17;
    
    // Solve for T (in °C, then convert to K)
    const T_C = B / (A - Math.log10(P_bar * 750.06)) - C;
    return T_C + 273.15;
};

/**
 * Approximate saturation pressure from temperature (water/steam)
 * 
 * @param {number} temperature - Temperature (K)
 * @returns {number} Saturation pressure (Pa)
 */
const saturationPressure = (temperature) => {
    const T_C = temperature - 273.15;
    
    // Antoine coefficients
    const A = 5.11564;
    const B = 1687.537;
    const C = 230.17;
    
    const log10P_mmHg = A - B / (C + T_C);
    const P_bar = Math.pow(10, log10P_mmHg) / 750.06;
    
    return P_bar * 1e5;
};

/**
 * Approximate steam enthalpy (simplified correlation)
 * 
 * TODO: Replace with full IAPWS-IF97 implementation for production use
 * 
 * @param {number} temperature - Temperature (K)
 * @param {number} pressure - Pressure (Pa)
 * @param {string} phase - 'liquid', 'vapor', or 'superheated'
 * @returns {number} Specific enthalpy (kJ/kg)
 */
const steamEnthalpy = (temperature, pressure, phase = 'vapor') => {
    const T_C = temperature - 273.15;
    const P_bar = pressure / 1e5;
    
    if (phase === 'liquid') {
        // Approximate: h ≈ 4.18 * T (for compressed liquid at low pressures)
        return 4.18 * T_C;
    } else if (phase === 'vapor' || phase === 'saturated') {
        // Approximate saturated vapor enthalpy
        // Linear fit valid for 1-50 bar
        return 2675 + 0.5 * T_C;
    } else {
        // Superheated steam: h ≈ h_sat + cp * (T - T_sat)
        const T_sat = saturationTemperature(pressure) - 273.15;
        const h_sat = 2675 + 0.5 * T_sat;
        const cp = 2.0; // kJ/kg·K approximate for superheated steam
        return h_sat + cp * (T_C - T_sat);
    }
};

// ============================================================================
// REACTION KINETICS HELPERS
// ============================================================================

/**
 * Arrhenius rate constant
 * 
 * k = A * exp(-Ea / (R * T))
 * 
 * @param {number} A - Pre-exponential factor
 * @param {number} Ea - Activation energy (J/mol)
 * @param {number} T - Temperature (K)
 * @returns {number} Rate constant (units depend on reaction order)
 */
const arrheniusRate = (A, Ea, T) => {
    const R = GAS_CONSTANTS.SI;
    return A * Math.exp(-Ea / (R * T));
};

/**
 * First-order reaction rate
 * 
 * r = k * C
 * 
 * @param {number} k - Rate constant (1/s)
 * @param {number} C - Concentration (mol/m³)
 * @returns {number} Reaction rate (mol/m³/s)
 */
const firstOrderRate = (k, C) => {
    return k * C;
};

/**
 * Second-order reaction rate
 * 
 * r = k * CA * CB (for A + B → products)
 * 
 * @param {number} k - Rate constant (m³/mol/s)
 * @param {number} CA - Concentration of A (mol/m³)
 * @param {number} CB - Concentration of B (mol/m³)
 * @returns {number} Reaction rate (mol/m³/s)
 */
const secondOrderRate = (k, CA, CB) => {
    return k * CA * CB;
};

/**
 * Calculate conversion from concentrations
 * 
 * X = (C0 - C) / C0
 * 
 * @param {number} C0 - Initial concentration
 * @param {number} C - Current concentration
 * @returns {number} Conversion (0-1)
 */
const conversion = (C0, C) => {
    return (C0 - C) / C0;
};

// ============================================================================
// MASS BALANCE HELPERS
// ============================================================================

/**
 * Simple steady-state mass balance
 * 
 * Accumulation = In - Out + Generation - Consumption
 * At steady state: In + Generation = Out + Consumption
 * 
 * @param {number} massIn - Mass flow in (kg/s)
 * @param {number} generation - Mass generation rate (kg/s)
 * @param {number} consumption - Mass consumption rate (kg/s)
 * @returns {number} Mass flow out (kg/s)
 */
const steadyStateMassBalance = (massIn, generation = 0, consumption = 0) => {
    return massIn + generation - consumption;
};

/**
 * Component mass balance with reaction
 * 
 * @param {number} flowIn - Molar flow in (mol/s)
 * @param {number} flowOut - Molar flow out (mol/s)
 * @param {number} stoichCoef - Stoichiometric coefficient (negative for reactants)
 * @param {number} reactionExtent - Extent of reaction (mol/s)
 * @returns {number} Accumulation rate (mol/s)
 */
const componentBalance = (flowIn, flowOut, stoichCoef, reactionExtent) => {
    return flowIn - flowOut + stoichCoef * reactionExtent;
};

// ============================================================================
// UNIT CONVERSIONS
// ============================================================================

const UnitConversions = {
    // Temperature
    celsiusToKelvin: (C) => C + 273.15,
    kelvinToCelsius: (K) => K - 273.15,
    fahrenheitToKelvin: (F) => (F - 32) * 5/9 + 273.15,
    kelvinToFahrenheit: (K) => (K - 273.15) * 9/5 + 32,
    
    // Pressure
    barToPascal: (bar) => bar * 1e5,
    pascalToBar: (Pa) => Pa / 1e5,
    atmToPascal: (atm) => atm * 101325,
    pascalToAtm: (Pa) => Pa / 101325,
    psiToPascal: (psi) => psi * 6894.76,
    pascalToPsi: (Pa) => Pa / 6894.76,
    kPaToPascal: (kPa) => kPa * 1000,
    pascalToKPa: (Pa) => Pa / 1000,
    MPaToPascal: (MPa) => MPa * 1e6,
    pascalToMPa: (Pa) => Pa / 1e6,
    
    // Energy
    kJToJ: (kJ) => kJ * 1000,
    JToKJ: (J) => J / 1000,
    kWhToJ: (kWh) => kWh * 3.6e6,
    JToKWh: (J) => J / 3.6e6,
    
    // Power
    kWToW: (kW) => kW * 1000,
    WToKW: (W) => W / 1000,
    hpToW: (hp) => hp * 745.7,
    WToHp: (W) => W / 745.7,
    
    // Mass flow
    kgPerHrToKgPerS: (kgph) => kgph / 3600,
    kgPerSToKgPerHr: (kgps) => kgps * 3600
};

// ============================================================================
// INTERPOLATION UTILITIES
// ============================================================================

/**
 * Linear interpolation
 * 
 * @param {number} x - Input value
 * @param {number} x0 - Lower bound x
 * @param {number} x1 - Upper bound x
 * @param {number} y0 - Value at x0
 * @param {number} y1 - Value at x1
 * @returns {number} Interpolated y value
 */
const lerp = (x, x0, x1, y0, y1) => {
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
};

/**
 * Find value in table using linear interpolation
 * 
 * @param {number} x - Input value
 * @param {Array} xArray - Array of x values (must be sorted ascending)
 * @param {Array} yArray - Array of corresponding y values
 * @returns {number} Interpolated y value
 */
const tableInterpolate = (x, xArray, yArray) => {
    if (x <= xArray[0]) return yArray[0];
    if (x >= xArray[xArray.length - 1]) return yArray[yArray.length - 1];
    
    for (let i = 0; i < xArray.length - 1; i++) {
        if (x >= xArray[i] && x <= xArray[i + 1]) {
            return lerp(x, xArray[i], xArray[i + 1], yArray[i], yArray[i + 1]);
        }
    }
    
    return yArray[yArray.length - 1];
};

// ============================================================================
// STATISTICS AND DATA PROCESSING
// ============================================================================

/**
 * Calculate mean of an array
 */
const mean = (arr) => arr.reduce((a, b) => a + b, 0) / arr.length;

/**
 * Calculate standard deviation
 */
const stdDev = (arr) => {
    const avg = mean(arr);
    const squareDiffs = arr.map(value => Math.pow(value - avg, 2));
    return Math.sqrt(mean(squareDiffs));
};

/**
 * Find minimum value in array
 */
const min = (arr) => Math.min(...arr);

/**
 * Find maximum value in array
 */
const max = (arr) => Math.max(...arr);

// ============================================================================
// NUMERICAL METHODS FOR EQUATIONS
// ============================================================================

/**
 * Newton-Raphson method for finding roots
 * 
 * Finds x such that f(x) = 0
 * 
 * @param {Function} f - Function to find root of
 * @param {Function} df - Derivative of f
 * @param {number} x0 - Initial guess
 * @param {number} tolerance - Convergence tolerance
 * @param {number} maxIter - Maximum iterations
 * @returns {Object} {root: number, iterations: number, converged: boolean}
 */
const newtonRaphson = (f, df, x0, tolerance = 1e-10, maxIter = 100) => {
    let x = x0;
    
    for (let i = 0; i < maxIter; i++) {
        const fx = f(x);
        const dfx = df(x);
        
        if (Math.abs(dfx) < 1e-15) {
            return { root: x, iterations: i, converged: false };
        }
        
        const xNew = x - fx / dfx;
        
        if (Math.abs(xNew - x) < tolerance) {
            return { root: xNew, iterations: i + 1, converged: true };
        }
        
        x = xNew;
    }
    
    return { root: x, iterations: maxIter, converged: false };
};

/**
 * Bisection method for finding roots
 * 
 * More robust than Newton-Raphson but slower
 * 
 * @param {Function} f - Function to find root of
 * @param {number} a - Lower bound
 * @param {number} b - Upper bound
 * @param {number} tolerance - Convergence tolerance
 * @param {number} maxIter - Maximum iterations
 * @returns {Object} {root: number, iterations: number, converged: boolean}
 */
const bisection = (f, a, b, tolerance = 1e-10, maxIter = 100) => {
    if (f(a) * f(b) > 0) {
        console.warn('Bisection: f(a) and f(b) must have opposite signs');
        return { root: NaN, iterations: 0, converged: false };
    }
    
    for (let i = 0; i < maxIter; i++) {
        const c = (a + b) / 2;
        const fc = f(c);
        
        if (Math.abs(fc) < tolerance || (b - a) / 2 < tolerance) {
            return { root: c, iterations: i + 1, converged: true };
        }
        
        if (f(a) * fc < 0) {
            b = c;
        } else {
            a = c;
        }
    }
    
    return { root: (a + b) / 2, iterations: maxIter, converged: false };
};

// ============================================================================
// EXPORT ALL UTILITIES
// ============================================================================

const MathUtils = {
    // Integrators
    euler,
    eulerSystem,
    rk4,
    rk4System,
    
    // Constants
    GAS_CONSTANTS,
    GAS_PROPERTIES,
    
    // Thermodynamics
    isentropicTemperature,
    actualTemperature,
    shaftWork,
    polytropicExponent,
    polytropicEfficiency,
    idealGasDensity,
    specificVolume,
    entropyChange,
    
    // Steam properties
    saturationTemperature,
    saturationPressure,
    steamEnthalpy,
    
    // Kinetics
    arrheniusRate,
    firstOrderRate,
    secondOrderRate,
    conversion,
    
    // Mass balance
    steadyStateMassBalance,
    componentBalance,
    
    // Units
    UnitConversions,
    
    // Interpolation
    lerp,
    tableInterpolate,
    
    // Statistics
    mean,
    stdDev,
    min,
    max,
    
    // Root finding
    newtonRaphson,
    bisection
};

// Export for ES6 modules (if using module bundler)
// export default MathUtils;

// Make available globally for browser usage
if (typeof window !== 'undefined') {
    window.MathUtils = MathUtils;
}
