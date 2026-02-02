/**
 * ============================================================================
 * DISTILLATION COLUMN SIMULATION - ADVANCED
 * ============================================================================
 * 
 * Complete distillation column design with:
 *   - McCabe-Thiele graphical method
 *   - Non-ideal VLE (Activity Coefficients: Wilson, Margules, Van Laar)
 *   - Ponchon-Savarit method (Enthalpy-Concentration)
 *   - Multi-component FUG (Fenske-Underwood-Gilliland)
 *   - Hydraulic calculations (Pressure drop, flooding, column sizing)
 * 
 * @version 2.0.0-advanced
 */

'use strict';

const DistillationModule = {
    name: 'distillation',
    
    defaults: {
        // Basic parameters
        feedComposition: 0.5,
        distillateComposition: 0.95,
        bottomsComposition: 0.05,
        feedCondition: 1.0,
        relativeVolatility: 2.5,
        refluxRatio: 2.0,
        columnType: 'tray',
        hetp: 0.5,
        feedFlowRate: 100,
        trayEfficiency: 0.7,
        
        // VLE Model
        vleModel: 'ideal',  // 'ideal', 'margules', 'vanlaar', 'wilson', 'nrtl'
        
        // Activity coefficient parameters (for binary)
        A12: 0.5,  // Margules/Van Laar parameter
        A21: 0.7,
        
        // Enthalpy parameters (Ponchon-Savarit)
        hL_A: 0,      // Liquid enthalpy pure A (kJ/mol)
        hL_B: 0,
        hV_A: 35,     // Vapor enthalpy pure A (kJ/mol) 
        hV_B: 40,
        hMix: 2,      // Mixing enthalpy parameter
        
        // Hydraulics
        pressure: 101.325,    // kPa
        vaporDensity: 2.5,    // kg/m³
        liquidDensity: 800,   // kg/m³
        surfaceTension: 0.02, // N/m
        traySpacing: 0.6,     // m
        weirHeight: 0.05,     // m
        holeArea: 0.1,        // fraction of active area
        
        // Multi-component (FUG)
        multiComponent: false,
        components: ['A', 'B', 'C'],
        feedMoleFractions: [0.3, 0.4, 0.3],
        alphas: [3.0, 1.5, 1.0],  // relative to heavy key
        lightKey: 0,   // index
        heavyKey: 2,   // index
        lkRecovery: 0.99,
        hkRecovery: 0.01
    },
    
    params: null,
    results: null,
    charts: {}
};

// ============================================================================
// VLE CALCULATIONS - IDEAL & NON-IDEAL
// ============================================================================

/**
 * Calculate activity coefficients using various models
 */
const ActivityModels = {
    /**
     * Margules 2-suffix equation
     * ln γ1 = x2² [A12 + 2(A21 - A12)x1]
     * ln γ2 = x1² [A21 + 2(A12 - A21)x2]
     */
    margules: function(x1, A12, A21) {
        const x2 = 1 - x1;
        const lnGamma1 = x2 * x2 * (A12 + 2 * (A21 - A12) * x1);
        const lnGamma2 = x1 * x1 * (A21 + 2 * (A12 - A21) * x2);
        return {
            gamma1: Math.exp(lnGamma1),
            gamma2: Math.exp(lnGamma2)
        };
    },
    
    /**
     * Van Laar equation
     * ln γ1 = A12 / (1 + A12*x1/(A21*x2))²
     * ln γ2 = A21 / (1 + A21*x2/(A12*x1))²
     */
    vanlaar: function(x1, A12, A21) {
        const x2 = 1 - x1;
        if (x1 < 0.001 || x2 < 0.001) {
            return { gamma1: Math.exp(A12), gamma2: Math.exp(A21) };
        }
        const term1 = 1 + (A12 * x1) / (A21 * x2);
        const term2 = 1 + (A21 * x2) / (A12 * x1);
        const lnGamma1 = A12 / (term1 * term1);
        const lnGamma2 = A21 / (term2 * term2);
        return {
            gamma1: Math.exp(lnGamma1),
            gamma2: Math.exp(lnGamma2)
        };
    },
    
    /**
     * Wilson equation (simplified for binary)
     * ln γ1 = -ln(x1 + Λ12*x2) + x2*(Λ12/(x1+Λ12*x2) - Λ21/(x2+Λ21*x1))
     */
    wilson: function(x1, Lambda12, Lambda21) {
        const x2 = 1 - x1;
        const sum1 = x1 + Lambda12 * x2;
        const sum2 = x2 + Lambda21 * x1;
        
        const lnGamma1 = -Math.log(sum1) + x2 * (Lambda12/sum1 - Lambda21/sum2);
        const lnGamma2 = -Math.log(sum2) - x1 * (Lambda12/sum1 - Lambda21/sum2);
        return {
            gamma1: Math.exp(lnGamma1),
            gamma2: Math.exp(lnGamma2)
        };
    },
    
    /**
     * NRTL equation (simplified, α = 0.3 typical)
     */
    nrtl: function(x1, tau12, tau21, alpha = 0.3) {
        const x2 = 1 - x1;
        const G12 = Math.exp(-alpha * tau12);
        const G21 = Math.exp(-alpha * tau21);
        
        const sum1 = x1 + x2 * G21;
        const sum2 = x2 + x1 * G12;
        
        const lnGamma1 = x2 * x2 * (tau21 * Math.pow(G21/sum1, 2) + 
                         tau12 * G12 / (sum2 * sum2));
        const lnGamma2 = x1 * x1 * (tau12 * Math.pow(G12/sum2, 2) + 
                         tau21 * G21 / (sum1 * sum1));
        return {
            gamma1: Math.exp(lnGamma1),
            gamma2: Math.exp(lnGamma2)
        };
    }
};

/**
 * Calculate equilibrium vapor composition (modified Raoult's law)
 * y1 = γ1 * x1 * P1sat / P
 * For relative volatility form with activity coefficients:
 * y = (γ1/γ2) * α * x / (1 + ((γ1/γ2)*α - 1) * x)
 */
function equilibriumY(x, alpha, params) {
    if (params.vleModel === 'ideal') {
        return (alpha * x) / (1 + (alpha - 1) * x);
    }
    
    // Non-ideal: get activity coefficients
    let gammas;
    switch (params.vleModel) {
        case 'margules':
            gammas = ActivityModels.margules(x, params.A12, params.A21);
            break;
        case 'vanlaar':
            gammas = ActivityModels.vanlaar(x, params.A12, params.A21);
            break;
        case 'wilson':
            gammas = ActivityModels.wilson(x, params.A12, params.A21);
            break;
        case 'nrtl':
            gammas = ActivityModels.nrtl(x, params.A12, params.A21);
            break;
        default:
            gammas = { gamma1: 1, gamma2: 1 };
    }
    
    // Modified relative volatility
    const alphaEff = alpha * gammas.gamma1 / gammas.gamma2;
    return (alphaEff * x) / (1 + (alphaEff - 1) * x);
}

/**
 * Generate equilibrium curve with activity coefficients
 */
function generateEquilibriumCurve(params) {
    const points = [];
    const alpha = params.relativeVolatility;
    
    for (let i = 0; i <= 100; i++) {
        const x = i / 100;
        const y = equilibriumY(x, alpha, params);
        
        // Also calculate activity coefficients for display
        let gamma1 = 1, gamma2 = 1;
        if (params.vleModel !== 'ideal') {
            const gammas = ActivityModels[params.vleModel](x, params.A12, params.A21);
            gamma1 = gammas.gamma1;
            gamma2 = gammas.gamma2;
        }
        
        points.push({ x, y, gamma1, gamma2 });
    }
    return points;
}

/**
 * Find x given y on equilibrium curve (inverse, for non-ideal requires iteration)
 */
function findEquilibriumX(yTarget, alpha, params) {
    // Bisection method
    let xLow = 0.001, xHigh = 0.999;
    for (let i = 0; i < 50; i++) {
        const xMid = (xLow + xHigh) / 2;
        const yMid = equilibriumY(xMid, alpha, params);
        if (yMid < yTarget) {
            xLow = xMid;
        } else {
            xHigh = xMid;
        }
        if (Math.abs(yMid - yTarget) < 1e-6) break;
    }
    return (xLow + xHigh) / 2;
}

// ============================================================================
// PONCHON-SAVARIT METHOD (Enthalpy-Concentration)
// ============================================================================

/**
 * Calculate liquid enthalpy as function of composition
 * hL = x*hL_A + (1-x)*hL_B + x*(1-x)*hMix
 */
function liquidEnthalpy(x, params) {
    const { hL_A, hL_B, hMix } = params;
    return x * hL_A + (1 - x) * hL_B + x * (1 - x) * hMix;
}

/**
 * Calculate vapor enthalpy as function of composition
 * hV = y*hV_A + (1-y)*hV_B
 */
function vaporEnthalpy(y, params) {
    const { hV_A, hV_B } = params;
    return y * hV_A + (1 - y) * hV_B;
}

/**
 * Generate H-x-y diagram data for Ponchon-Savarit
 */
function generateEnthalpyDiagram(params) {
    const liquidLine = [];
    const vaporLine = [];
    const tieLines = [];
    
    for (let i = 0; i <= 100; i += 5) {
        const x = i / 100;
        const y = equilibriumY(x, params.relativeVolatility, params);
        
        const hL = liquidEnthalpy(x, params);
        const hV = vaporEnthalpy(y, params);
        
        liquidLine.push({ x, h: hL });
        vaporLine.push({ x: y, h: hV });
        
        // Tie line connects (x, hL) to (y, hV)
        if (i % 10 === 0) {
            tieLines.push({ x, hL, y, hV });
        }
    }
    
    return { liquidLine, vaporLine, tieLines };
}

/**
 * Ponchon-Savarit stage calculation
 * Uses enthalpy balance and operating poles
 */
function ponchonSavarit(params) {
    const {
        feedComposition, distillateComposition, bottomsComposition,
        feedCondition, relativeVolatility, refluxRatio, feedFlowRate
    } = params;
    
    const zF = feedComposition;
    const xD = distillateComposition;
    const xB = bottomsComposition;
    const R = refluxRatio;
    
    // Material balance
    const F = feedFlowRate;
    const D = F * (zF - xB) / (xD - xB);
    const B = F - D;
    
    // Condenser duty: QC = D * (R + 1) * (hV_D - hL_D)
    const hV_D = vaporEnthalpy(xD, params);
    const hL_D = liquidEnthalpy(xD, params);
    const QC = D * (R + 1) * (hV_D - hL_D);
    
    // Rectifying pole: (xD, hD + QC/D)
    const rectPoleH = hL_D + QC / D;
    
    // Feed enthalpy (based on q)
    const hF = feedCondition * liquidEnthalpy(zF, params) + 
               (1 - feedCondition) * vaporEnthalpy(zF, params);
    
    // Reboiler duty: QB = QC + F * (hF - D*hL_D - B*hL_B) / B approx
    const hL_B = liquidEnthalpy(xB, params);
    const QB = QC + F * (hF - hL_D * D / F - hL_B * B / F);
    
    // Stripping pole: (xB, hB - QB/B)
    const stripPoleH = hL_B - Math.abs(QB) / B;
    
    // Step off stages (simplified - full implementation would trace lines)
    const stages = [];
    let x = xD;
    let n = 0;
    let feedStage = 0;
    
    while (x > xB && n < 50) {
        n++;
        const y = equilibriumY(x, relativeVolatility, params);
        
        // Determine pole based on position
        const isRectifying = x > zF;
        
        stages.push({
            n, x, y,
            hL: liquidEnthalpy(x, params),
            hV: vaporEnthalpy(y, params),
            section: isRectifying ? 'rectifying' : 'stripping'
        });
        
        if (!isRectifying && feedStage === 0) feedStage = n;
        
        // Next stage liquid (simplified linear interpolation)
        x = x - (xD - xB) / 20;  // Approximate
    }
    
    return {
        theoreticalStages: n,
        feedStage,
        QC, QB,
        poles: {
            rectifying: { x: xD, h: rectPoleH },
            stripping: { x: xB, h: stripPoleH }
        },
        stages,
        enthalpyDiagram: generateEnthalpyDiagram(params),
        materialBalance: { F, D, B }
    };
}

// ============================================================================
// MULTI-COMPONENT DISTILLATION (FUG Method)
// ============================================================================

/**
 * Fenske equation for minimum stages at total reflux
 * Nmin = ln[(xLK/xHK)_D * (xHK/xLK)_B] / ln(αLK/αHK)
 */
function fenskeMinStages(params) {
    const { components, alphas, lightKey, heavyKey, lkRecovery, hkRecovery, feedMoleFractions } = params;
    
    const zLK = feedMoleFractions[lightKey];
    const zHK = feedMoleFractions[heavyKey];
    const alphaLK = alphas[lightKey];
    const alphaHK = alphas[heavyKey];
    
    // Recoveries to compositions
    // Assume distillate primarily LK, bottoms primarily HK
    const xLK_D = lkRecovery * zLK / (lkRecovery * zLK + (1 - hkRecovery) * zHK);
    const xHK_D = 1 - xLK_D;
    const xLK_B = (1 - lkRecovery) * zLK / ((1 - lkRecovery) * zLK + hkRecovery * zHK);
    const xHK_B = 1 - xLK_B;
    
    // Fenske equation
    const Nmin = Math.log((xLK_D / xHK_D) * (xHK_B / xLK_B)) / Math.log(alphaLK / alphaHK);
    
    return {
        Nmin: Math.max(1, Nmin),
        xLK_D, xHK_D, xLK_B, xHK_B
    };
}

/**
 * Underwood equations for minimum reflux
 * Sum(αi * zi / (αi - θ)) = 1 - q  (find θ between αLK and αHK)
 * Rmin + 1 = Sum(αi * xD_i / (αi - θ))
 */
function underwoodMinReflux(params) {
    const { components, alphas, feedMoleFractions, feedCondition, lightKey, heavyKey, lkRecovery, hkRecovery } = params;
    
    const q = feedCondition;
    const alphaLK = alphas[lightKey];
    const alphaHK = alphas[heavyKey];
    
    // Find theta by bisection (between αHK and αLK)
    let thetaLow = alphaHK + 0.001;
    let thetaHigh = alphaLK - 0.001;
    let theta = (thetaLow + thetaHigh) / 2;
    
    for (let iter = 0; iter < 50; iter++) {
        theta = (thetaLow + thetaHigh) / 2;
        
        let sum = 0;
        for (let i = 0; i < components.length; i++) {
            sum += alphas[i] * feedMoleFractions[i] / (alphas[i] - theta);
        }
        
        if (sum > 1 - q) {
            thetaLow = theta;
        } else {
            thetaHigh = theta;
        }
        
        if (Math.abs(sum - (1 - q)) < 1e-6) break;
    }
    
    // Estimate distillate compositions based on recoveries
    const xD = feedMoleFractions.map((z, i) => {
        if (i === lightKey) return lkRecovery * z;
        if (i === heavyKey) return (1 - hkRecovery) * z;
        // Distribute others based on relative volatility
        return alphas[i] > 1 ? 0.9 * z : 0.1 * z;
    });
    const sumXD = xD.reduce((a, b) => a + b, 0);
    const xD_norm = xD.map(x => x / sumXD);
    
    // Calculate Rmin + 1
    let RminPlus1 = 0;
    for (let i = 0; i < components.length; i++) {
        RminPlus1 += alphas[i] * xD_norm[i] / (alphas[i] - theta);
    }
    
    return {
        theta,
        Rmin: Math.max(0.1, RminPlus1 - 1),
        xD: xD_norm
    };
}

/**
 * Gilliland correlation for actual stages
 * (N - Nmin)/(N + 1) = f((R - Rmin)/(R + 1))
 */
function gillilandStages(Nmin, Rmin, R) {
    const X = (R - Rmin) / (R + 1);
    
    // Gilliland correlation (Molokanov form)
    const Y = 1 - Math.exp((1 + 54.4 * X) / (11 + 117.2 * X) * (X - 1) / Math.sqrt(X));
    
    // N from Y = (N - Nmin)/(N + 1)
    const N = (Nmin + Y) / (1 - Y);
    
    return Math.ceil(N);
}

/**
 * Complete FUG method calculation
 */
function fugMethod(params) {
    const fenske = fenskeMinStages(params);
    const underwood = underwoodMinReflux(params);
    const R = params.refluxRatio;
    
    if (R < underwood.Rmin) {
        return {
            error: true,
            message: `Reflux (${R.toFixed(2)}) < Rmin (${underwood.Rmin.toFixed(2)})`
        };
    }
    
    const N = gillilandStages(fenske.Nmin, underwood.Rmin, R);
    
    // Feed stage (Kirkbride correlation)
    const { feedMoleFractions, lightKey, heavyKey } = params;
    const zLK = feedMoleFractions[lightKey];
    const zHK = feedMoleFractions[heavyKey];
    const B_D = (1 - zLK) / zLK;  // Approximate
    
    const logRatio = 0.206 * Math.log10((zHK / zLK) * Math.pow(fenske.xLK_B / fenske.xHK_D, 2) * B_D);
    const feedStage = Math.round(N / (1 + Math.pow(10, logRatio)));
    
    return {
        Nmin: fenske.Nmin,
        Rmin: underwood.Rmin,
        theta: underwood.theta,
        theoreticalStages: N,
        feedStage,
        RoverRmin: R / underwood.Rmin,
        distillateCompositions: underwood.xD,
        components: params.components
    };
}

// ============================================================================
// HYDRAULIC CALCULATIONS
// ============================================================================

/**
 * Calculate column diameter based on flooding velocity
 * Using Fair correlation for sieve trays
 */
function calculateColumnDiameter(params, vaporFlow, liquidFlow) {
    const { vaporDensity, liquidDensity, surfaceTension, traySpacing } = params;
    
    const rhoV = vaporDensity;
    const rhoL = liquidDensity;
    const sigma = surfaceTension;
    const ts = traySpacing;
    
    // Flow parameter
    const Flv = (liquidFlow / vaporFlow) * Math.sqrt(rhoV / rhoL);
    
    // Capacity parameter from Fair correlation (simplified)
    // Csb_flood = f(Flv, ts)
    const Csb_flood = 0.1 * Math.pow(ts, 0.5) * Math.exp(-1.463 * Math.pow(Flv, 0.842));
    
    // Surface tension correction
    const Csb = Csb_flood * Math.pow(sigma / 0.02, 0.2);
    
    // Flooding velocity
    const Uf = Csb * Math.sqrt((rhoL - rhoV) / rhoV);
    
    // Design velocity (80% of flooding)
    const U = 0.8 * Uf;
    
    // Column area
    const Ac = vaporFlow / (rhoV * U);
    
    // Column diameter
    const Dc = Math.sqrt(4 * Ac / Math.PI);
    
    return {
        floodingVelocity: Uf,
        designVelocity: U,
        columnArea: Ac,
        columnDiameter: Dc,
        percentFlood: 80,
        flowParameter: Flv
    };
}

/**
 * Calculate tray pressure drop
 * ΔP_tray = ΔP_dry + ΔP_static + ΔP_residual
 */
function calculatePressureDrop(params, vaporVelocity) {
    const { vaporDensity, liquidDensity, weirHeight, holeArea, traySpacing } = params;
    
    const rhoV = vaporDensity;
    const rhoL = liquidDensity;
    const hw = weirHeight;
    const Ah_At = holeArea;
    
    // Hole velocity
    const Uh = vaporVelocity / Ah_At;
    
    // Dry tray pressure drop (orifice equation)
    // ΔP_dry = K * (rhoV * Uh²) / 2
    const Co = 0.85;  // Orifice coefficient
    const dP_dry = (1 / (Co * Co) - 1) * rhoV * Uh * Uh / 2;
    
    // Static head (liquid on tray)
    // how = weir crest, assume 25mm typical
    const how = 0.025;
    const dP_static = rhoL * 9.81 * (hw + how);
    
    // Residual (bubble formation) ~12.5mm liquid
    const dP_residual = rhoL * 9.81 * 0.0125;
    
    // Total per tray
    const dP_tray = dP_dry + dP_static + dP_residual;
    
    // Convert to mmHg for display
    const dP_mmHg = dP_tray / 133.322;
    
    return {
        dP_dry: dP_dry / 1000,      // kPa
        dP_static: dP_static / 1000,
        dP_residual: dP_residual / 1000,
        dP_tray: dP_tray / 1000,
        dP_mmHg_per_tray: dP_mmHg
    };
}

/**
 * Check for weeping (minimum vapor velocity)
 */
function checkWeeping(params, vaporVelocity) {
    const { vaporDensity, liquidDensity, holeArea, weirHeight } = params;
    
    const Ah_At = holeArea;
    const Uh = vaporVelocity / Ah_At;
    
    // Minimum hole velocity to prevent weeping (simplified)
    const Uh_min = 0.5 * Math.sqrt(9.81 * weirHeight * liquidDensity / vaporDensity);
    
    return {
        holeVelocity: Uh,
        minHoleVelocity: Uh_min,
        isWeeping: Uh < Uh_min,
        weepingMargin: (Uh - Uh_min) / Uh_min * 100
    };
}

/**
 * Complete hydraulic analysis
 */
function hydraulicAnalysis(params, vaporFlow_kmol_h, liquidFlow_kmol_h) {
    // Convert molar to volumetric (assume MW ~100 g/mol average)
    const MW = 100;
    const vaporFlow = vaporFlow_kmol_h * MW / 3600 / params.vaporDensity;  // m³/s
    const liquidFlow = liquidFlow_kmol_h * MW / 3600 / params.liquidDensity;
    
    const sizing = calculateColumnDiameter(params, vaporFlow * params.vaporDensity, liquidFlow * params.liquidDensity);
    const pressureDrop = calculatePressureDrop(params, sizing.designVelocity);
    const weeping = checkWeeping(params, sizing.designVelocity);
    
    return {
        ...sizing,
        pressureDrop,
        weeping,
        vaporFlowVolumetric: vaporFlow,
        liquidFlowVolumetric: liquidFlow
    };
}

// ============================================================================
// MCCABE-THIELE METHOD
// ============================================================================

/**
 * Calculate minimum reflux ratio at pinch point
 */
function calculateMinReflux(params) {
    const { feedComposition, distillateComposition, feedCondition, relativeVolatility } = params;
    
    const zF = feedComposition;
    const xD = distillateComposition;
    const q = feedCondition;
    const alpha = relativeVolatility;
    
    // At x = zF (for q=1)
    const yEq = equilibriumY(zF, alpha, params);
    
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
        // For non-ideal, need to solve iteratively
        let xEq;
        if (params.vleModel === 'ideal') {
            // Direct solution for ideal case
            xEq = y / (alpha - y * (alpha - 1));
        } else {
            // Bisection for non-ideal
            xEq = findEquilibriumX(y, alpha, params);
        }
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
    // Check if multi-component mode
    if (params.multiComponent) {
        return calculateMultiComponent(params);
    }
    
    // Material balance
    const F = params.feedFlowRate;
    const zF = params.feedComposition;
    const xD = params.distillateComposition;
    const xB = params.bottomsComposition;
    
    // D and B from material balance
    const D = F * (zF - xB) / (xD - xB);
    const B = F - D;
    
    // Select method
    let mcResult;
    if (params.method === 'ponchon') {
        mcResult = ponchonSavarit(params);
    } else {
        mcResult = params.columnType === 'packed' 
            ? calculatePackedHeight(params) 
            : mcCabeThiele(params);
    }
    
    if (mcResult.error) {
        return mcResult;
    }
    
    // Actual stages with tray efficiency
    const actualStages = params.columnType === 'tray'
        ? Math.ceil(mcResult.theoreticalStages / params.trayEfficiency)
        : mcResult.theoreticalStages;
    
    // Hydraulic calculations
    const R = params.refluxRatio;
    const L = R * D;           // Liquid flow in rectifying section
    const V = L + D;           // Vapor flow
    const hydraulics = hydraulicAnalysis(params, V, L);
    
    // Total column pressure drop
    const totalPressureDrop = hydraulics.pressureDrop.dP_tray * actualStages;
    
    const calculationSteps = [
        {
            description: 'Bilan matière global',
            formula: `F = D + B → ${F.toFixed(1)} = D + B`,
            result: `D = ${D.toFixed(2)} kmol/h, B = ${B.toFixed(2)} kmol/h`
        },
        {
            description: 'Modèle VLE',
            formula: params.vleModel === 'ideal' ? 'Loi de Raoult idéale' : `Coefficients d'activité (${params.vleModel})`,
            result: params.vleModel === 'ideal' ? 'γ₁ = γ₂ = 1' : `A₁₂ = ${params.A12}, A₂₁ = ${params.A21}`
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
            formula: params.method === 'ponchon' ? 'Méthode de Ponchon-Savarit' : 'Méthode de McCabe-Thiele',
            result: `N = ${mcResult.theoreticalStages} plateaux`
        },
        {
            description: 'Étage d\'alimentation optimal',
            formula: 'Intersection des droites opératoires',
            result: `Étage n°${mcResult.feedStage} (depuis le haut)`
        },
        {
            description: 'Diamètre de colonne',
            formula: `Corrélation de Fair (80% engorgement)`,
            result: `D = ${hydraulics.columnDiameter.toFixed(2)} m`
        },
        {
            description: 'Perte de charge totale',
            formula: `ΔP = ${hydraulics.pressureDrop.dP_tray.toFixed(3)} kPa/plateau × ${actualStages}`,
            result: `ΔP total = ${totalPressureDrop.toFixed(2)} kPa`
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
        hydraulics,
        totalPressureDrop,
        calculationSteps,
        equilibriumCurve: generateEquilibriumCurve(params),
        timestamp: new Date().toISOString()
    };
}

/**
 * Multi-component distillation using FUG method
 */
function calculateMultiComponent(params) {
    const fugResult = fugMethod(params);
    
    if (fugResult.error) return fugResult;
    
    const F = params.feedFlowRate;
    const D = F * 0.5;  // Approximate split
    const B = F - D;
    
    const R = params.refluxRatio;
    const L = R * D;
    const V = L + D;
    const hydraulics = hydraulicAnalysis(params, V, L);
    
    const actualStages = Math.ceil(fugResult.theoreticalStages / params.trayEfficiency);
    
    const calculationSteps = [
        {
            description: 'Méthode FUG (Multi-composants)',
            formula: 'Fenske-Underwood-Gilliland',
            result: `${params.components.length} composants`
        },
        {
            description: 'Nmin (Fenske)',
            formula: 'À reflux total',
            result: `Nmin = ${fugResult.Nmin.toFixed(2)}`
        },
        {
            description: 'Rmin (Underwood)',
            formula: `θ = ${fugResult.theta.toFixed(3)}`,
            result: `Rmin = ${fugResult.Rmin.toFixed(3)}`
        },
        {
            description: 'N théoriques (Gilliland)',
            formula: `R/Rmin = ${fugResult.RoverRmin.toFixed(2)}`,
            result: `N = ${fugResult.theoreticalStages} plateaux`
        },
        {
            description: 'Étage alimentation (Kirkbride)',
            formula: 'Position optimale',
            result: `Étage n°${fugResult.feedStage}`
        },
        {
            description: 'Diamètre colonne',
            formula: 'Corrélation de Fair',
            result: `D = ${hydraulics.columnDiameter.toFixed(2)} m`
        }
    ];
    
    return {
        ...fugResult,
        actualStages,
        hydraulics,
        calculationSteps,
        materialBalance: { F, D, B }
    };
}

// ============================================================================
// UI RENDERING - ADVANCED
// ============================================================================

DistillationModule.render = function() {
    this.params = { ...this.defaults };
    
    return `
        <style>
            .tabs-container { border-bottom: 1px solid #ddd; margin-bottom: 15px; display: flex; gap: 10px; }
            .tab-btn { padding: 10px 20px; border: none; background: transparent; cursor: pointer; 
                       border-bottom: 2px solid transparent; transition: all 0.2s; font-weight: 500; }
            .tab-btn:hover { background-color: #f0f0f0; }
            .tab-btn.active { border-bottom-color: var(--primary); color: var(--primary); font-weight: bold; background-color: #e8f0fe; }
            .tab-panel { display: none; }
            .tab-panel.active { display: block; }
            .param-section { background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 15px; border: 1px solid #e9ecef; }
            .param-section h4 { margin: 0 0 10px 0; font-size: 0.9em; color: #495057; text-transform: uppercase; font-weight: 600; border-bottom: 1px solid #dee2e6; padding-bottom: 5px; }
            .hydraulics-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 10px; }
            .hydraulic-item { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                              color: white; padding: 12px; border-radius: 8px; text-align: center; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
            .hydraulic-item .value { font-size: 1.4em; font-weight: bold; }
            .hydraulic-item .label { font-size: 0.75em; opacity: 0.9; }
            .vle-curve { stroke-width: 2; fill: none; }
            .operating-line { stroke-width: 1.5; stroke-dasharray: 5,3; }
        </style>
        
        <div class="distillation-module" id="distillation-content" style="padding: 10px;">
            <!-- Header -->
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;">
                <div>
                    <h2 style="margin: 0; color: #2c3e50;">Simulation de Distillation Avancée v2.0</h2>
                    <p class="text-muted" style="margin: 5px 0 0; font-size: 0.9em;">McCabe-Thiele • Ponchon-Savarit • FUG • Hydraulique</p>
                </div>
                <div class="status-badge success" id="simulation-status" style="display: flex; align-items: center; gap: 5px; padding: 5px 10px; background-color: #d4edda; color: #155724; border-radius: 20px; font-size: 0.85em;">
                    <span class="status-dot" style="height: 8px; width: 8px; background-color: #28a745; border-radius: 50%; display: inline-block;"></span>
                    <span>Modèle Complet</span>
                </div>
            </div>
            
            <!-- Method Selection Tabs -->
            <div class="tabs-container">
                <button class="tab-btn active" data-tab="binary">Binaire (McCabe)</button>
                <button class="tab-btn" data-tab="ponchon">Ponchon-Savarit</button>
                <button class="tab-btn" data-tab="multicomp">Multi-composants (FUG)</button>
                <button class="tab-btn" data-tab="hydraulics">Hydraulique</button>
            </div>

            <div style="display: grid; grid-template-columns: 1fr 1.2fr; gap: 20px;">
                <!-- LEFT: Controls -->
                <div>
                    <!-- BINARY TAB -->
                    <div class="tab-panel active" id="panel-binary">
                        <div class="card mb-4" style="border: 1px solid #ddd; border-radius: 8px; overflow: hidden; margin-bottom: 20px;">
                            <div class="card-header" style="background-color: #f8f9fa; padding: 10px 15px; border-bottom: 1px solid #ddd; font-weight: bold;">
                                <h3 class="card-title" style="margin: 0; font-size: 1.1em;">Configuration</h3>
                            </div>
                            <div class="card-body" style="padding: 15px;">
                                <div class="mode-selector mb-4" style="display: flex; gap: 10px; margin-bottom: 15px;">
                                    <button class="mode-btn active" data-column="tray" style="flex: 1; cursor: pointer;">Plateaux</button>
                                    <button class="mode-btn" data-column="packed" style="flex: 1; cursor: pointer;">Garnie</button>
                                </div>
                                
                                <div class="param-section">
                                    <h4>Modèle VLE</h4>
                                    <select class="form-select" id="vle-model" style="width: 100%; padding: 8px; margin-bottom: 10px; border-radius: 4px; border: 1px solid #ccc;">
                                        <option value="ideal">Idéal (Raoult)</option>
                                        <option value="margules">Margules</option>
                                        <option value="vanlaar">Van Laar</option>
                                        <option value="wilson">Wilson</option>
                                        <option value="nrtl">NRTL</option>
                                    </select>
                                    <div id="activity-params" style="display:none; margin-top:10px;">
                                        <div class="param-grid" style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px;">
                                            <div class="form-group">
                                                <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">A₁₂</label>
                                                <input type="number" class="form-input" id="a12" value="0.5" step="0.1" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                            </div>
                                            <div class="form-group">
                                                <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">A₂₁</label>
                                                <input type="number" class="form-input" id="a21" value="0.7" step="0.1" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                            </div>
                                        </div>
                                    </div>
                                </div>
                                
                                <div class="param-section">
                                    <h4>Compositions</h4>
                                    <div class="form-group" style="margin-bottom: 10px;">
                                        <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">zF (alimentation)</label>
                                        <div style="display: flex; align-items: center; gap: 10px;">
                                            <input type="range" class="form-range" id="feed-comp" value="50" min="10" max="90" style="flex: 1;">
                                            <span class="range-value" style="width: 40px; text-align: right; font-weight: bold;">0.50</span>
                                        </div>
                                    </div>
                                    <div class="form-group" style="margin-bottom: 10px;">
                                        <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">xD (distillat)</label>
                                        <div style="display: flex; align-items: center; gap: 10px;">
                                            <input type="range" class="form-range" id="distillate-comp" value="95" min="80" max="99" style="flex: 1;">
                                            <span class="range-value" style="width: 40px; text-align: right; font-weight: bold;">0.95</span>
                                        </div>
                                    </div>
                                    <div class="form-group" style="margin-bottom: 10px;">
                                        <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">xB (résidu)</label>
                                        <div style="display: flex; align-items: center; gap: 10px;">
                                            <input type="range" class="form-range" id="bottoms-comp" value="5" min="1" max="20" style="flex: 1;">
                                            <span class="range-value" style="width: 40px; text-align: right; font-weight: bold;">0.05</span>
                                        </div>
                                    </div>
                                </div>
                                
                                <div class="param-section">
                                    <h4>Opération</h4>
                                    <div class="param-grid" style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px;">
                                        <div class="form-group">
                                            <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">α (volatilité)</label>
                                            <input type="number" class="form-input" id="alpha" value="2.5" min="1.1" step="0.1" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                        </div>
                                        <div class="form-group">
                                            <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">R (reflux)</label>
                                            <input type="number" class="form-input" id="reflux" value="2.0" min="0.5" step="0.1" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                        </div>
                                        <div class="form-group">
                                            <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">F (kmol/h)</label>
                                            <input type="number" class="form-input" id="feed-flow" value="100" min="10" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                        </div>
                                        <div class="form-group">
                                            <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">η (efficacité)</label>
                                            <input type="number" class="form-input" id="efficiency" value="0.7" min="0.3" max="1" step="0.05" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                    
                    <!-- PONCHON TAB -->
                    <div class="tab-panel" id="panel-ponchon">
                         <div class="card mb-4" style="border: 1px solid #ddd; border-radius: 8px; overflow: hidden; margin-bottom: 20px;">
                            <div class="card-header" style="background-color: #f8f9fa; padding: 10px 15px; border-bottom: 1px solid #ddd; font-weight: bold;">
                                <h3 class="card-title" style="margin: 0; font-size: 1.1em;">Données Enthalpiques</h3>
                            </div>
                            <div class="card-body" style="padding: 15px;">
                                <div class="param-section">
                                    <h4>Enthalpies Liquide (kJ/mol)</h4>
                                    <div class="param-grid" style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px;">
                                        <div class="form-group">
                                            <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">hL (A pur)</label>
                                            <input type="number" class="form-input" id="hl-a" value="0" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                        </div>
                                        <div class="form-group">
                                            <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">hL (B pur)</label>
                                            <input type="number" class="form-input" id="hl-b" value="0" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                        </div>
                                    </div>
                                </div>
                                <div class="param-section">
                                    <h4>Enthalpies Vapeur (kJ/mol)</h4>
                                    <div class="param-grid" style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px;">
                                        <div class="form-group">
                                            <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">hV (A pur)</label>
                                            <input type="number" class="form-input" id="hv-a" value="35" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                        </div>
                                        <div class="form-group">
                                            <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">hV (B pur)</label>
                                            <input type="number" class="form-input" id="hv-b" value="40" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                        </div>
                                    </div>
                                </div>
                                <div class="form-group">
                                    <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">Enthalpie de mélange (kJ/mol)</label>
                                    <input type="number" class="form-input" id="h-mix" value="2" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                </div>
                            </div>
                        </div>
                    </div>
                    
                    <!-- MULTICOMP TAB -->
                    <div class="tab-panel" id="panel-multicomp">
                         <div class="card mb-4" style="border: 1px solid #ddd; border-radius: 8px; overflow: hidden; margin-bottom: 20px;">
                            <div class="card-header" style="background-color: #f8f9fa; padding: 10px 15px; border-bottom: 1px solid #ddd; font-weight: bold;">
                                <h3 class="card-title" style="margin: 0; font-size: 1.1em;">Distillation Multi-composants</h3>
                            </div>
                            <div class="card-body" style="padding: 15px;">
                                <div class="param-section">
                                    <h4>Composants & Volatilités</h4>
                                    <table class="mini-table" style="width:100%; font-size:0.9em; border-collapse: collapse;">
                                        <thead>
                                            <tr style="border-bottom: 1px solid #ddd;"><th style="padding: 5px;">Composant</th><th style="padding: 5px;">zF</th><th style="padding: 5px;">α</th><th style="padding: 5px;">Clé</th></tr>
                                        </thead>
                                        <tbody id="comp-table">
                                            <tr>
                                                <td style="padding: 5px;"><input type="text" class="form-input" value="A" style="width:60px; padding: 4px; border: 1px solid #ccc; border-radius: 3px;"></td>
                                                <td style="padding: 5px;"><input type="number" class="form-input" value="0.3" step="0.1" style="width:60px; padding: 4px; border: 1px solid #ccc; border-radius: 3px;"></td>
                                                <td style="padding: 5px;"><input type="number" class="form-input" value="3.0" step="0.1" style="width:60px; padding: 4px; border: 1px solid #ccc; border-radius: 3px;"></td>
                                                <td style="padding: 5px;"><select class="form-select" style="width:70px; padding: 4px; border: 1px solid #ccc; border-radius: 3px;"><option>LK</option><option>-</option><option>HK</option></select></td>
                                            </tr>
                                            <tr>
                                                <td style="padding: 5px;"><input type="text" class="form-input" value="B" style="width:60px; padding: 4px; border: 1px solid #ccc; border-radius: 3px;"></td>
                                                <td style="padding: 5px;"><input type="number" class="form-input" value="0.4" step="0.1" style="width:60px; padding: 4px; border: 1px solid #ccc; border-radius: 3px;"></td>
                                                <td style="padding: 5px;"><input type="number" class="form-input" value="1.5" step="0.1" style="width:60px; padding: 4px; border: 1px solid #ccc; border-radius: 3px;"></td>
                                                <td style="padding: 5px;"><select class="form-select" style="width:70px; padding: 4px; border: 1px solid #ccc; border-radius: 3px;"><option>LK</option><option selected>-</option><option>HK</option></select></td>
                                            </tr>
                                            <tr>
                                                <td style="padding: 5px;"><input type="text" class="form-input" value="C" style="width:60px; padding: 4px; border: 1px solid #ccc; border-radius: 3px;"></td>
                                                <td style="padding: 5px;"><input type="number" class="form-input" value="0.3" step="0.1" style="width:60px; padding: 4px; border: 1px solid #ccc; border-radius: 3px;"></td>
                                                <td style="padding: 5px;"><input type="number" class="form-input" value="1.0" step="0.1" style="width:60px; padding: 4px; border: 1px solid #ccc; border-radius: 3px;"></td>
                                                <td style="padding: 5px;"><select class="form-select" style="width:70px; padding: 4px; border: 1px solid #ccc; border-radius: 3px;"><option>LK</option><option>-</option><option selected>HK</option></select></td>
                                            </tr>
                                        </tbody>
                                    </table>
                                </div>
                                <div class="param-section">
                                    <h4>Récupérations</h4>
                                    <div class="param-grid" style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px;">
                                        <div class="form-group">
                                            <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">LK dans D (%)</label>
                                            <input type="number" class="form-input" id="lk-recovery" value="99" min="90" max="99.9" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                        </div>
                                        <div class="form-group">
                                            <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">HK dans B (%)</label>
                                            <input type="number" class="form-input" id="hk-recovery" value="99" min="90" max="99.9" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                    
                    <!-- HYDRAULICS TAB -->
                    <div class="tab-panel" id="panel-hydraulics">
                        <div class="card mb-4" style="border: 1px solid #ddd; border-radius: 8px; overflow: hidden; margin-bottom: 20px;">
                            <div class="card-header" style="background-color: #f8f9fa; padding: 10px 15px; border-bottom: 1px solid #ddd; font-weight: bold;">
                                <h3 class="card-title" style="margin: 0; font-size: 1.1em;">Propriétés Physiques</h3>
                            </div>
                            <div class="card-body" style="padding: 15px;">
                                <div class="param-grid" style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px;">
                                    <div class="form-group">
                                        <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">ρ vapeur (kg/m³)</label>
                                        <input type="number" class="form-input" id="rho-v" value="2.5" step="0.1" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                    </div>
                                    <div class="form-group">
                                        <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">ρ liquide (kg/m³)</label>
                                        <input type="number" class="form-input" id="rho-l" value="800" step="10" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                    </div>
                                    <div class="form-group">
                                        <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">σ tension surf. (N/m)</label>
                                        <input type="number" class="form-input" id="sigma" value="0.02" step="0.001" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                    </div>
                                    <div class="form-group">
                                        <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">Espacement plateaux (m)</label>
                                        <input type="number" class="form-input" id="tray-spacing" value="0.6" step="0.05" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                    </div>
                                    <div class="form-group">
                                        <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">Hauteur déversoir (m)</label>
                                        <input type="number" class="form-input" id="weir-height" value="0.05" step="0.01" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                    </div>
                                    <div class="form-group">
                                        <label style="display: block; margin-bottom: 5px; font-size: 0.9em;">Fraction aire trous</label>
                                        <input type="number" class="form-input" id="hole-area" value="0.1" step="0.01" style="width: 100%; padding: 6px; border: 1px solid #ccc; border-radius: 4px;">
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                    
                    <button class="btn btn-primary w-full" id="run-simulation" style="width: 100%; padding: 12px; background-color: #007bff; color: white; border: none; border-radius: 4px; font-size: 1.1em; cursor: pointer; transition: background-color 0.2s;">▶️ Lancer le Calcul</button>
                </div>
                
                <!-- RIGHT: Results -->
                <div>
                    <!-- Key Results -->
                    <div class="card mb-4" style="border: 1px solid #ddd; border-radius: 8px; overflow: hidden; margin-bottom: 20px;">
                        <div class="card-header d-flex justify-between" style="background-color: #f8f9fa; padding: 10px 15px; border-bottom: 1px solid #ddd; display: flex; justify-content: space-between; align-items: center;">
                            <h3 class="card-title" style="margin: 0; font-size: 1.1em;">Résultats Principaux</h3>
                            <span class="badge" id="method-badge" style="background-color: #6c757d; color: white; padding: 3px 8px; border-radius: 12px; font-size: 0.8em;">McCabe-Thiele</span>
                        </div>
                        <div class="card-body" style="padding: 15px;">
                            <div class="results-grid" style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 15px;">
                                <div class="result-item" style="text-align: center; padding: 10px; background-color: #f8f9fa; border-radius: 6px;">
                                    <div class="result-label" style="font-size: 0.85em; color: #666; margin-bottom: 5px;">N théoriques</div>
                                    <div class="result-value" id="result-stages" style="font-size: 1.2em; font-weight: bold; color: #007bff;">--</div>
                                </div>
                                <div class="result-item" style="text-align: center; padding: 10px; background-color: #f8f9fa; border-radius: 6px;">
                                    <div class="result-label" style="font-size: 0.85em; color: #666; margin-bottom: 5px;">N réels</div>
                                    <div class="result-value" id="result-actual" style="font-size: 1.2em; font-weight: bold; color: #007bff;">--</div>
                                </div>
                                <div class="result-item" style="text-align: center; padding: 10px; background-color: #f8f9fa; border-radius: 6px;">
                                    <div class="result-label" style="font-size: 0.85em; color: #666; margin-bottom: 5px;">Étage alim.</div>
                                    <div class="result-value" id="result-feed-stage" style="font-size: 1.2em; font-weight: bold; color: #007bff;">--</div>
                                </div>
                                <div class="result-item" style="text-align: center; padding: 10px; background-color: #f8f9fa; border-radius: 6px;">
                                    <div class="result-label" style="font-size: 0.85em; color: #666; margin-bottom: 5px;">Rmin</div>
                                    <div class="result-value" id="result-rmin" style="font-size: 1.2em; font-weight: bold; color: #28a745;">--</div>
                                </div>
                                <div class="result-item" style="text-align: center; padding: 10px; background-color: #f8f9fa; border-radius: 6px;">
                                    <div class="result-label" style="font-size: 0.85em; color: #666; margin-bottom: 5px;">R/Rmin</div>
                                    <div class="result-value" id="result-r-ratio" style="font-size: 1.2em; font-weight: bold; color: #28a745;">--</div>
                                </div>
                                <div class="result-item" style="text-align: center; padding: 10px; background-color: #f8f9fa; border-radius: 6px;">
                                    <div class="result-label" style="font-size: 0.85em; color: #666; margin-bottom: 5px;">D (kmol/h)</div>
                                    <div class="result-value" id="result-d" style="font-size: 1.2em; font-weight: bold; color: #17a2b8;">--</div>
                                </div>
                            </div>
                        </div>
                    </div>
                    
                    <!-- Hydraulics Results -->
                    <div class="card mb-4" id="hydraulics-results" style="border: 1px solid #ddd; border-radius: 8px; overflow: hidden; margin-bottom: 20px;">
                        <div class="card-header" style="background-color: #f8f9fa; padding: 10px 15px; border-bottom: 1px solid #ddd; font-weight: bold;">
                            <h3 class="card-title" style="margin: 0; font-size: 1.1em;">Dimensionnement Hydraulique</h3>
                        </div>
                        <div class="card-body" style="padding: 15px;">
                            <div class="hydraulics-grid">
                                <div class="hydraulic-item">
                                    <div class="value" id="res-diameter">--</div>
                                    <div class="label">Diamètre (m)</div>
                                </div>
                                <div class="hydraulic-item">
                                    <div class="value" id="res-flood">--</div>
                                    <div class="label">% Engorgement</div>
                                </div>
                                <div class="hydraulic-item">
                                    <div class="value" id="res-dp">--</div>
                                    <div class="label">ΔP total (kPa)</div>
                                </div>
                            </div>
                            <div class="mt-3" style="font-size: 0.85em; color: #666; margin-top: 10px;">
                                <span id="weeping-status">⏳ En attente de calcul...</span>
                            </div>
                        </div>
                    </div>
                    
                    <!-- Chart Area -->
                    <div class="card" style="border: 1px solid #ddd; border-radius: 8px; overflow: hidden; margin-bottom: 20px;">
                        <div class="card-header d-flex justify-between" style="background-color: #f8f9fa; padding: 10px 15px; border-bottom: 1px solid #ddd; display: flex; justify-content: space-between; align-items: center;">
                            <h3 class="card-title" style="margin: 0; font-size: 1.1em;">Diagramme</h3>
                            <div style="display: flex; gap: 5px;">
                                <button class="btn-sm" data-chart="xy" style="font-size:0.8em; padding: 2px 8px; border: 1px solid #ccc; border-radius: 4px; background: white; cursor: pointer;">x-y</button>
                                <button class="btn-sm" data-chart="hxy" style="font-size:0.8em; padding: 2px 8px; border: 1px solid #ccc; border-radius: 4px; background: white; cursor: pointer;">H-x-y</button>
                                <button class="btn-sm" data-chart="gamma" style="font-size:0.8em; padding: 2px 8px; border: 1px solid #ccc; border-radius: 4px; background: white; cursor: pointer;">γ</button>
                            </div>
                        </div>
                        <div class="card-body" style="padding: 15px; background: white;">
                            <canvas id="distill-canvas" style="width:100%; height:300px;"></canvas>
                        </div>
                    </div>
                    
                    <!-- Calculation Steps -->
                    <div class="card mt-4" style="border: 1px solid #ddd; border-radius: 8px; overflow: hidden;">
                        <div class="card-header" style="background-color: #f8f9fa; padding: 10px 15px; border-bottom: 1px solid #ddd;">
                            <!-- Detailed calculations button removed -->
                        </div>
                        <div class="card-body hidden" id="div-calcs-dist" style="padding: 15px;"></div>
                    </div>
                </div>
            </div>
        </div>
    `;
};

DistillationModule.init = async function() {
    setupDistillationEventListeners();
    runDistillationSimulation();
    console.log("Distillation Module Initialized with Event Listeners");

    // Load diagram
    const container = document.getElementById('distillation-diagram');
    if (container) {
        try {
            const response = await fetch('assets/distillation.svg');
            if (response.ok) {
                container.innerHTML = await response.text();
            }
        } catch (e) {
            container.innerHTML = '<p class="text-muted">Schéma non disponible</p>';
        }
    }
    
    // setupDistillationEventListeners(); // Removed duplicate call
    // runDistillationSimulation(); // Removed duplicate call
};

function setupDistillationEventListeners() {
    // Main simulation button
    document.getElementById('run-simulation')?.addEventListener('click', runDistillationSimulation);
    
    // Tab navigation
    document.querySelectorAll('.tab-btn[data-tab]').forEach(btn => {
        btn.addEventListener('click', (e) => {
            const tab = e.target.dataset.tab;
            // Update tab buttons
            document.querySelectorAll('.tab-btn[data-tab]').forEach(b => b.classList.remove('active'));
            e.target.classList.add('active');
            // Update panels
            document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
            const panel = document.getElementById('panel-' + tab);
            if (panel) panel.classList.add('active');
            // Update method badge
            const badge = document.getElementById('method-badge');
            if (badge) {
                const methodNames = {binary: 'McCabe-Thiele', ponchon: 'Ponchon-Savarit', multicomp: 'FUG', hydraulics: 'Hydraulique'};
                badge.textContent = methodNames[tab] || tab;
            }
            DistillationModule.currentTab = tab;
        });
    });
    
    // Column type
    document.querySelectorAll('[data-column]').forEach(btn => {
        btn.addEventListener('click', (e) => {
            document.querySelectorAll('[data-column]').forEach(b => b.classList.remove('active'));
            e.target.classList.add('active');
            DistillationModule.params.columnType = e.target.dataset.column;
            runDistillationSimulation();
        });
    });
    
    // VLE model selection
    document.getElementById('vle-model')?.addEventListener('change', (e) => {
        const showParams = e.target.value !== 'ideal';
        const paramsEl = document.getElementById('activity-params');
        if (paramsEl) paramsEl.style.display = showParams ? 'block' : 'none';
    });
    
    // Range sliders
    ['feed-comp', 'distillate-comp', 'bottoms-comp'].forEach(id => {
        const el = document.getElementById(id);
        if (el) {
            el.addEventListener('input', (e) => {
                const sibling = e.target.nextElementSibling;
                if (sibling) sibling.textContent = (e.target.value / 100).toFixed(2);
            });
        }
    });
    
    // Toggle calculations listener removed
    
    // Chart type buttons
    document.querySelectorAll('[data-chart]').forEach(btn => {
        btn.addEventListener('click', (e) => {
            drawDistillationChart(e.target.dataset.chart);
        });
    });
}

function runDistillationSimulation() {
    const currentTab = DistillationModule.currentTab || 'binary';
    
    // Parse multi-component table
    const comps = [];
    const zs = [];
    const alph = [];
    let lk = 0;
    let hk = 1;
    
    const rows = document.querySelectorAll('#comp-table tr');
    if (rows.length > 0) {
        rows.forEach((tr, idx) => {
            const inputs = tr.querySelectorAll('input, select');
            if (inputs.length >= 4) {
                 comps.push(inputs[0].value);
                 zs.push(parseFloat(inputs[1].value) || 0);
                 alph.push(parseFloat(inputs[2].value) || 1);
                 if (inputs[3].value === 'LK') lk = idx;
                 if (inputs[3].value === 'HK') hk = idx;
            }
        });
    }

    const params = {
        // Basic params
        feedComposition: parseFloat(document.getElementById('feed-comp')?.value || 50) / 100,
        distillateComposition: parseFloat(document.getElementById('distillate-comp')?.value || 95) / 100,
        bottomsComposition: parseFloat(document.getElementById('bottoms-comp')?.value || 5) / 100,
        feedCondition: 1.0,
        relativeVolatility: parseFloat(document.getElementById('alpha')?.value || 2.5),
        refluxRatio: parseFloat(document.getElementById('reflux')?.value || 2.0),
        columnType: document.querySelector('[data-column].active')?.dataset.column || 'tray',
        hetp: 0.5,
        feedFlowRate: parseFloat(document.getElementById('feed-flow')?.value || 100),
        trayEfficiency: parseFloat(document.getElementById('efficiency')?.value || 0.7),
        
        // VLE model
        vleModel: document.getElementById('vle-model')?.value || 'ideal',
        A12: parseFloat(document.getElementById('a12')?.value || 0.5),
        A21: parseFloat(document.getElementById('a21')?.value || 0.7),
        
        // Method
        method: currentTab === 'ponchon' ? 'ponchon' : 'mccabe',
        multiComponent: currentTab === 'multicomp',
        
        // Enthalpy (Ponchon)
        hL_A: parseFloat(document.getElementById('hl-a')?.value || 0),
        hL_B: parseFloat(document.getElementById('hl-b')?.value || 0),
        hV_A: parseFloat(document.getElementById('hv-a')?.value || 35),
        hV_B: parseFloat(document.getElementById('hv-b')?.value || 40),
        hMix: parseFloat(document.getElementById('h-mix')?.value || 2),
        
        // Hydraulics
        vaporDensity: parseFloat(document.getElementById('rho-v')?.value || 2.5),
        liquidDensity: parseFloat(document.getElementById('rho-l')?.value || 800),
        surfaceTension: parseFloat(document.getElementById('sigma')?.value || 0.02),
        traySpacing: parseFloat(document.getElementById('tray-spacing')?.value || 0.6),
        weirHeight: parseFloat(document.getElementById('weir-height')?.value || 0.05),
        holeArea: parseFloat(document.getElementById('hole-area')?.value || 0.1),
        
        // Multi-component 
        components: comps.length ? comps : ['A', 'B', 'C'],
        feedMoleFractions: zs.length ? zs : [0.3, 0.4, 0.3],
        alphas: alph.length ? alph : [3.0, 1.5, 1.0],
        lightKey: lk,
        heavyKey: hk,
        lkRecovery: parseFloat(document.getElementById('lk-recovery')?.value || 99) / 100,
        hkRecovery: parseFloat(document.getElementById('hk-recovery')?.value || 99) / 100
    };
    
    DistillationModule.params = params;
    const results = calculateDistillation(params);
    DistillationModule.results = results;
    
    if (results.error) {
        if (typeof showToast === 'function') showToast(results.message, 'error');
        return;
    }
    
    // Update main results
    const setVal = (id, val) => {
        const el = document.getElementById(id);
        if (el) el.textContent = val;
    };
    
    setVal('result-stages', results.theoreticalStages);
    setVal('result-actual', results.actualStages);
    setVal('result-feed-stage', results.feedStage);
    setVal('result-rmin', results.Rmin?.toFixed(3) || '--');
    setVal('result-r-ratio', results.RoverRmin?.toFixed(2) || '--');
    setVal('result-d', results.materialBalance?.D?.toFixed(1) || '--');
    
    // Hydraulics results
    if (results.hydraulics) {
        setVal('res-diameter', results.hydraulics.columnDiameter?.toFixed(2) || '--');
        setVal('res-flood', results.hydraulics.percentFlood || '80');
        setVal('res-dp', results.totalPressureDrop?.toFixed(2) || '--');
        
        const weepingEl = document.getElementById('weeping-status');
        if (weepingEl && results.hydraulics.weeping) {
            if (results.hydraulics.weeping.isWeeping) {
                weepingEl.innerHTML = '⚠️ <span style="color:red;">Risque de pleurage détecté!</span>';
            } else {
                weepingEl.innerHTML = '✅ Pas de pleurage (marge: ' + results.hydraulics.weeping.weepingMargin?.toFixed(0) + '%)';
            }
        }
    }
    
    // Calculation steps population removed
    
    // Draw chart
    drawDistillationChart('xy');
    
    if (typeof showToast === 'function') showToast('Calcul de distillation terminé', 'success');
}

/**
 * Draw distillation diagrams on canvas
 */
function drawDistillationChart(type) {
    const canvas = document.getElementById('distill-canvas');
    if (!canvas || !DistillationModule.results) return;
    
    const ctx = canvas.getContext('2d');
    const W = canvas.width = canvas.parentElement.clientWidth;
    const H = canvas.height = 300;
    
    ctx.clearRect(0, 0, W, H);
    
    const padding = 40;
    const params = DistillationModule.params;
    const results = DistillationModule.results;
    
    // Coordinate transforms
    const scaleX = (x) => padding + x * (W - 2 * padding);
    const scaleY = (y) => H - padding - y * (H - 2 * padding);
    
    // Draw axes
    ctx.strokeStyle = '#333';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(padding, padding);
    ctx.lineTo(padding, H - padding);
    ctx.lineTo(W - padding, H - padding);
    ctx.stroke();
    
    // Labels
    ctx.fillStyle = '#333';
    ctx.font = '12px sans-serif';
    ctx.fillText('x', W - padding - 10, H - padding + 15);
    ctx.fillText('y', padding - 15, padding + 10);
    
    if (type === 'xy') {
        // Draw diagonal
        ctx.strokeStyle = '#ccc';
        ctx.setLineDash([5, 5]);
        ctx.beginPath();
        ctx.moveTo(scaleX(0), scaleY(0));
        ctx.lineTo(scaleX(1), scaleY(1));
        ctx.stroke();
        ctx.setLineDash([]);
        
        // Draw equilibrium curve
        if (results.equilibriumCurve) {
            ctx.strokeStyle = '#2ecc71';
            ctx.lineWidth = 2;
            ctx.beginPath();
            results.equilibriumCurve.forEach((pt, i) => {
                const x = scaleX(pt.x);
                const y = scaleY(pt.y);
                if (i === 0) ctx.moveTo(x, y);
                else ctx.lineTo(x, y);
            });
            ctx.stroke();
            ctx.fillStyle = '#2ecc71';
            ctx.fillText('Équilibre', scaleX(0.7), scaleY(0.9));
        }
        
        // Draw operating lines
        if (results.operatingLines) {
            const { rectifying, stripping, qLine } = results.operatingLines;
            
            // Rectifying line
            ctx.strokeStyle = '#3498db';
            ctx.lineWidth = 1.5;
            ctx.beginPath();
            ctx.moveTo(scaleX(params.distillateComposition), scaleY(params.distillateComposition));
            ctx.lineTo(scaleX(qLine.x), scaleY(qLine.y));
            ctx.stroke();
            
            // Stripping line
            ctx.strokeStyle = '#e74c3c';
            ctx.beginPath();
            ctx.moveTo(scaleX(qLine.x), scaleY(qLine.y));
            ctx.lineTo(scaleX(params.bottomsComposition), scaleY(params.bottomsComposition));
            ctx.stroke();
            
            // q-line (vertical for q=1)
            ctx.strokeStyle = '#9b59b6';
            ctx.setLineDash([3, 3]);
            ctx.beginPath();
            ctx.moveTo(scaleX(params.feedComposition), scaleY(0));
            ctx.lineTo(scaleX(params.feedComposition), scaleY(1));
            ctx.stroke();
            ctx.setLineDash([]);
        }
        
        // Draw stages as steps
        if (results.stages) {
            ctx.strokeStyle = '#f39c12';
            ctx.lineWidth = 1;
            ctx.beginPath();
            let x = params.distillateComposition;
            let y = params.distillateComposition;
            ctx.moveTo(scaleX(x), scaleY(y));
            
            results.stages.forEach(stage => {
                // Horizontal to equilibrium
                ctx.lineTo(scaleX(stage.x), scaleY(stage.y));
                // Vertical to operating line
                const nextY = stage.x > params.feedComposition 
                    ? results.operatingLines.rectifying.slope * stage.x + results.operatingLines.rectifying.intercept
                    : results.operatingLines.stripping.slope * (stage.x - params.bottomsComposition) + params.bottomsComposition;
                ctx.lineTo(scaleX(stage.x), scaleY(Math.max(params.bottomsComposition, nextY)));
            });
            ctx.stroke();
        }
        
    } else if (type === 'gamma') {
        // Activity coefficient plot
        if (results.equilibriumCurve && params.vleModel !== 'ideal') {
            ctx.strokeStyle = '#e74c3c';
            ctx.lineWidth = 2;
            ctx.beginPath();
            const maxGamma = Math.max(...results.equilibriumCurve.map(p => Math.max(p.gamma1, p.gamma2)));
            results.equilibriumCurve.forEach((pt, i) => {
                const x = scaleX(pt.x);
                const y = H - padding - (pt.gamma1 / maxGamma) * (H - 2 * padding);
                if (i === 0) ctx.moveTo(x, y);
                else ctx.lineTo(x, y);
            });
            ctx.stroke();
            ctx.fillStyle = '#e74c3c';
            ctx.fillText('γ₁', scaleX(0.1), scaleY(0.8));
            
            ctx.strokeStyle = '#3498db';
            ctx.beginPath();
            results.equilibriumCurve.forEach((pt, i) => {
                const x = scaleX(pt.x);
                const y = H - padding - (pt.gamma2 / maxGamma) * (H - 2 * padding);
                if (i === 0) ctx.moveTo(x, y);
                else ctx.lineTo(x, y);
            });
            ctx.stroke();
            ctx.fillStyle = '#3498db';
            ctx.fillText('γ₂', scaleX(0.9), scaleY(0.8));
        } else {
            ctx.fillStyle = '#999';
            ctx.fillText('Sélectionnez un modèle non-idéal pour voir les γ', W/2 - 100, H/2);
        }
    } else if (type === 'hxy') {
        // Enthalpy-concentration diagram (simplified)
        ctx.fillStyle = '#999';
        ctx.fillText('Diagramme H-x-y (Ponchon-Savarit)', W/2 - 80, H/2);
        ctx.fillText('Sélectionnez l\'onglet Ponchon-Savarit', W/2 - 80, H/2 + 20);
    }
}

DistillationModule.getExplanation = function() {
    return {
        title: 'Simulation de Distillation Avancée',
        description: `
            Module complet de dimensionnement de colonnes de distillation avec 
            plusieurs méthodes de calcul et modèles thermodynamiques.
        `,
        theory: `
            <p><strong>Modèles VLE:</strong></p>
            <ul>
                <li>Idéal (Loi de Raoult)</li>
                <li>Margules, Van Laar (coefficients d'activité)</li>
                <li>Wilson, NRTL (modèles de composition locale)</li>
            </ul>
            <p><strong>Méthodes de calcul:</strong></p>
            <ul>
                <li>McCabe-Thiele (binaire, CMO)</li>
                <li>Ponchon-Savarit (enthalpie-concentration)</li>
                <li>FUG (multi-composants)</li>
            </ul>
        `,
        formulas: `
            <p><strong>Équilibre non-idéal:</strong></p>
            $$ y_i = \\frac{\\gamma_i x_i P_i^{sat}}{P} $$
            
            <p><strong>Margules:</strong></p>
            $$ \\ln\\gamma_1 = x_2^2[A_{12} + 2(A_{21}-A_{12})x_1] $$
            
            <p><strong>Fenske (Nmin):</strong></p>
            $$ N_{min} = \\frac{\\ln[(x_{LK}/x_{HK})_D \\cdot (x_{HK}/x_{LK})_B]}{\\ln(\\alpha_{LK}/\\alpha_{HK})} $$
        `,
        references: [
            'McCabe, Thiele "Graphical Design of Fractionating Columns" (1925)',
            'Ponchon, Savarit "Enthalpy-Concentration Method"',
            'Fenske, Underwood, Gilliland "Shortcut Methods"',
            'Fair "Tray Hydraulics Correlation"'
        ]
    };
};

if (typeof registerModule !== 'undefined') {
    registerModule('distillation', DistillationModule);
}

window.DistillationModule = DistillationModule;
