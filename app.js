/**
 * ============================================================================
 * MAIN APPLICATION - ROUTING AND CORE FUNCTIONALITY
 * ============================================================================
 * 
 * This is the main bootstrap and routing logic for the Chemical Engineering
 * Lab Simulation Platform. It handles:
 * - Module loading and switching
 * - Theme management (light/dark)
 * - Language switching (French/English)
 * - Mode management (Student/Assignment/Teacher)
 * - Modal dialogs
 * - Toast notifications
 * - Export functionality
 * - Assignment generator
 * 
 * @author Chemical Engineering Lab Simulation Platform
 * @version 1.0.1
 */

'use strict';

// ============================================================================
// APPLICATION STATE
// ============================================================================

const App = {
    // Current state
    state: {
        currentModule: null,
        mode: 'student',        // 'student', 'assignment', 'teacher'
        language: 'fr',          // 'fr', 'en'
        theme: 'light',          // 'light', 'dark'
        panelOpen: true,
        sidebarOpen: false       // For mobile
    },
    
    // Available modules registry
    modules: {},
    
    // Current simulation results
    currentResults: null,
    
    // Chart.js availability flag
    chartsAvailable: typeof Chart !== 'undefined'
};

// ============================================================================
// INTERNATIONALIZATION (i18n)
// ============================================================================

const i18n = {
    fr: {
        // Navigation
        'nav.compressor': 'Simulation d\'un compresseur',
        'nav.rankine': 'Centrale thermique (cycle de Hirn)',
        'nav.combined': 'Réacteurs + valorisation énergétique',
        'nav.reactors': 'Types de réacteurs (Batch, CSTR, PFR, PBR)',
        'nav.distillation': 'Colonnes de distillation',
        'nav.ester': 'Production d\'ester',
        'nav.chlorobenzene': 'Production de chlorobenzène',
        'nav.generator': 'Générateur d\'exercices',
        'nav.export': 'Exporter le rapport',
        
        // Common
        'common.run': 'Lancer la simulation',
        'common.reset': 'Réinitialiser',
        'common.example': 'Exemple',
        'common.showCalc': 'Afficher les calculs',
        'common.hideCalc': 'Masquer les calculs',
        'common.showSolution': 'Voir la solution',
        'common.export': 'Exporter',
        'common.loading': 'Chargement...',
        'common.inlet': 'Entrée',
        'common.outlet': 'Sortie',
        'common.power': 'Puissance',
        'common.efficiency': 'Rendement',
        'common.temperature': 'Température',
        'common.pressure': 'Pression',
        'common.massFlow': 'Débit massique',
        
        // Compressor specific
        'compressor.title': 'Simulation d\'un compresseur',
        'compressor.description': 'Modélisation d\'un compresseur centrifuge ou à piston avec calcul des conditions de sortie et de la puissance requise.',
        'compressor.inletTemp': 'Température d\'entrée',
        'compressor.inletPressure': 'Pression d\'entrée',
        'compressor.pressureRatio': 'Taux de compression',
        'compressor.isentropicEff': 'Rendement isentropique',
        'compressor.gasType': 'Type de gaz',
        'compressor.runTransient': 'Simulation transitoire',
        
        // Panel
        'panel.title': 'Explications',
        'panel.about': 'À propos',
        'panel.aboutText': 'Cette plateforme permet de simuler différents procédés de génie chimique industriel. Chaque module inclut des explications théoriques, des calculs détaillés et des visualisations interactives.',
        'panel.modes': 'Modes disponibles',
        'panel.formula': 'Formules principales',
        
        // Toasts
        'toast.success': 'Succès',
        'toast.error': 'Erreur',
        'toast.simulationComplete': 'Simulation terminée avec succès',
        'toast.exportComplete': 'Export terminé',
        'toast.copied': 'Copié dans le presse-papiers'
    },
    
    en: {
        'nav.compressor': 'Compressor simulation',
        'nav.rankine': 'Thermal power plant (Hirn cycle)',
        'nav.combined': 'Reactors + energy recovery',
        'nav.reactors': 'Reactor types (Batch, CSTR, PFR, PBR)',
        'nav.distillation': 'Distillation columns',
        'nav.ester': 'Ester production',
        'nav.chlorobenzene': 'Chlorobenzene production',
        'nav.generator': 'Assignment generator',
        'nav.export': 'Export report',
        
        'common.run': 'Run simulation',
        'common.reset': 'Reset',
        'common.example': 'Example',
        'common.showCalc': 'Show calculations',
        'common.hideCalc': 'Hide calculations',
        'common.showSolution': 'Show solution',
        'common.export': 'Export',
        'common.loading': 'Loading...',
        'common.inlet': 'Inlet',
        'common.outlet': 'Outlet',
        'common.power': 'Power',
        'common.efficiency': 'Efficiency',
        'common.temperature': 'Temperature',
        'common.pressure': 'Pressure',
        'common.massFlow': 'Mass flow',
        
        'compressor.title': 'Compressor Simulation',
        'compressor.description': 'Modeling of a centrifugal or piston compressor with outlet conditions and power requirement calculations.',
        'compressor.inletTemp': 'Inlet temperature',
        'compressor.inletPressure': 'Inlet pressure',
        'compressor.pressureRatio': 'Pressure ratio',
        'compressor.isentropicEff': 'Isentropic efficiency',
        'compressor.gasType': 'Gas type',
        'compressor.runTransient': 'Transient simulation',
        
        // Panel
        'panel.title': 'Explanations',
        'panel.about': 'About',
        'panel.aboutText': 'This platform allows you to simulate various industrial chemical engineering processes. Each module includes theoretical explanations, detailed calculations, and interactive visualizations.',
        'panel.modes': 'Available modes',
        'panel.formula': 'Main formulas',
        
        'toast.success': 'Success',
        'toast.error': 'Error',
        'toast.simulationComplete': 'Simulation completed successfully',
        'toast.exportComplete': 'Export complete',
        'toast.copied': 'Copied to clipboard'
    }
};

/**
 * Get translated string for current language
 * @param {string} key - Translation key
 * @returns {string} Translated string or key if not found
 */
function t(key) {
    const lang = App.state.language;
    return i18n[lang]?.[key] || i18n['fr']?.[key] || key;
}

/**
 * Update all translatable elements in the DOM
 */
function updateTranslations() {
    document.querySelectorAll('[data-i18n]').forEach(el => {
        const key = el.getAttribute('data-i18n');
        el.textContent = t(key);
    });
}

// ============================================================================
// THEME MANAGEMENT
// ============================================================================

/**
 * Set the application theme
 * @param {string} theme - 'light' or 'dark'
 */
function setTheme(theme) {
    App.state.theme = theme;
    document.documentElement.setAttribute('data-theme', theme);
    localStorage.setItem('theme', theme);
    
    // Update theme toggle
    const toggle = document.getElementById('theme-toggle');
    if (toggle) {
        toggle.checked = theme === 'dark';
    }
}

/**
 * Initialize theme from localStorage or system preference
 */
function initTheme() {
    const saved = localStorage.getItem('theme');
    if (saved) {
        setTheme(saved);
    } else if (window.matchMedia('(prefers-color-scheme: dark)').matches) {
        setTheme('dark');
    }
}

// ============================================================================
// TOAST NOTIFICATIONS
// ============================================================================

/**
 * Show a toast notification
 * @param {string} message - The message to display
 * @param {string} type - 'success', 'error', 'warning', 'info'
 * @param {number} duration - Duration in ms (default: 3000)
 */
function showToast(message, type = 'success', duration = 3000) {
    const container = document.getElementById('toast-container');
    
    const toast = document.createElement('div');
    toast.className = `toast ${type}`;
    
    const icons = {
        success: '✓',
        error: '✕',
        warning: '⚠',
        info: 'ℹ'
    };
    
    toast.innerHTML = `
        <span class="toast-icon">${icons[type] || icons.info}</span>
        <span class="toast-message">${message}</span>
        <button class="toast-close" aria-label="Fermer">&times;</button>
    `;
    
    container.appendChild(toast);
    
    // Close button handler
    toast.querySelector('.toast-close').addEventListener('click', () => {
        toast.remove();
    });
    
    // Auto-remove after duration
    setTimeout(() => {
        toast.style.animation = 'slideIn 0.3s ease reverse';
        setTimeout(() => toast.remove(), 300);
    }, duration);
}

// ============================================================================
// MODAL MANAGEMENT
// ============================================================================

/**
 * Open a modal dialog
 * @param {string} modalId - The ID of the modal element
 */
function openModal(modalId) {
    const modal = document.getElementById(modalId);
    if (modal) {
        modal.classList.add('open');
        modal.setAttribute('aria-hidden', 'false');
        
        // Focus first focusable element
        const focusable = modal.querySelector('button, input, select, textarea');
        if (focusable) focusable.focus();
        
        // Trap focus within modal
        modal.addEventListener('keydown', trapFocus);
    }
}

/**
 * Close a modal dialog
 * @param {string} modalId - The ID of the modal element
 */
function closeModal(modalId) {
    const modal = document.getElementById(modalId);
    if (modal) {
        modal.classList.remove('open');
        modal.setAttribute('aria-hidden', 'true');
        modal.removeEventListener('keydown', trapFocus);
    }
}

/**
 * Close all open modals
 */
function closeAllModals() {
    document.querySelectorAll('.modal-overlay.open').forEach(modal => {
        modal.classList.remove('open');
        modal.setAttribute('aria-hidden', 'true');
    });
}

/**
 * Trap focus within modal for accessibility
 */
function trapFocus(e) {
    if (e.key !== 'Tab') return;
    
    const modal = e.currentTarget;
    const focusable = modal.querySelectorAll('button, input, select, textarea, [tabindex]:not([tabindex="-1"])');
    const first = focusable[0];
    const last = focusable[focusable.length - 1];
    
    if (e.shiftKey && document.activeElement === first) {
        e.preventDefault();
        last.focus();
    } else if (!e.shiftKey && document.activeElement === last) {
        e.preventDefault();
        first.focus();
    }
}

// ============================================================================
// MODULE LOADING AND ROUTING
// ============================================================================

/**
 * Register a simulation module
 * @param {string} name - Module identifier
 * @param {Object} module - Module object with render() and init() methods
 */
function registerModule(name, module) {
    App.modules[name] = module;
}

/**
 * Import modules that were registered before app.js loaded
 */
function importPreRegisteredModules() {
    // 1. Try registry
    if (window.ModuleRegistry) {
        Object.keys(window.ModuleRegistry).forEach(name => {
            App.modules[name] = window.ModuleRegistry[name];
        });
    }
    
    // 2. Fallback to globals if missing (for robustness)
    const fallbackMap = {
        'compressor': window.CompressorModule,
        'rankine': window.RankineModule,
        'combined': window.CombinedModule,
        'reactors': window.ReactorsModule,
        'distillation': window.DistillationModule,
        'ester': window.EsterModule,
        'chlorobenzene': window.ChlorobenzeneModule
    };
    
    Object.keys(fallbackMap).forEach(name => {
        if (!App.modules[name] && fallbackMap[name]) {
            console.warn(`Module ${name} recovered from global scope`);
            App.modules[name] = fallbackMap[name];
        }
    });
}

/**
 * Load and display a module
 * @param {string} moduleName - The name of the module to load
 */
async function loadModule(moduleName) {
    console.log(` Attempting to load module: ${moduleName}`);
    
    const container = document.getElementById('module-container');
    const loadingState = document.getElementById('loading-state');
    
    // Show loading
    loadingState.style.display = 'block';
    document.getElementById('welcome-content')?.remove();
    
    // Update navigation state
    document.querySelectorAll('.nav-item').forEach(item => {
        const isActive = item.dataset.module === moduleName;
        item.classList.toggle('active', isActive);
        item.setAttribute('aria-current', isActive ? 'page' : 'false');
    });
    
    // Close mobile sidebar
    if (window.innerWidth < 768) {
        document.getElementById('sidebar')?.classList.remove('open');
    }
    
    try {
        // Ensure registry is up to date
        if (window.ModuleRegistry && window.ModuleRegistry[moduleName]) {
             App.modules[moduleName] = window.ModuleRegistry[moduleName];
        }

        let module = App.modules[moduleName];

        // Fallback: Check global objects (e.g. RankineModule) if registry failed
        if (!module) {
             const pascalCase = moduleName.charAt(0).toUpperCase() + moduleName.slice(1) + 'Module';
             if (window[pascalCase]) {
                 console.log(`Recovered module "${moduleName}" from global scope: ${pascalCase}`);
                 module = window[pascalCase];
                 App.modules[moduleName] = module;
             }
        }
        
        if (!module) {
            throw new Error(`Module "${moduleName}" not found in registry or global scope`);
        }
        
        if (typeof module.render !== 'function') {
             // Fallback for modules without render (scaffolds)
             if (module.defaults) { 
                 console.warn(`Module "${moduleName}" has no render method. Using default maintenance view.`);
                 module.render = () => `
                    <div class="card">
                        <div class="card-body" style="text-align: center; padding: 40px;">
                            <h3>🚧 ${module.name || moduleName} under construction</h3>
                            <p>This module exists but has no render method.</p>
                        </div>
                    </div>`;
             } else {
                 throw new Error(`Module "${moduleName}" is valid but missing render() method`);
             }
        }

        // Hide loading
        loadingState.style.display = 'none';
        
        // Render module content
        container.innerHTML = module.render();
        
        // Initialize module
        if (module.init) {
            // Await init to ensure event listeners are bound
            await new Promise(resolve => setTimeout(resolve, 0)); 
            await module.init();
        }
        
        // Update right panel with module-specific content
        updateExplanationPanel(moduleName);
        
        // Update state
        App.state.currentModule = moduleName;
        
        // Update URL hash for bookmarking
        window.location.hash = moduleName;
        
        // Re-render MathJax if available
        if (window.MathJax?.typesetPromise) {
            await MathJax.typesetPromise();
        }
        
        // Update translations for the newly loaded module content
        updateTranslations();
        
    } catch (error) {
        console.error('Error loading module:', error);
        if (typeof showToast === 'function') {
            showToast(`Error: ${error.message}`, 'error');
        }
        loadingState.style.display = 'none';
        container.innerHTML = `
            <div class="card">
                <div class="card-body" style="text-align: center; padding: var(--space-8);">
                    <h2 style="color: var(--danger);">⚠️ Erreur de chargement</h2>
                    <p class="text-danger">${error.message}</p>
                    <p class="text-muted small">Check console for details (F12)</p>
                    <button class="btn btn-primary mt-4" onclick="loadModule('compressor')">
                        Retour au compresseur
                    </button>
                </div>
            </div>
        `;
    }
}

/**
 * Update the explanation panel based on the current module
 * @param {string} moduleName - Current module name
 */
function updateExplanationPanel(moduleName) {
    const module = App.modules[moduleName];
    const content = document.getElementById('explanation-content');
    const formulas = document.getElementById('current-formulas');
    
    if (module?.getExplanation) {
        const explanation = module.getExplanation(App.state.language);
        
        content.innerHTML = `
            <div class="explanation-section">
                <h4>${explanation.title}</h4>
                <p class="text-muted">${explanation.description}</p>
            </div>
            
            ${explanation.theory ? `
            <div class="explanation-section mt-5">
                <h4>📖 Théorie</h4>
                <div class="text-muted">${explanation.theory}</div>
            </div>
            ` : ''}
            
            <div class="explanation-section mt-5">
                <h4>📐 Formules principales</h4>
                <div class="equations-box">
                    ${explanation.formulas}
                </div>
            </div>
            
            ${explanation.references ? `
            <div class="explanation-section mt-5">
                <h4>📚 Références</h4>
                <ul class="text-muted" style="padding-left: var(--space-5); font-size: var(--font-size-sm);">
                    ${explanation.references.map(ref => `<li>${ref}</li>`).join('')}
                </ul>
            </div>
            ` : ''}
        `;
    }
}

// ============================================================================
// EXPORT FUNCTIONALITY
// ============================================================================

/**
 * Export simulation results to CSV
 * @param {Object} data - Data to export
 * @param {string} filename - Filename without extension
 */
function exportToCSV(data, filename = 'simulation_results') {
    if (!data || !data.headers || !data.rows) {
        showToast('Aucune donnée à exporter', 'warning');
        return;
    }
    
    // Build CSV content
    let csv = data.headers.join(',') + '\n';
    data.rows.forEach(row => {
        csv += row.map(cell => {
            // Escape quotes and wrap in quotes if contains comma
            if (typeof cell === 'string' && (cell.includes(',') || cell.includes('"'))) {
                return `"${cell.replace(/"/g, '""')}"`;
            }
            return cell;
        }).join(',') + '\n';
    });
    
    // Download
    downloadFile(csv, `${filename}.csv`, 'text/csv');
    showToast(t('toast.exportComplete'), 'success');
}

/**
 * Export simulation report as HTML
 * @param {Object} reportData - Report data
 * @param {string} filename - Filename without extension
 */
function exportToHTML(reportData, filename = 'simulation_report') {
    const html = `
<!DOCTYPE html>
<html lang="fr">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Rapport de Simulation - ${reportData.title}</title>
    <style>
        body { font-family: Arial, sans-serif; max-width: 800px; margin: 40px auto; padding: 20px; }
        h1 { color: #008080; border-bottom: 2px solid #008080; padding-bottom: 10px; }
        h2 { color: #333; margin-top: 30px; }
        table { width: 100%; border-collapse: collapse; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 10px; text-align: left; }
        th { background: #008080; color: white; }
        .param { background: #f5f5f5; }
        .result { font-weight: bold; color: #008080; }
        .formula { background: #f0f8ff; padding: 15px; border-left: 4px solid #008080; margin: 20px 0; font-family: monospace; }
        .footer { margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd; color: #666; font-size: 12px; }
        @media print { body { margin: 0; } }
    </style>
</head>
<body>
    <h1>📊 ${reportData.title}</h1>
    <p><strong>Date:</strong> ${new Date().toLocaleString('fr-FR')}</p>
    <p><strong>Module:</strong> ${reportData.module}</p>
    
    <h2>Paramètres d'entrée</h2>
    <table>
        <tr><th>Paramètre</th><th>Valeur</th><th>Unité</th></tr>
        ${reportData.inputs.map(p => `
            <tr class="param">
                <td>${p.name}</td>
                <td>${p.value}</td>
                <td>${p.unit}</td>
            </tr>
        `).join('')}
    </table>
    
    <h2>Résultats</h2>
    <table>
        <tr><th>Résultat</th><th>Valeur</th><th>Unité</th></tr>
        ${reportData.outputs.map(r => `
            <tr>
                <td>${r.name}</td>
                <td class="result">${r.value}</td>
                <td>${r.unit}</td>
            </tr>
        `).join('')}
    </table>
    
    ${reportData.calculations ? `
    <h2>Calculs détaillés</h2>
    ${reportData.calculations.map((step, i) => `
        <div class="formula">
            <strong>Étape ${i + 1}:</strong> ${step.description}<br>
            ${step.formula}<br>
            <em>Résultat: ${step.result}</em>
        </div>
    `).join('')}
    ` : ''}
    
    <div class="footer">
        <p>Généré par le Laboratoire de Génie Chimique - Plateforme de Simulation</p>
        <p>Ce rapport est fourni à des fins éducatives uniquement.</p>
    </div>
</body>
</html>
    `;
    
    downloadFile(html, `${filename}.html`, 'text/html');
    showToast(t('toast.exportComplete'), 'success');
}

/**
 * Download a file
 * @param {string} content - File content
 * @param {string} filename - Filename with extension
 * @param {string} mimeType - MIME type
 */
function downloadFile(content, filename, mimeType) {
    const blob = new Blob([content], { type: mimeType });
    const url = URL.createObjectURL(blob);
    
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}

// ============================================================================
// ASSIGNMENT GENERATOR
// ============================================================================

/**
 * Generate an assignment JSON
 */
function generateAssignment() {
    const form = document.getElementById('assignment-form');
    const module = document.getElementById('assignment-module').value;
    const objective = document.getElementById('assignment-objective').value;
    const constraints = document.getElementById('assignment-constraints').value;
    const difficulty = form.querySelector('.mode-btn.active')?.dataset.difficulty || 'medium';
    
    const assignment = {
        id: `assignment_${Date.now()}`,
        module: module,
        difficulty: difficulty,
        objective: objective,
        constraints: constraints.split('\n').filter(c => c.trim()),
        createdAt: new Date().toISOString(),
        parameters: App.modules[module]?.getDefaultParameters?.() || {}
    };
    
    const json = JSON.stringify(assignment, null, 2);
    downloadFile(json, `assignment_${module}_${difficulty}.json`, 'application/json');
    closeModal('assignment-modal');
    showToast('Exercice exporté avec succès', 'success');
}

// ============================================================================
// CHART HELPERS
// ============================================================================

/**
 * Create a Chart.js chart or fallback to SVG
 * @param {string} containerId - Container element ID
 * @param {Object} config - Chart configuration
 * @returns {Object|null} Chart instance or null
 */
function createChart(containerId, config) {
    const container = document.getElementById(containerId);
    if (!container) return null;
    
    // If Chart.js is available, use it
    if (App.chartsAvailable) {
        // Clear any existing chart
        const existingCanvas = container.querySelector('canvas');
        if (existingCanvas) {
            Chart.getChart(existingCanvas)?.destroy();
            existingCanvas.remove();
        }
        
        const canvas = document.createElement('canvas');
        container.appendChild(canvas);
        
        return new Chart(canvas, config);
    }
    
    // Fallback to simple SVG chart
    return createSVGChart(container, config);
}

/**
 * Create a simple SVG line chart (fallback when Chart.js unavailable)
 * @param {HTMLElement} container - Container element
 * @param {Object} config - Chart configuration
 */
function createSVGChart(container, config) {
    const width = container.clientWidth || 400;
    const height = container.clientHeight || 300;
    const padding = { top: 30, right: 30, bottom: 40, left: 60 };
    
    const rawData = config.data.datasets[0]?.data || [];
    
    if (rawData.length === 0) {
        container.innerHTML = '<p class="text-muted text-center">Aucune donnée</p>';
        return null;
    }
    
    // Normalize data to points {x, y}
    const points = rawData.map((d, i) => {
        if (typeof d === 'number') return { x: i, y: d };
        // Handle {x,y} or {t,y} or fallback to index
        const x = d.x !== undefined ? d.x : (d.t !== undefined ? d.t : i);
        const y = d.y !== undefined ? d.y : d;
        return { x, y };
    });

    const xValues = points.map(p => p.x);
    const yValues = points.map(p => p.y);
    
    const xMin = Math.min(...xValues);
    const xMax = Math.max(...xValues);
    // Ensure we don't have NaN for min/max
    const yMinCalc = Math.min(...yValues);
    const yMaxCalc = Math.max(...yValues);
    
    const yMin = isNaN(yMinCalc) ? 0 : yMinCalc * 0.9;
    const yMax = isNaN(yMaxCalc) ? 100 : yMaxCalc * 1.1;
    
    const xRange = (xMax - xMin) || 1;
    const yRange = (yMax - yMin) || 1;
    
    const scaleX = (v) => padding.left + ((v - xMin) / xRange) * (width - padding.left - padding.right);
    const scaleY = (v) => height - padding.bottom - ((v - yMin) / yRange) * (height - padding.top - padding.bottom);
    
    const pathD = points.map((p, i) => `${i === 0 ? 'M' : 'L'} ${scaleX(p.x)} ${scaleY(p.y)}`).join(' ');
    
    const svg = `
        <svg viewBox="0 0 ${width} ${height}" style="width: 100%; height: 100%;">
            <!-- Grid lines -->
            <g stroke="#e0e0e0" stroke-width="1">
                ${[0, 0.25, 0.5, 0.75, 1].map(f => {
                    const y = padding.top + f * (height - padding.top - padding.bottom);
                    return `<line x1="${padding.left}" y1="${y}" x2="${width - padding.right}" y2="${y}"/>`;
                }).join('')}
            </g>
            
            <!-- Axes -->
            <line x1="${padding.left}" y1="${height - padding.bottom}" x2="${width - padding.right}" y2="${height - padding.bottom}" stroke="#333" stroke-width="2"/>
            <line x1="${padding.left}" y1="${padding.top}" x2="${padding.left}" y2="${height - padding.bottom}" stroke="#333" stroke-width="2"/>
            
            <!-- Data line -->
            <path d="${pathD}" fill="none" stroke="#008080" stroke-width="2"/>
            
            <!-- Data points -->
            ${points.map(p => `<circle cx="${scaleX(p.x)}" cy="${scaleY(p.y)}" r="4" fill="#008080"/>`).join('')}
            
            <!-- Y axis labels -->
            ${[0, 0.5, 1].map(f => {
                const y = padding.top + f * (height - padding.top - padding.bottom);
                const val = yMax - f * (yMax - yMin);
                return `<text x="${padding.left - 10}" y="${y + 4}" text-anchor="end" font-size="11" fill="#666">${val.toFixed(1)}</text>`;
            }).join('')}
            
            <!-- Title -->
            <text x="${width / 2}" y="20" text-anchor="middle" font-size="14" font-weight="bold" fill="#333">
                ${config.options?.plugins?.title?.text || 'Graphique'}
            </text>
        </svg>
    `;
    
    container.innerHTML = svg;
    return { svg: true, container };
}

/**
 * Update an existing chart with new data
 * @param {Object} chart - Chart instance
 * @param {Array} labels - New labels
 * @param {Array} data - New data
 */
function updateChart(chart, labels, data) {
    if (!chart) return;
    
    if (chart.svg) {
        // Re-render SVG chart
        createSVGChart(chart.container, {
            data: { labels, datasets: [{ data }] }
        });
    } else {
        // Update Chart.js chart
        chart.data.labels = labels;
        chart.data.datasets[0].data = data;
        chart.update('none');
    }
}

// ============================================================================
// EVENT HANDLERS
// ============================================================================

/**
 * Toggle the right panel open/closed
 */
function toggleRightPanel() {
    const panel = document.getElementById('right-panel');
    const content = document.querySelector('.content-area');
    const openBtn = document.getElementById('panel-open-btn');
    
    if (!panel || !content) return;
    
    panel.classList.toggle('open');
    content.classList.toggle('panel-open');
    App.state.panelOpen = panel.classList.contains('open');
    
    // Show/hide the open button
    if (openBtn) {
        openBtn.style.display = App.state.panelOpen ? 'none' : 'flex';
    }
}

function initEventListeners() {
    // Theme toggle
    document.getElementById('theme-toggle')?.addEventListener('change', (e) => {
        setTheme(e.target.checked ? 'dark' : 'light');
    });
    
    // Language selector
    document.getElementById('language-select')?.addEventListener('change', (e) => {
        App.state.language = e.target.value;
        localStorage.setItem('language', e.target.value);
        updateTranslations();
        
        // Reload current module to update content
        if (App.state.currentModule) {
            loadModule(App.state.currentModule);
        }
    });
    
    // Mode selector
    document.querySelectorAll('.mode-selector .mode-btn').forEach(btn => {
        btn.addEventListener('click', (e) => {
            const parent = e.target.closest('.mode-selector');
            parent.querySelectorAll('.mode-btn').forEach(b => {
                b.classList.remove('active');
                b.setAttribute('aria-selected', 'false');
            });
            e.target.classList.add('active');
            e.target.setAttribute('aria-selected', 'true');
            
            const mode = e.target.dataset.mode;
            if (mode) {
                App.state.mode = mode;
                document.body.setAttribute('data-mode', mode);
                
                // Notify current module of mode change
                const module = App.modules[App.state.currentModule];
                if (module?.onModeChange) {
                    module.onModeChange(mode);
                }
            }
        });
    });
    
    // Navigation items (Handled by onclick in HTML for robustness)
    // document.querySelectorAll('.nav-item[data-module]').forEach(item => {
    //     item.addEventListener('click', () => {
    //         loadModule(item.dataset.module);
    //     });
    // });
    
    // Welcome button
    document.querySelector('[data-module="compressor"]:not(.nav-item)')?.addEventListener('click', () => {
        loadModule('compressor');
    });
    
    // Sidebar toggle (Desktop & Mobile)
    document.getElementById('menu-toggle')?.addEventListener('click', () => {
        const sidebar = document.getElementById('sidebar');
        const content = document.getElementById('main-content');
        
        if (window.innerWidth < 768) {
             // Mobile
            sidebar.classList.toggle('open');
            App.state.sidebarOpen = sidebar.classList.contains('open');
        } else {
             // Tablet/Desktop
            sidebar.classList.toggle('closed');
            content.classList.toggle('sidebar-closed');
            App.state.sidebarOpen = !sidebar.classList.contains('closed');
        }
    });
    
    // Right panel toggle
    document.getElementById('panel-toggle')?.addEventListener('click', toggleRightPanel);
    document.getElementById('panel-open-btn')?.addEventListener('click', toggleRightPanel);
    
    // Modal close buttons
    document.querySelectorAll('.modal-close, [data-action="close-modal"]').forEach(btn => {
        btn.addEventListener('click', (e) => {
            const modal = e.target.closest('.modal-overlay');
            if (modal) closeModal(modal.id);
        });
    });
    
    // Modal overlay click to close
    document.querySelectorAll('.modal-overlay').forEach(overlay => {
        overlay.addEventListener('click', (e) => {
            if (e.target === overlay) {
                closeModal(overlay.id);
            }
        });
    });
    
    // Assignment generator
    document.querySelector('[data-action="assignment-generator"]')?.addEventListener('click', () => {
        openModal('assignment-modal');
    });
    
    document.getElementById('generate-assignment')?.addEventListener('click', generateAssignment);
    
    // Export report
    document.querySelector('[data-action="export-report"]')?.addEventListener('click', () => {
        openModal('export-modal');
    });
    
    document.getElementById('do-export')?.addEventListener('click', () => {
        const format = document.querySelector('input[name="export-format"]:checked')?.value || 'csv';
        const module = App.modules[App.state.currentModule];
        
        if (module?.getExportData) {
            const data = module.getExportData();
            
            if (format === 'csv') {
                exportToCSV(data.csv, `${App.state.currentModule}_results`);
            } else {
                exportToHTML(data.report, `${App.state.currentModule}_report`);
            }
        } else {
            showToast('Aucune donnée à exporter', 'warning');
        }
        
        closeModal('export-modal');
    });
    
    // Print solution
    document.getElementById('print-solution')?.addEventListener('click', () => {
        window.print();
    });
    
    // Keyboard navigation
    document.addEventListener('keydown', (e) => {
        // Escape closes modals
        if (e.key === 'Escape') {
            closeAllModals();
        }
        
        // Ctrl+E exports
        if (e.ctrlKey && e.key === 'e') {
            e.preventDefault();
            openModal('export-modal');
        }
    });
    
    // Handle URL hash for direct linking
    window.addEventListener('hashchange', () => {
        const module = window.location.hash.slice(1);
        if (module && App.modules[module] && module !== App.state.currentModule) {
            loadModule(module);
        }
    });
    
    // Close sidebar when clicking outside on mobile
    document.addEventListener('click', (e) => {
        if (window.innerWidth < 768) {
            const sidebar = document.getElementById('sidebar');
            const menuToggle = document.getElementById('menu-toggle');
            
            if (!sidebar.contains(e.target) && !menuToggle.contains(e.target)) {
                sidebar.classList.remove('open');
            }
        }
    });
}

// ============================================================================
// INITIALIZATION
// ============================================================================

/**
 * Initialize the application
 */
async function initApp() {
    console.log('🔬 Chemical Engineering Lab Simulation Platform v1.0.1');
    console.log('📚 Starting initialization...');
    
    // Import modules that were registered before app.js loaded
    importPreRegisteredModules();
    
    // Initialize theme
    initTheme();
    
    // Initialize language
    const savedLang = localStorage.getItem('language') || 'fr';
    App.state.language = savedLang;
    document.getElementById('language-select').value = savedLang;
    updateTranslations();
    
    // Setup event listeners
    initEventListeners();
    
    // Check Chart.js availability
    if (!App.chartsAvailable) {
        console.warn('⚠️ Chart.js not available. Using SVG fallback for charts.');
    }

    // Initialize Web Worker
    if (window.Worker) {
        try {
            App.worker = new Worker('workers/simulation-worker.js');
            App.worker.onmessage = function(e) {
                console.log('⚙️ Worker:', e.data);
            };
            // Ping the worker to ensure it's alive
            App.worker.postMessage({ type: 'ping', id: 'init' });
            console.log('✅ Simulation Worker initialized');
        } catch (err) {
            console.error('Failed to initialize worker:', err);
        }
    }
    
    // Load module from URL hash or default to compressor
    const hash = window.location.hash.slice(1);
    
    // Attempt to force-load all modules into App.modules before initial routing
    importPreRegisteredModules();

    // Determine initial module
    let moduleToLoad = 'compressor';
    
    // If hash exists and is a valid module, use it
    if (hash && (App.modules[hash] || window.ModuleRegistry?.[hash])) {
        moduleToLoad = hash;
    }

    console.log(`🚀 Booting into module: ${moduleToLoad}`);
    
    // Wait a brief moment for scripts to settle if needed, then load
    setTimeout(() => {
        loadModule(moduleToLoad).catch(e => console.error("Initial load failed", e));
    }, 100);
    
    console.log('✅ Application initialized successfully');
}

// Start application when DOM is ready
document.addEventListener('DOMContentLoaded', initApp);

// ============================================================================
// GLOBAL EXPORTS
// ============================================================================

// Make functions available globally for module usage
window.App = App;
window.registerModule = registerModule;
window.loadModule = loadModule;
window.showToast = showToast;
window.openModal = openModal;
window.closeModal = closeModal;
window.createChart = createChart;
window.updateChart = updateChart;
window.exportToCSV = exportToCSV;
window.exportToHTML = exportToHTML;
window.toggleRightPanel = toggleRightPanel;
window.t = t;
