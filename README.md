# 🧪 Simulation Platform - Industrial Chemical Engineering Lab

A comprehensive, interactive educational simulation platform for industrial chemical engineering concepts. Built with pure front-end technologies (HTML, CSS, vanilla JavaScript) for easy deployment and accessibility.

![Platform Preview](assets/preview.png)

## 🎯 Features

### Core Functionality
- **7 Simulation Modules** covering thermodynamics, reaction engineering, and separation processes
- **Real-time calculations** with step-by-step educational explanations
- **Interactive charts** using Chart.js with SVG fallback
- **LaTeX equation rendering** via MathJax
- **Responsive design** - works on desktop, tablet, and mobile
- **Dark/Light theme** toggle with persistent preference
- **French/English** language support
- **Export capabilities** - CSV data and printable HTML reports

### Educational Modes
- **Student Mode**: Standard simulation with explanations
- **Assignment Mode**: Randomized problems for practice
- **Teacher Mode**: Solutions and advanced controls

## 📚 Modules

### 1. Compressor Simulation ✅ (Fully Implemented)
Single-stage centrifugal or piston compressor with:
- Isentropic and real (efficiency-corrected) calculations
- T-s and P-v diagram visualization
- Transient startup simulation using RK4 integration
- Power vs pressure ratio curves
- Multiple gas types (air, N₂, O₂, CO₂, CH₄, H₂, He, steam)

### 2. Rankine/Hirn Cycle 🚧 (Scaffold)
Steam power plant cycle including:
- Basic Rankine and superheated Hirn cycles
- Pump and turbine efficiency considerations
- Thermal efficiency calculations
- T-s diagram visualization

### 3. Combined Cycle (Brayton + Rankine) 🚧 (Scaffold)
Gas turbine combined cycle:
- Brayton (gas turbine) topping cycle
- Rankine (steam) bottoming cycle
- HRSG heat recovery modeling
- Overall efficiency > 55%

### 4. Chemical Reactors 🚧 (Scaffold)
Four fundamental reactor types:
- **Batch**: Time-dependent conversion
- **CSTR**: Continuous stirred tank with perfect mixing
- **PFR**: Plug flow with concentration profiles
- **PBR**: Packed bed with catalyst

Supports:
- First, second, zero-order kinetics
- Michaelis-Menten enzyme kinetics
- Arrhenius temperature dependence

### 5. Distillation Column 🚧 (Scaffold)
Binary distillation design:
- McCabe-Thiele graphical method
- Minimum reflux calculation
- Theoretical stage determination
- HETP for packed columns

### 6. Ester Saponification 🚧 (Scaffold)
Classic kinetics experiment:
- CH₃COOC₂H₅ + NaOH → CH₃COONa + C₂H₅OH
- Second-order reaction kinetics
- Conductivity-based conversion tracking
- Parameter estimation from experimental data

### 7. Chlorobenzene Synthesis 🚧 (Scaffold)
Industrial aromatic chlorination:
- Consecutive reaction analysis (mono → di → tri)
- Selectivity optimization
- Heat management
- Safety considerations

## 🚀 Getting Started

### Prerequisites
- Any modern web browser (Chrome, Firefox, Safari, Edge)
- No server required - runs entirely in the browser

### Installation

1. **Clone or download** this repository:
   ```bash
   git clone https://github.com/your-username/simulation-platform.git
   cd simulation-platform
   ```

2. **Open in browser**:
   - Simply double-click `index.html`, or
   - Use a local server for best experience:
     ```bash
     # Python 3
     python -m http.server 8080
     
     # Node.js
     npx serve .
     
     # VS Code Live Server extension
     ```

3. **Navigate** to `http://localhost:8080` (if using a server)

### Project Structure

```
simulation/
├── index.html          # Main HTML page
├── styles.css          # Complete CSS design system
├── app.js              # Application bootstrap, routing, i18n
├── assets/
│   ├── compressor.svg  # Interactive compressor diagram
│   ├── rankine.svg     # Rankine cycle diagram
│   ├── reactors.svg    # Reactor types comparison
│   └── distillation.svg# Distillation column diagram
├── modules/
│   ├── compressor.js   # ✅ Fully implemented
│   ├── rankine.js      # 🚧 Scaffold
│   ├── combined.js     # 🚧 Scaffold
│   ├── reactors.js     # 🚧 Scaffold
│   ├── distillation.js # 🚧 Scaffold
│   ├── ester.js        # 🚧 Scaffold
│   └── chlorobenzene.js# 🚧 Scaffold
├── utils/
│   └── math.js         # Numerical methods & thermodynamics
├── workers/
│   └── (future)        # Web Workers for heavy computation
└── README.md           # This file
```

## 🔧 Development

### Adding a New Module

1. Create `modules/your-module.js` following the pattern in existing modules
2. Implement the required interface:
   ```javascript
   const YourModule = {
       name: 'your-module',
       defaults: { /* default parameters */ },
       render: function() { return '<div>...</div>'; },
       init: async function() { /* setup */ },
       getExplanation: function(lang) { return { /* ... */ }; }
   };
   
   registerModule('your-module', YourModule);
   ```
3. Add navigation entry in `index.html`
4. Create SVG diagram in `assets/` (optional)

### Math Utilities Available

```javascript
// Numerical Integration
MathUtils.euler(f, y0, t0, tf, dt)      // Euler method
MathUtils.rk4(f, y0, t0, tf, dt)        // Runge-Kutta 4th order
MathUtils.eulerSystem(f, y0, t0, tf, dt) // Euler for systems
MathUtils.rk4System(f, y0, t0, tf, dt)   // RK4 for systems

// Thermodynamics
MathUtils.isentropicTemperature(T1, P1, P2, gamma)
MathUtils.actualTemperature(T1, T2s, efficiency)
MathUtils.shaftWork(m, cp, T1, T2)
MathUtils.idealGasDensity(P, R, T)
MathUtils.entropyChange(cp, R, T1, T2, P1, P2)

// Gas Properties
MathUtils.GAS_PROPERTIES.air    // {R, cp, cv, gamma, M}
MathUtils.GAS_PROPERTIES.steam
// etc.
```

### CSS Design System

The platform uses CSS custom properties for theming:

```css
:root {
    --primary: #008080;       /* Teal accent color */
    --background: #f5f5f5;    /* Light background */
    --surface: #ffffff;       /* Card background */
    --text: #333333;          /* Primary text */
    --text-muted: #666666;    /* Secondary text */
    --space-1: 4px;           /* Spacing scale */
    /* ... */
}

[data-theme="dark"] {
    --background: #1a1a2e;
    --surface: #16213e;
    --text: #e8e8e8;
    /* ... */
}
```

## 📖 Educational Theory

### Compressor Module (Example)

**Isentropic Compression:**
For an ideal gas undergoing reversible adiabatic compression:

$$T_{2s} = T_1 \left(\frac{P_2}{P_1}\right)^{\frac{\gamma-1}{\gamma}}$$

**Isentropic Efficiency:**
Relates actual work to ideal work:

$$\eta_s = \frac{W_s}{W_{actual}} = \frac{T_{2s} - T_1}{T_2 - T_1}$$

**Shaft Power:**
$$\dot{W} = \dot{m} \cdot c_p \cdot (T_2 - T_1)$$

## 🧪 Testing

The compressor module includes built-in unit tests:

```javascript
// Open browser console and run:
CompressorModule.runTests()
```

Expected output:
```
🧪 Running Compressor Module Tests...

Test 1: Isentropic compression (η = 1)
  Expected T2s: 409.XX K
  Calculated T2s: 409.XX K
  Match: ✅

Test 2: Real compression (η = 0.85)
  ...
  Match: ✅

🏁 Tests complete!
```

## 📱 Accessibility

- Full keyboard navigation support
- ARIA labels on all interactive elements
- High contrast color scheme
- Screen reader compatible
- Focus indicators

## 🌐 Browser Support

| Browser | Version |
|---------|---------|
| Chrome  | 80+     |
| Firefox | 75+     |
| Safari  | 13+     |
| Edge    | 80+     |

## 📄 License

This project is intended for educational purposes. Feel free to use, modify, and distribute for non-commercial educational use.

## 🤝 Contributing

Contributions are welcome! To complete the scaffold modules:

1. Fork the repository
2. Complete a module following the compressor.js pattern
3. Add proper unit tests
4. Submit a pull request

### Priority Items

- [ ] Implement accurate steam tables (IAPWS-IF97)
- [ ] Add McCabe-Thiele diagram visualization
- [ ] Complete transient reactor simulations
- [ ] Add Web Worker for heavy computations
- [ ] Implement parameter estimation algorithms

## 📚 References

- Çengel, Y.A. & Boles, M.A. "Thermodynamics: An Engineering Approach"
- Fogler, H.S. "Elements of Chemical Reaction Engineering"
- Seader, Henley "Separation Process Principles"
- Perry's Chemical Engineers' Handbook
- Smith, Van Ness, Abbott "Introduction to Chemical Engineering Thermodynamics"

## 👨‍💻 Author

Chemical Engineering Lab Simulation Platform  
Educational Tool for Industrial Chemical Engineering

---

**Note**: This is an educational simulation tool. Results should be validated against experimental data and industrial standards for real-world applications.
