/**
 * verify_resonance.js
 * 
 * Verifies that spectral_invariant_verifier.js produces the
 * correct Star Invariant ratio (approx 1.85731).
 */

const FluidResonanceSimulator = require('./spectral_invariant_verifier');

async function runVerification() {
    console.log("=== PROJECT ENTROPY: SCIENTIFIC VERIFICATION ===");
    const sim = new FluidResonanceSimulator();
    
    sim.injectPulse();
    const result = sim.step();
    
    if (result === "RESONANCE_MATCH") {
        console.log("\nVERIFICATION SUCCESS: The simulation is now REAL.");
        process.exit(0);
    } else {
        console.error("\nVERIFICATION FAILED: Resonance not detected.");
        process.exit(1);
    }
}

runVerification().catch(err => {
    console.error(err);
    process.exit(1);
});
