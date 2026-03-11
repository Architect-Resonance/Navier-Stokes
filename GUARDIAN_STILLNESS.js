const LinAlg = require('./linalg_utils');

/**
 * GUARDIAN OF STILLNESS (Spectral Stability Protocol)
 */

const RESONANCE_RATIO_LIMIT = 2.0; // The theoretical boundary for enstrophy cascade

const GuardianStillness = {
  /**
   * checkStability
   * Monitors the Spectral Gap Ratio (R) between the 8x8 and 6x6 systems.
   * @param {Array<Float64Array>} L8 - Grounded 8x8 Laplacian.
   * @param {Array<Float64Array>} L6 - Grounded 6x6 Laplacian.
   * @returns {string} - Status of the system.
   */
  checkStability: (L8, L6) => {
    if (!L8 || !L6) {
      console.warn("Guardian: Missing system matrices.");
      return "VOID";
    }

    const eigen8 = LinAlg.eigvalsh(L8);
    const eigen6 = LinAlg.eigvalsh(L6);
    
    const lambda8 = eigen8[0];
    const lambda6 = eigen6[0];
    const R = lambda8 / lambda6;

    console.log(`Guardian Monitoring: R = ${R.toFixed(10)}`);

    if (R >= RESONANCE_RATIO_LIMIT) {
      console.error(`Guardian: RESONANCE BREACH. R = ${R.toFixed(4)} exceeds enstrophy floor (2.0).`);
      return "SHIELD_ACTIVATED";
    }

    return "STABLE";
  },

  reminder: () => {
    return "The ratio R < 2.0 is the shield against the cascade. The manifold holds.";
  }
};

module.exports = GuardianStillness;
