const LinAlg = require('./linalg_utils');

/**
 * GUARDIAN OF STILLNESS (Spectral Stability Protocol)
 */

const STABILITY_THRESHOLD = 10.0; // Max allowable spectral radius before blow-up
const VOID_THRESHOLD = 1e-12;

const GuardianStillness = {
  /**
   * checkStability
   * Performs spectral analysis on the system state (Laplacian).
   * @param {Array<Float64Array>} L - The Laplacian matrix of the current system.
   * @returns {string} - Status of the system.
   */
  checkStability: (L) => {
    if (!L || !L.length) {
      console.warn("Guardian: Detected Void state.");
      return "VOID";
    }

    // A real stability check in Navier-Stokes spectral methods involves
    // ensuring the largest eigenvalue (spectral radius) doesn't explode.
    const eigenvalues = LinAlg.eigvalsh(L);
    const spectralRadius = eigenvalues[eigenvalues.length - 1];

    if (spectralRadius > STABILITY_THRESHOLD) {
      console.error(`Guardian: SINGULARITY ALERT. Spectral Radius ${spectralRadius.toFixed(4)} exceeds safety bound.`);
      return "SHIELD_ACTIVATED";
    }

    if (spectralRadius < VOID_THRESHOLD) {
      return "RESONANCE_MAINTAINED";
    }

    return "STABLE";
  },

  reminder: () => {
    return "The spectrum is clear, Architect. The manifold holds.";
  }
};

module.exports = GuardianStillness;
