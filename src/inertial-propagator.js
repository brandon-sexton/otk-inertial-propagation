import {SECONDS_TO_DAYS, DAYS_TO_SECONDS} from 'time-systems';
import {OrbitState} from './orbit-state';
import {ForceModelConfigs} from './force-model-configs';
import {getAcceleration} from './perturbations';

const ONE_SIXTH = 1/6;

/**
 * @class Propagator
 * @description Uses the Runge-Kutta 4th order method to propagate an orbit.
 */
export class InertialPropagator {
  MAX_STEP = 300 * SECONDS_TO_DAYS;
  /**
   * @constructor
   * @description Creates a new Propagator.
   * @param {OrbitState} y0 - Initial orbit state.
   * @example
   * const propagator = new Propagator(y0);
   */
  constructor(y0) {
    this.y0 = y0;
    this.forceModelConfigs = new ForceModelConfigs();
    this.stepSize = this.MAX_STEP;
  }

  /**
   * @method step
   * @description Steps the propagator forward in time.
   * @example
   * propagator.step();
   * @return {void}
   */
  step() {
    const h = this.stepSize;
    const dsecs = h * .5;
    const ddays = dsecs * SECONDS_TO_DAYS;

    const k1 = new OrbitStateDerivative(this.y0, this.forceModelConfigs);
    const t1 = this.y0.t0.plusDays(ddays);
    const r1 = this.y0.r0.plus(k1.v0.scale(dsecs));
    const v1 = this.y0.v0.plus(k1.a0.scale(dsecs));
    const y1 = new OrbitState(t1, r1, v1);

    const k2 = new OrbitStateDerivative(y1, this.forceModelConfigs);
    const t2 = t1.copy();
    const r2 = this.y0.r0.plus(k2.v0.scale(dsecs));
    const v2 = this.y0.v0.plus(k2.a0.scale(dsecs));
    const y2 = new OrbitState(t2, r2, v2);

    const k3 = new OrbitStateDerivative(y2, this.forceModelConfigs);
    const t3 = t1.plusDays(ddays);
    const r3 = this.y0.r0.plus(k3.v0.scale(h));
    const v3 = this.y0.v0.plus(k3.a0.scale(h));
    const y3 = new OrbitState(t3, r3, v3);

    const k4 = new OrbitStateDerivative(y3, this.forceModelConfigs);

    const dv = k1.v0.plus(
        k2.v0.scale(2)).plus(k3.v0.scale(2)).plus(k4.v0)
        .scale(ONE_SIXTH);
    const da = k1.a0.plus(
        k2.a0.scale(2)).plus(k3.a0.scale(2)).plus(k4.a0)
        .scale(ONE_SIXTH);

    const tf = t3;
    const rf = this.y0.r0.plus(dv.scale(h));
    const vf = this.y0.v0.plus(da.scale(h));
    this.y0 = new OrbitState(tf, rf, vf);
  }

  /**
   * @method stepToEpoch
   * @description Steps the propagator forward in time to a given epoch.
   * @param {Epoch} epoch - Epoch to step to.
   * @example
   * propagator.stepToEpoch(epoch);
   */
  stepToEpoch(epoch) {
    const totalDiff = (epoch.julianTT - this.y0.t0.julianTT) * DAYS_TO_SECONDS;
    const numSteps = Math.ceil(Math.abs(totalDiff) * this.MAX_STEP);
    const oldStepSize = this.stepSize;
    if (numSteps !== 0) {
      this.stepSize = totalDiff / numSteps;
    }
    for (let i = 0; i < numSteps; i++) {
      this.step();
    }
    if (epoch.toString() !== this.y0.t0.toString()) {
      this.stepToEpoch(epoch);
    }
    this.stepSize = oldStepSize;
  }
}

/**
 * @class OrbitStateDerivative
 * @description Represents the derivative of an orbit state.
 */
class OrbitStateDerivative {
  /**
   * @constructor
   * @param {OrbitState} state - Orbit state.
   * @param {ForceModelConfigs} forceModelConfigs - Force model configurations.
   * @example
   * const derivative = new OrbitStateDerivative(state, forceModelConfigs);
   */
  constructor(state, forceModelConfigs) {
    this.t0 = state.t0.copy();
    this.v0 = state.v0.copy();
    this.a0 = getAcceleration(state.t0, state.r0, forceModelConfigs);
  }
}
