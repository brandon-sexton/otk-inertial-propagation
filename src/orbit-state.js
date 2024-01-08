/**
 * @class OrbitState
 * @description Represents time, position, and velocity of an orbiting object
 */
export class OrbitState {
  /**
   * @constructor
   * @param {Epoch} t0 - The initial epoch.
   * @param {Vector} r0 - The initial position vector.
   * @param {Vector} v0 - The initial velocity vector.
   */
  constructor(t0, r0, v0) {
    this.t0 = t0.copy();
    this.r0 = r0.copy();
    this.v0 = v0.copy();
  }

  /**
   * @method copy
   * @description Returns a copy of this state.
   * @return {OrbitState} A copy of this state.
   */
  copy() {
    return new OrbitState(this.t0, this.r0, this.v0);
  }

  /**
   * Returns a string representing this state.
   *
   * @return {string} A string representing this state.
   */
  toString() {
    return `OrbitState(${this.t0}, ${this.r0}, ${this.v0})`;
  }
}
