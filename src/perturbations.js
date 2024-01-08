/**
 * @param {Vector} position - Position of the satellite in the inertial frame.
 * @param {number} mu - Gravitational parameter of the central body.
 * @return {Vector} - Acceleration due to the central body in an inertial frame.
 */
function getCentralBodyAcceleration(position, mu) {
  const r = position.magnitude();
  return position.scale(-mu / (r * r * r));
}

/**
 * @param {Vector} position - Position of the satellite in the inertial frame.
 * @param {number} mu - Gravitational parameter of the third body.
 * @param {Vector} bodyPosition - Inertial position of the third body.
 * @return {Vector} - Inertial acceleration due to the third body
 */
function getThirdBodyAcceleration(position, mu, bodyPosition) {
  const s = bodyPosition;
  const r = s.minus(position);
  const sMag = s.magnitude();
  const rMag = r.magnitude();
  const v1 = r.scale(1 / (rMag * rMag * rMag));
  const v2 = s.scale(1 / (sMag * sMag * sMag));
  return v1.minus(v2).scale(mu);
}

/**
 * @param {Epoch} epoch - Epoch of the satellite state.
 * @param {Vector} position - Position of the satellite in the inertial frame.
 * @param {ForceModelConfigs} configs - Force model configurations.
 * @return {Vector} - Acceleration due to the force model in an inertial frame.
 */
export function getAcceleration(epoch, position, configs) {
  let a = getCentralBodyAcceleration(position, configs.centralBody.mu);
  for (const body of configs.thirdBodies) {
    const bodyPosition = body.position(epoch);
    const mu = body.mu;
    a = a.plus(getThirdBodyAcceleration(position, mu, bodyPosition));
  }
  a = a.plus(
      configs.centralBody.getPotentialAcceleration(
          epoch,
          position,
          configs.centralBody,
      ),
  );
  return a;
}

/**
 * @return {Object} - Unit test exports.
 */
export function unitTestExports() {
  return {
    getCentralBodyAcceleration,
    getThirdBodyAcceleration,
  };
}
