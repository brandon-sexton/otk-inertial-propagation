import {getSunPosition, getMoonPosition} from 'otk-celestial-bodies';
import {
  getSphericalFromCartesian,
  getFixedFromInertial,
  getInertialFromFixed,
} from 'otk-coordinate-transforms';
import {Vector} from 'otk-linear-algebra';

/**
 * @class ForceModelConfigs
 */
export class ForceModel {
  DEFAULT_CENTRAL_BODY = {
    mu: 3.986004418e5,
    radius: 6378.1363,
    C: [
      [1],
      [0, 0],
      [
        -0.484165143790815e-3 / Math.sqrt(0.2),
        -0.206615509074176e-9 / Math.sqrt(0.6),
        0.243938357328313e-5 / Math.sqrt(2.4),
      ],
      [
        0.957161207093473e-6 / Math.sqrt(1.0 / 7.0),
        0.203046201047864e-5 / Math.sqrt(6.0 / 7.0),
        0.904787894809528e-6 / Math.sqrt(60.0 / 7.0),
        0.721321757121568e-6 / Math.sqrt(360.0 / 7.0),
      ],
      [
        0.539965866638991e-6 / Math.sqrt(24.0 / 216.0),
        -0.536157389388867e-6 / Math.sqrt(10.0 / 9.0),
        0.350501623962649e-6 / Math.sqrt(20.0),
        0.990856766672321e-6 / Math.sqrt(280.0),
        -0.188519633023033e-6 / Math.sqrt(2240.0),
      ],
    ],
    S: [
      [0],
      [0, 0],
      [
        0.0,
        0.138441389137979e-8 / Math.sqrt(0.6),
        -0.140027370385934e-5 / Math.sqrt(2.4),
      ],
      [
        0.0,
        0.248200415856872e-6 / Math.sqrt(6.0 / 7.0),
        -0.619005475177618e-6 / Math.sqrt(60.0 / 7.0),
        0.141434926192941e-5 / Math.sqrt(360.0 / 7.0),
      ],
      [
        0.0,
        -0.473567346518086e-6 / Math.sqrt(10.0 / 9.0),
        0.662480026275829e-6 / Math.sqrt(20.0),
        -0.200956723567452e-6 / Math.sqrt(280.0),
        0.308803882149194e-6 / Math.sqrt(2240.0),
      ],
    ],
    getGravityAcceleration: getGravityAcceleration,
  };
  DEFAULT_SRP_AREA = 17;
  DEFAULT_SRP_CR = 1.8;
  DEFAULT_DRAG_AREA = 17;
  DEFAULT_DRAG_CD = 2.2;
  DEFAULT_THIRD_BODIES = [
    {
      mu: 4902.800066,
      position: getMoonPosition,
    },
    {
      mu: 1.327124400419e11,
      position: getSunPosition,
    },
  ];
  /**
   * @constructor
   */
  constructor() {
    this.centralBody = this.DEFAULT_CENTRAL_BODY;
    this.thirdBodies = this.DEFAULT_THIRD_BODIES;
    this.srpArea = this.DEFAULT_SRP_AREA;
    this.srpCr = this.DEFAULT_SRP_CR;
    this.dragArea = this.DEFAULT_DRAG_AREA;
    this.dragCd = this.DEFAULT_DRAG_CD;
  }

  /**
   * @param {Epoch} epoch - Epoch of the satellite state.
   * @param {Vector} position - Position of the satellite in the inertial frame.
   * @return {Vector} - Net acceleration in the inertial frame.
   */
  getAcceleration(epoch, position) {
    let a = getCentralBodyAcceleration(position, this.centralBody.mu);
    for (const body of this.thirdBodies) {
      const bodyPosition = body.position(epoch);
      const mu = body.mu;
      a = a.plus(getThirdBodyAcceleration(position, mu, bodyPosition));
    }
    a = a.plus(
        this.centralBody.getGravityAcceleration(
            epoch,
            position,
            this.centralBody,
        ),
    );
    return a;
  }
}

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
 * @description Gets the legendre polynomal for a given latitude
 * @param {number} phi - latitude in radians
 * @return {Array<Array<number>>}
 */
function getLegendrePolynomial(phi) {
  const cosPhi = Math.cos(phi);
  const sinPhi = Math.sin(phi);
  const cosPhiSquared = cosPhi * cosPhi;
  const sinPhiSquared = sinPhi * sinPhi;

  const p = [
    [1.0, 0.0],
    [sinPhi, cosPhi, 0.0],
    [
      (3.0 * sinPhiSquared - 1.0) * 0.5,
      3.0 * sinPhi * cosPhi,
      3.0 * cosPhiSquared,
      0.0,
    ],
    [
      sinPhi * (5.0 * sinPhiSquared - 3.0) * 0.5,
      (15.0 * sinPhiSquared - 3.0) * cosPhi * 0.5,
      15.0 * sinPhi * cosPhiSquared,
      15.0 * cosPhiSquared * cosPhi,
      0.0,
    ],
    [
      0.125 * (35 * sinPhiSquared * sinPhiSquared - 30 * sinPhiSquared + 3.0),
      2.5 * (7.0 * sinPhiSquared * sinPhi - 3.0 * sinPhi) * cosPhi,
      (7.0 * sinPhiSquared - 1.0) * cosPhiSquared * 7.5,
      105.0 * cosPhi * cosPhiSquared * sinPhi,
      105.0 * cosPhiSquared * cosPhiSquared,
      0.0,
    ],
  ];

  return p;
}

/**
 * @description Gets the acceleration due to the body potential
 * @param {Epoch} epoch - Epoch of the satellite state.
 * @param {Vector} position - Position of the satellite in the inertial frame.
 * @param {object} centralBody
 * @return {Vector} - Acceleration due to potential in an inertial frame.
 */
function getGravityAcceleration(epoch, position, centralBody) {
  const ecef = getFixedFromInertial(epoch, position);
  const sphrPos = getSphericalFromCartesian(ecef);
  const p = getLegendrePolynomial(sphrPos[2]);

  let partialR = 0.0;
  let partialPhi = 0.0;
  let partialLamb = 0.0;
  const recipR = 1.0 / position.magnitude();
  const muOverR = centralBody.mu * recipR;
  const rOverR = centralBody.radius * recipR;
  let rExponent = 0.0;
  let cLam = 0.0;
  let sLam = 0.0;
  const recipRoot = 1.0 / Math.sqrt(ecef[0] * ecef[0] + ecef[1] * ecef[1]);
  const rzOverRoot = ecef[2] * recipR * recipR * recipRoot;

  let m = 0;
  let n = 2;
  while (n < 5) {
    m = 0;
    rExponent = Math.pow(rOverR, n);
    while (m <= n) {
      cLam = Math.cos(m * sphrPos[1]);
      sLam = Math.sin(m * sphrPos[1]);
      partialR += rExponent * (n + 1) * p[n][m] *
        (centralBody.C[n][m] * cLam + centralBody.S[n][m] * sLam);
      partialPhi += rExponent * (p[n][m + 1] - m *
        Math.tan(sphrPos[2]) * p[n][m]) *
        (centralBody.C[n][m] * cLam + centralBody.S[n][m] * sLam);
      partialLamb += rExponent * m * p[n][m] *
        (centralBody.S[n][m] * cLam - centralBody.C[n][m] * sLam);
      m++;
    }
    n++;
  }

  partialR *= -recipR * muOverR;
  partialPhi *= muOverR;
  partialLamb *= muOverR;

  return getInertialFromFixed(
      epoch,
      new Vector(
          (recipR * partialR - rzOverRoot * partialPhi) *
            ecef[0] - recipRoot * recipRoot * partialLamb * ecef[1],
          (recipR * partialR - rzOverRoot * partialPhi) *
            ecef[1] + recipRoot * recipRoot * partialLamb * ecef[0],
          recipR * partialR * ecef[2] + 1.0 / recipRoot *
            recipR * recipR * partialPhi,
      ),
  );
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
