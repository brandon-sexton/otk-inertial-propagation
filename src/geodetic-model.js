import {
  getSphericalFromCartesian,
  getFixedFromInertial,
  getInertialFromFixed,
} from 'openspace-coordinates';
import {Vector} from 'vector-matrix-math';
/**
 * @description Oblate spheroid model of the Earth
 */
export class GeodeticModel {
  /**
   * @description Creates a new Geodetic Model
   * @example
   * const geodeticModel = new GeodeticModel();
   */
  constructor() {
    this.mu = 3.986004418e5;
    this.radius = 6378.1363;
    this.C = [
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
    ];
    this.S = [
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
    ];
  }

  /**
   * @description Gets the legendre polynomal for a given latitude
   * @param {number} phi - latitude in radians
   * @return {Array<Array<number>>}
   */
  getLegendrePolynomial(phi) {
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
   * @param {object} centralBody - Oblate model of Earth
   * @return {Vector} - Acceleration due to potential in an inertial frame.
   */
  getPotentialAcceleration(epoch, position) {
    const ecef = getFixedFromInertial(epoch, position);
    const sphrPos = getSphericalFromCartesian(ecef);
    const p = this.getLegendrePolynomial(sphrPos[2]);

    let partialR = 0.0;
    let partialPhi = 0.0;
    let partialLamb = 0.0;
    const recipR = 1.0 / position.magnitude();
    const muOverR = this.mu * recipR;
    const rOverR = this.radius * recipR;
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
          (this.C[n][m] * cLam + this.S[n][m] * sLam);
        partialPhi += rExponent * (p[n][m + 1] - m *
          Math.tan(sphrPos[2]) * p[n][m]) *
          (this.C[n][m] * cLam + this.S[n][m] * sLam);
        partialLamb += rExponent * m * p[n][m] *
          (this.S[n][m] * cLam - this.C[n][m] * sLam);
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
}
