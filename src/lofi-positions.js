import {Vector} from 'vector-matrix-math';
import {degreesToRadians} from 'unit-conversions-js';
import {
  W_PLUS_W,
  COS_OBLIQUITY,
  SIN_OBLIQUITY,
  OBLIQUITY_OF_ECLIPTIC,
} from './constants';

/**
 * Get the sun's position in the J2000 ecliptic frame.
 * @param {Epoch} epoch - The epoch to calculate the sun's position for.
 * @return {Vector} - The sun's position in the J2000 ecliptic frame.
 */
export function getSunPosition(epoch) {
  const a = 0.0334133589;
  const b = 0.0003490659;
  const t = epoch.getJulianCenturiesPastJ2000();
  const ma = degreesToRadians(357.5256 + 35999.049 * t);
  const sma = Math.sin(ma);
  const cma = Math.cos(ma);
  const c2ma = Math.cos(2 * ma);
  const lam = W_PLUS_W + ma + a * sma + b * 2 * sma * cma;
  const r = (149.619 - 2.499 * cma - 0.021 * c2ma) * 1e6;
  const x = r * Math.cos(lam);
  const slam = Math.sin(lam);
  const y = r * slam * COS_OBLIQUITY;
  const z = r * slam * SIN_OBLIQUITY;
  return new Vector(x, y, z);
}

/**
 * Get the moon's position in the J2000 ecliptic frame.
 * @param {Epoch} epoch - The epoch to calculate the moon's position for.
 * @return {Vector} - The moon's position in the J2000 ecliptic frame.
 */
export function getMoonPosition(epoch) {
  const a0 = 0.109839480287776;
  const a1 = 0.0037322220424506;
  const a2 = 0.0221806468707036;
  const a3 = 0.0114797042944426;
  const a4 = 0.0032417852998172;
  const a5 = 0.0019977929354784;
  const a6 = 0.0010265025328168;
  const a7 = 0.0009997601892124;
  const a8 = 0.0009312755020032;
  const a9 = 0.0008002981275842;
  const a10 = 0.0007170000723212;
  const a11 = 0.000606017101918;
  const a12 = 0.0005332974922114;
  const a13 = 0.0002666487461057;

  const b0 = 0.0898933928087776;
  const b1 = 0.0019977929354784;
  const b2 = 0.0026218341701444;
  const b3 = 0.0025555590348276;
  const b4 = 0.0002131418144236;
  const b5 = 0.0001502389689262;
  const b6 = 0.000121203693847;
  const b7 = 0.000111229249161;
  const b8 = 0.000101254804475;
  const b9 = 0.000053137093141;

  const t = epoch.getJulianCenturiesPastJ2000();
  const l0 = degreesToRadians(218.31617 + 481267.88088 * t - 1.3972 * t);
  const l = degreesToRadians(134.96292 + 477198.86753 * t);
  const lp = degreesToRadians(357.52543 + 35999.04944 * t);
  const f = degreesToRadians(93.27283 + 483202.01873 * t);
  const d = degreesToRadians(297.85027 + 445267.11135 * t);
  const lam =
    l0 +
    a0 * Math.sin(l) +
    a1 * Math.sin(2 * l) -
    a2 * Math.sin(l - 2 * d) +
    a3 * Math.sin(2 * d) -
    a4 * Math.sin(lp) -
    a5 * Math.sin(2 * f) -
    a6 * Math.sin(2 * l - 2 * d) -
    a7 * Math.sin(l + lp - 2 * d) +
    a8 * Math.sin(l + 2 * d) -
    a9 * Math.sin(lp - 2 * d) +
    a10 * Math.sin(l - lp) -
    a11 * Math.sin(d) -
    a12 * Math.sin(l + lp) -
    a13 * Math.sin(2 * f - 2 * d);

  const beta =
    b0 * Math.sin(f + lam - l0 + b1 * Math.sin(2 * f + b2 * Math.sin(lp))) -
    b3 * Math.sin(f - 2 * d) +
    b4 * Math.sin(l + f - 2 * d) -
    b5 * Math.sin(-l + f - 2 * d) -
    b6 * Math.sin(-2 * l + f) -
    b7 * Math.sin(lp + f - 2 * d) +
    b8 * Math.sin(-l + f) +
    b9 * Math.sin(-lp + f - 2 * d);

  const r =
    385000 -
    20905 * Math.cos(l) -
    3699 * Math.cos(2 * d - l) -
    2956 * Math.cos(2 * d) -
    570 * Math.cos(2 * l) +
    246 * Math.cos(2 * l - 2 * d) -
    205 * Math.cos(lp - 2 * d) -
    171 * Math.cos(l + 2 * d) -
    152 * Math.cos(l + lp - 2 * d);

  const x = r * Math.cos(lam);
  const y = r * Math.sin(lam);
  const z = r * Math.sin(beta);

  const xAxis = new Vector(1, 0, 0);
  return new Vector(x, y, z).rotateAboutAxis(xAxis, OBLIQUITY_OF_ECLIPTIC);
}
