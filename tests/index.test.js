import {Vector} from 'otk-linear-algebra';
import {Epoch} from 'otk-time-systems';
import {OrbitState} from '../src/orbit-state';
import {InertialPropagator} from '../src/inertial-propagator';

test('Test mefi propagation', () => {
  const epoch = new Epoch('2022-12-20T00:00:00.000Z');
  const position = new Vector(42164, 0, 0);
  const velocity = new Vector(0, 3.07375, 0);
  const state = new OrbitState(epoch, position, velocity);
  const endEpoch = new Epoch('2022-12-20T19:00:00.000Z');
  const propagator = new InertialPropagator(state);
  propagator.stepToEpoch(endEpoch, state);
  const propagatedState = propagator.y0;
  expect(propagatedState.r0[0]).toBeCloseTo(11701.163084, -1);
  expect(propagatedState.r0[1]).toBeCloseTo(-40487.016256, -1);
  expect(propagatedState.r0[2]).toBeCloseTo(-2.099302, -1);
  expect(propagatedState.v0[0]).toBeCloseTo(2.954853, 4);
  expect(propagatedState.v0[1]).toBeCloseTo(0.851923, 4);
  expect(propagatedState.v0[2]).toBeCloseTo(-0.000141, 4);
});

// const zero = new Vector(0, 0, 0);
// const epoch = new Epoch('2021-12-25T04:42:42.424Z');
// const state = new OrbitState(epoch, new Vector(10000, 40000, -5000), zero);
// const surface = new OrbitState(epoch, new Vector(6378, 0, 0), zero);

// test('Test acceleration from Sun', () => {
//   const a = getThirdBody(state.r0, epoch);
//   expect(a.magnitude()).toBeCloseTo(3.00696687920662e-9, 9);
// });

// test('Test acceleration from Earth', () => {
//   const a = getAccelerationFromEarth(surface.r0);
//   expect(a.magnitude()).toBeCloseTo(0.00979870641977297, 9);
// });

// test('Test acceleration from Moon', () => {
//   const a = getAccelerationFromMoon(state.r0, epoch);
//   expect(a.magnitude()).toBeCloseTo(3.4064259486775476e-9, 9);
// });
