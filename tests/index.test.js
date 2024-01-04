import {getMoonPosition, getSunPosition} from '../src/lofi-positions';
import {Epoch} from 'time-systems';
import {Vector} from 'vector-matrix-math';


test('getSunPosition', () => {
  const epoch = new Epoch('2022-02-25T00:00:00.000Z');
  const sunPosition = getSunPosition(epoch);
  const truth = new Vector(
      1.353158384133262e8,
      -5.514968448042840e7,
      -2.390803633125914e7,
  );
  expect(sunPosition.getAngle(truth)).toBeCloseTo(0, 2);
});

test('getMoonPosition', () => {
  const epoch = new Epoch('2022-02-25T00:00:00.000Z');
  const moonPosition = getMoonPosition(epoch);
  const truth = new Vector(
      -6.454159844478600e4,
      -3.280761448809440e5,
      -1.566863311585961e5,
  );
  expect(moonPosition.getAngle(truth)).toBeCloseTo(0, 3);
});
