import {getSunPosition, getMoonPosition} from 'celestial-body-positions';
import {GeodeticModel} from './geodetic-model';
/**
 * @class ForceModelConfigs
 */
export class ForceModelConfigs {
  DEFAULT_CENTRAL_BODY = new GeodeticModel();
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
}
