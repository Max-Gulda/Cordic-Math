#ifndef PID_CONTROLLER_H
#define PID_CONTROLLER_H

typedef struct {

	/* Controller gains */
	int Kp;
	int Ki;
	int Kd;

	/* Output limits */
	int limMin;
	int limMax;
	
	/* Integrator limits */
	int limMinInt;
	int limMaxInt;

	/* Controller "memory" */
	int proportional;
	int integrator;
	int prevError;			/* Required for integrator */
	int differentiator;
	int prevMeasurement;		/* Required for differentiator */

	/* Controller output */
	int out;

} PIDController;

void  PIDControllerFP_Init(PIDController *pid);
int PIDControllerFP_Update(PIDController *pid, int setpoint, int measurement, int deltaTime);
int millis(void);

#endif