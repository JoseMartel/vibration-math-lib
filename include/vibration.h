#ifndef vibration_h /* Prevent loading library twice */
#define vibration_h

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

/* Custom constants */
#define MIN_FREQ 5 // -> Related to spectral resolution freq/samples
#define MIN_POWER_ON_VIBRATION 40
#define FFT_FORWARD 0x01
#define FFT_REVERSE 0x00
#define MAX_PEAKS 5

/* Constant for gravity */
#define CONST_GRAVITY 9.80665

/* Windowing type */
#define FFT_WIN_TYP_RECTANGLE 			0x00 /* rectangle (Box car) */
#define FFT_WIN_TYP_HAMMING 			0x01 /* hamming */
#define FFT_WIN_TYP_HANN 				0x02 /* hann */
#define FFT_WIN_TYP_TRIANGLE 			0x03 /* triangle (Bartlett) */
#define FFT_WIN_TYP_NUTTALL 			0x04 /* nuttall */
#define FFT_WIN_TYP_BLACKMAN 			0x05 /* blackman */
#define FFT_WIN_TYP_BLACKMAN_NUTTALL 	0x06 /* blackman nuttall */
#define FFT_WIN_TYP_BLACKMAN_HARRIS 	0x07 /* blackman harris*/
#define FFT_WIN_TYP_FLT_TOP 			0x08 /* flat top */
#define FFT_WIN_TYP_WELCH 				0x09 /* welch */

/* EXTI - DATA READY PIN */
#define LIS3DH_DRDY_PIN                 IOID_22

/* Anomalie */
#define NO_ANOMALIE     0x00
#define LOOSENESS       0x01
#define MISALIGNMENT    0x02
#define UNBALANCE       0x04
#define VIBRATION       0x08
#define RPM             0x10
#define BEHAVIOR        0x20

/* Typedef declaration */
typedef struct {
    float amp;
    float freq;
} Peaks_t;

typedef struct {
    float       aRMS; 
    float       vRMS;
    Peaks_t     Peak1;
    Peaks_t     Peak2;
    Peaks_t     Peak3;
    uint16_t    Speed;
    uint8_t     Anomalie;
    uint8_t     MotorState;
    uint8_t     SensorState;    
} Vibration_t;

typedef struct {
    uint16_t    maxVibration;
    uint16_t    nomRPM;
    uint16_t    nomPower;
    uint8_t     ThresholdSelect;
    uint8_t     Foundation;
    uint8_t     Motor;
    uint8_t     Poles;
    uint8_t     Freq;
} Parameters_t;

// Utils
void    getVibration(Vibration_t *_vibration);
void    saveVibration(Vibration_t  *_vibration);
void    getParameters(Parameters_t *_parameters);
void    saveParameters(Parameters_t *_parameters);
void    ParameterDefault(Parameters_t  *_parameters);
uint8_t	FFT_LibRevision(void);

// FFT Functions
void 	FFT_Init(float *vReal, float *vImag, uint16_t samples, float samplingFrequency);
void 	FFT_Windowing(uint8_t windowType, uint8_t dir);
void 	FFT_Compute(uint8_t dir);
void 	FFT_ComplexToMagnitude(void);
void 	FFT_MajorPeak(float* mag_out, float* freq_out, float magFact);
void 	FFT_MajorPeaks(Peaks_t* peaks, float magFact, int numPeaks);

// Signal processing
void 	removeOffset(float *array, size_t array_len);
void 	Remove3Offset(float *x, float *y, float *z, size_t len);
void	RMS_Value(float *data, uint16_t samples, float *RMSValue);
void	Velocity(float *a, float *v, size_t len, float dt);
void	rungeKutta(float *a, float *v, size_t len, int freq);
void	simpson13(float *a, float *v, size_t len, int freq);
void	resVector(float *x, float *y, float *z, uint16_t len);
void    getResultant(float *res, float *x, float *y, float *z, uint16_t len);
void    high_pass_filter_execute(float *input, int length);
void    detectAnomaly(float newMean, float normalMean, uint8_t anomaly);
void    mean(float *signal, uint16_t len ,float *mean);
void    deviation(float *signal, uint16_t len, float mean, float *deviation);

// Utils
void    vibAnalysis(Vibration_t *vibration, Parameters_t *parameters);
uint8_t getStatus(void);
void    setStatus(uint8_t sts);
uint8_t getAnBehavior( void );
void    setAnBehavior( uint8_t value );
uint8_t getEnLearning( void );
void    setEnLearning( uint8_t value );

#endif


