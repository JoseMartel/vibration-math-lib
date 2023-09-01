/**
 * @file vibration.c
 * @brief Utils functions to vibration analysis
 * @author Jose Martel & ChatGPT
 * @date 30/06/2023
 */

#include <stdio.h>
#include "vibration.h"
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

/* Mathematial constants */
#define twoPi 	6.28318531
#define fourPi 	12.56637061
#define sixPi 	18.84955593

#define FFT_LIB_REV 0x10

/* Private Variables */
static uint8_t	mb_status = 1; 
static uint16_t _samples;
static float 	_samplingFrequency;
static float*	_vReal;
static float*	_vImag;
static uint8_t 	_power;
static uint8_t	anBehavior;
static uint8_t	enLearning = 0;

static Vibration_t vibration;
static Parameters_t parameters;

/* Private Functions and Macros*/
static void 	Swap(float *x, float *y);
static uint8_t 	Exponent(uint16_t value);
#define sq(x) 	((x)*(x))

void getVibration(Vibration_t *_vibration)
{
	_vibration->aRMS = vibration.aRMS;
	_vibration->vRMS = vibration.vRMS;
	_vibration->Peak1 = vibration.Peak1;
	_vibration->Peak2 = vibration.Peak2;
	_vibration->Peak3 = vibration.Peak3;
	_vibration->Anomalie = vibration.Anomalie;
	_vibration->Speed = vibration.Speed;
	_vibration->MotorState = vibration.MotorState;
	_vibration->SensorState = vibration.SensorState;

}

void saveVibration(Vibration_t  *_vibration)
{
	vibration.aRMS = _vibration->aRMS;
	vibration.vRMS = _vibration->vRMS;
	vibration.Peak1 = _vibration->Peak1;
	vibration.Peak2 = _vibration->Peak2;
	vibration.Peak3 = _vibration->Peak3;
	vibration.Anomalie = _vibration->Anomalie;
	vibration.Speed = _vibration->Speed;
	vibration.MotorState = _vibration->MotorState;
	vibration.SensorState = _vibration->SensorState;
}

void    getParameters(Parameters_t *_parameters)
{
	_parameters->Foundation = parameters.Foundation;
	_parameters->Freq = parameters.Freq;
	_parameters->maxVibration = parameters.maxVibration;
	_parameters->Motor = parameters.Motor;
	_parameters->nomPower = parameters.nomPower;
	_parameters->nomRPM = parameters.nomRPM;
	_parameters->Poles = parameters.Poles;
	_parameters->ThresholdSelect = parameters.ThresholdSelect;
}

void    saveParameters(Parameters_t *_parameters)
{
	parameters.Foundation = _parameters->Foundation;
	parameters.Freq = _parameters->Freq;
	parameters.maxVibration = _parameters->maxVibration;
	parameters.Motor = _parameters->Motor;
	parameters.nomPower = _parameters->nomPower;
	parameters.nomRPM = _parameters->nomRPM;
	parameters.Poles = _parameters->Poles;
	parameters.ThresholdSelect = _parameters->ThresholdSelect;
}

void ParameterDefault(Parameters_t  *_parameters)
{
	_parameters->Foundation = 0;
	_parameters->Freq = 60;
	_parameters->maxVibration = 1;
	_parameters->Motor = 0;
	_parameters->nomPower = 500;
	_parameters->nomRPM = 3200;
	_parameters->Poles = 2;
	_parameters->ThresholdSelect = 0;
}

void FFT_Init(float *vReal, float *vImag, uint16_t samples, float samplingFrequency)
{
	_vReal = vReal;
	_vImag = vImag;
	_samples = samples;
	_samplingFrequency = samplingFrequency;
	_power = Exponent(samples);
}

void FFT_Compute(uint8_t dir)
{// Computes in-place complex-to-complex FFT /
	// Reverse bits /
	uint16_t j = 0;
	for (uint16_t i = 0; i < (_samples - 1); i++) {
		if (i < j) {
			Swap(&_vReal[i], &_vReal[j]);
			if(dir==FFT_REVERSE)
				Swap(&_vImag[i], &_vImag[j]);
		}
		uint16_t k = (_samples >> 1);
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
	// Compute the FFT  /
	float c1 = -1.0;
	float c2 = 0.0;
	uint16_t l2 = 1;
	for (uint8_t l = 0; (l < _power); l++) {
		uint16_t l1 = l2;
		l2 <<= 1;
		float u1 = 1.0;
		float u2 = 0.0;
		for (j = 0; j < l1; j++) {
			 for (uint16_t i = j; i < _samples; i += l2) {
					uint16_t i1 = i + l1;
					float t1 = u1 * _vReal[i1] - u2 * _vImag[i1];
					float t2 = u1 * _vImag[i1] + u2 * _vReal[i1];
					_vReal[i1] = _vReal[i] - t1;
					_vImag[i1] = _vImag[i] - t2;
					_vReal[i] += t1;
					_vImag[i] += t2;
			 }
			 float z = ((u1 * c1) - (u2 * c2));
			 u2 = ((u1 * c2) + (u2 * c1));
			 u1 = z;
		}
		c2 = sqrt((1.0 - c1) / 2.0);
		c1 = sqrt((1.0 + c1) / 2.0);
		if (dir == FFT_FORWARD) {
			c2 = -c2;
		}
	}
	// Scaling for reverse transform /
	if (dir != FFT_FORWARD) {
		for (uint16_t i = 0; i < _samples; i++) {
			 _vReal[i] /= _samples;
			 _vImag[i] /= _samples;
		}
	}
}

void FFT_ComplexToMagnitude(void)
{ // vM is half the size of vReal and vImag
	for (uint16_t i = 0; i < _samples; i++) {
		_vReal[i] = sqrt(sq(_vReal[i]) + sq(_vImag[i]));
	}
}

void FFT_Windowing(uint8_t windowType, uint8_t dir)
{// Weighing factors are computed once before multiple use of FFT
// The weighing function is symetric; half the weighs are recorded
	float samplesMinusOne = ((float)(_samples) - 1.0);
	for (uint16_t i = 0; i < (_samples >> 1); i++) {
		float indexMinusOne = (float)(i);
		float ratio = (indexMinusOne / samplesMinusOne);
		float weighingFactor = 1.0;
		// Compute and record weighting factor
		switch (windowType) {
		case FFT_WIN_TYP_RECTANGLE: // rectangle (box car)
			weighingFactor = 1.0;
			break;
		case FFT_WIN_TYP_HAMMING: // hamming
			weighingFactor = 0.54 - (0.46 * cos(twoPi * ratio));
			break;
		case FFT_WIN_TYP_HANN: // hann
			weighingFactor = 0.54 * (1.0 - cos(twoPi * ratio));
			break;
		case FFT_WIN_TYP_TRIANGLE: // triangle (Bartlett)
			weighingFactor = 1.0 - ((2.0 * fabs(indexMinusOne - (samplesMinusOne / 2.0))) / samplesMinusOne);
			break;
		case FFT_WIN_TYP_NUTTALL: // nuttall
			weighingFactor = 0.355768 - (0.487396 * (cos(twoPi * ratio))) + (0.144232 * (cos(fourPi * ratio))) - (0.012604 * (cos(sixPi * ratio)));
			break;
		case FFT_WIN_TYP_BLACKMAN: // blackman
			weighingFactor = 0.42323 - (0.49755 * (cos(twoPi * ratio))) + (0.07922 * (cos(fourPi * ratio)));
			break;
		case FFT_WIN_TYP_BLACKMAN_NUTTALL: // blackman nuttall
			weighingFactor = 0.3635819 - (0.4891775 * (cos(twoPi * ratio))) + (0.1365995 * (cos(fourPi * ratio))) - (0.0106411 * (cos(sixPi * ratio)));
			break;
		case FFT_WIN_TYP_BLACKMAN_HARRIS: // blackman harris
			weighingFactor = 0.35875 - (0.48829 * (cos(twoPi * ratio))) + (0.14128 * (cos(fourPi * ratio))) - (0.01168 * (cos(sixPi * ratio)));
			break;
		case FFT_WIN_TYP_FLT_TOP: // flat top
			weighingFactor = 0.2810639 - (0.5208972 * cos(twoPi * ratio)) + (0.1980399 * cos(fourPi * ratio));
			break;
		case FFT_WIN_TYP_WELCH: // welch
			weighingFactor = 1.0 - sq((indexMinusOne - samplesMinusOne / 2.0) / (samplesMinusOne / 2.0));
			break;
		}
		if (dir == FFT_FORWARD) {
			_vReal[i] *= weighingFactor;
			_vReal[_samples - (i + 1)] *= weighingFactor;
		}
		else {
			_vReal[i] /= weighingFactor;
			_vReal[_samples - (i + 1)] /= weighingFactor;
		}
	}
}

void FFT_MajorPeak(float* mag_out, float* freq_out, float magFact)
{
	float maxY = 0;
	uint16_t IndexOfMaxY = 0;

	for (uint16_t i = 1; i < ((_samples >> 1) + 1); i++) {
		if ((_vReal[i-1] < _vReal[i]) && (_vReal[i] > _vReal[i+1])) {
			if (_vReal[i] > maxY) {
				maxY = _vReal[i];
				IndexOfMaxY = i;
			}
		}
	}
	float delta = 0.5 * ((_vReal[IndexOfMaxY-1] - _vReal[IndexOfMaxY+1]) / (_vReal[IndexOfMaxY-1] - (2.0 * _vReal[IndexOfMaxY]) + _vReal[IndexOfMaxY+1]));
	float interpolatedX = ((IndexOfMaxY + delta)  * _samplingFrequency) / (_samples-1);
	if(IndexOfMaxY==(_samples >> 1))
		interpolatedX = ((IndexOfMaxY + delta)  * _samplingFrequency) / (_samples);

	*mag_out = _vReal[IndexOfMaxY]/magFact;
	*freq_out = interpolatedX;
}

void FFT_MajorPeaks(Peaks_t* peaks, float magFact, int numPeaks)
{	
	/* Delete low frequencies components */
	_vReal[0] = 0; /* Amplitude at 0Hz */
	_vReal[1] = 0; /* Amplitude at 1.58 Hz */
	_vReal[2] = 0; /* Amplitude at 3.17 Hz */
	_vReal[3] = 0; /* Amplitude at 4.75 Hz */
	_vReal[4] = 0; /* Amplitude at 6.34 Hz */
	_vReal[5] = 0; /* Amplitude at 7.86 Hz */
	_vReal[6] = 0; /* Amplitude at 9.34 Hz */

    // Inicializar los arreglos de salida
    for (int i = 0; i < numPeaks; i++) {
        peaks[i].amp = 0.0;
        peaks[i].freq = 0.0;
    }

    // Buscar los picos de la FFT
    uint16_t pos_out[numPeaks];
    for (int k = 0; k < numPeaks; k++) {
        float maxY = 0;
        uint16_t IndexOfMaxY = 0;
        for (uint16_t i = 1; i < ((_samples >> 1) + 1); i++) {
            if ((_vReal[i-1] < _vReal[i]) && (_vReal[i] > _vReal[i+1])) {
                bool PeakFound = true;
                for (int j = 0; j < k; j++) {
                    if (i == pos_out[j]) {
                        PeakFound = false;
                        break;
                    }
                }
                if (PeakFound) {
                    if ((_vReal[i] > maxY)) {
                        maxY = _vReal[i];
                        IndexOfMaxY = i;
                    }
                }
            }
        }
        // Interpolar la frecuencia del pico
        float delta = 0.5 * ((_vReal[IndexOfMaxY-1] - _vReal[IndexOfMaxY+1]) / (_vReal[IndexOfMaxY-1] - (2.0 * _vReal[IndexOfMaxY]) + _vReal[IndexOfMaxY+1]));
        float interpolatedX = ((IndexOfMaxY + delta)  * _samplingFrequency) / (_samples-1);
        if(IndexOfMaxY==(_samples >> 1))
            interpolatedX = ((IndexOfMaxY + delta)  * _samplingFrequency) / (_samples);

        // Guardar la magnitud y la frecuencia del pico en los arreglos de salida
		pos_out[k] = IndexOfMaxY;

        peaks[k].amp = _vReal[IndexOfMaxY]/magFact*0.031883;;
		if((interpolatedX < 1) || isnan(interpolatedX))
			peaks[k].freq = 0.0;
		else
			peaks[k].freq = interpolatedX;
    }
}


void removeOffset(float *array, size_t len)
{
	float mean = 0;
	float naem = 0;

	/* Calculate mean of acceleration vector */
	for (size_t i = 0; i < len; i++)
	{
		mean += array[i];
	}
	mean = mean / len;

	/* Subtract the mean from array */
	for (size_t j = 0; j < len; j++)
	{	
		array[j] = array[j] - mean;
		naem += fabs(array[j]); 
	}
	naem = naem / len;

	printf("Mean : %.2f \n", naem);

	/* If new mean is less than 0.035g it consider noise and motor off */
	if(naem < 35.0)
	{
		for(size_t i = 0; i < len; i++)
		{
			array[i] = 0;
		}
	}
}

void Remove3Offset(float *x, float *y, float *z, size_t len)
{
    float xmean = 0, ymean = 0, zmean = 0;

    // Calculate the mean of acceleration data
    for (size_t i = 0; i < len; i++)
    {
        xmean += x[i];
        ymean += y[i];
        zmean += z[i];
    }
    
    xmean /= len;
    ymean /= len;
    zmean /= len;

    // Subtract the mean from vibration data in 3 axis
    for (size_t i = 0; i < len; i++)
    {
        x[i] -= xmean;
        y[i] -= ymean;
        z[i] -= zmean;

        // Delete static noise values
        if (fabs(x[i]) < 25.0f)
            x[i] = 0.0f;

        if (fabs(y[i]) < 25.0f)
            y[i] = 0.0f;

        if (fabs(z[i]) < 25.0f)
            z[i] = 0.0f;
    }
}



void RMS_Value(float *data, uint16_t n, float *RMSValue)
{
    float sum = 0.0, rms = 0.0;
	uint16_t i;

    // Calcular la suma de los cuadrados
    for (i = 0; i < n; i++)
	{
        sum += data[i] * data[i];
    }

    // Calcular el valor rms
	rms = sqrt(sum / n);

    *RMSValue = rms;
}

void Velocity(float *a, float *v, size_t len, float dt)
{
	// Initial Conditions
	v[0] = 0.0;

    for (int i = 1; i < (len-1); i++) 
	{
        v[i] = (a[i-1]*CONST_GRAVITY + a[i]*CONST_GRAVITY) / 2 * dt ;
    }

	v[len-1] = 0;
}

void rungeKutta(float *a, float *v, size_t len, int freq)
{
    double delta_t = 1.0 / freq;
    double k1, k2, k3, k4;
    v[0] = 0.0; // Velocidad inicial asumida como cero

    for (int i = 1; i < len; i++) {
        // Cálculo de los coeficientes k1, k2, k3 y k4
        k1 = a[i - 1];
        k2 = a[i - 1] + 0.5 * delta_t * k1;
        k3 = a[i - 1] + 0.5 * delta_t * k2;
        k4 = a[i];

        // Cálculo de la velocidad en el tiempo actual utilizando los coeficientes
        v[i] = v[i - 1] + (delta_t / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
    }
}

void simpson13(float *a, float *v, size_t len, int freq)
{
    float dt = 1.0 / freq;

	// Initial Conditions
	v[0] = 0.0;

    // Integración numérica de Simpson 1/3
    for (int i = 0; i < (len-3); i++) 
	{
		v[i] = (2 * dt / 3) * (a[i]*CONST_GRAVITY + 4 * a[i + 1]*CONST_GRAVITY + a[i + 2]*CONST_GRAVITY);
    }

	v[len-1] = 0;
	v[len-2] = 0;
	v[len-3] = 0;
}

void resVector(float *x, float *y, float *z, uint16_t len){

	// for (uint16_t i = 800; i < 821; i++) 
	// {
	// 	printf("[$] x: %.2f | y: %.2f | z: %.2f \n", x[i], y[i], z[i]);
	// }

	for (uint16_t i = 0; i < len; i++) 
	{
		// printf("[$] x: %.2f | y: %.2f | z: %.2f |", x[i], y[i], z[i]);
		// x[i] = sqrt((x[i]*x[i]) + (y[i]*y[i]) ); //+ (z[i]*z[i]));
		// printf("[$] R: %.2f\n", x[i]);
		/* Vibracion en el plano radial */
		x[i] = (x[i] + y[i])/2;
	}
}

void getResultant(float *res, float *x, float *y, float *z, uint16_t len){

	for (uint16_t i = 0; i < len; i++) 
	{
		res[i] = sqrt((x[i]*x[i]) + (y[i]*y[i]) + (z[i]*z[i]));
		// printf("[$] R: %.2f | y: %.2f | z: %.2f \n", x[i], y[i], z[i]);
	}

	// for (uint16_t i = 753; i < 821; i++) 
	// {
	// 	printf("[$] R: %.2f | y: %.2f | z: %.2f \n", x[i], y[i], z[i]);
	// }

}

// ============================================

void vibAnalysis(Vibration_t *vibration, Parameters_t *parameters)
{
	if(parameters->ThresholdSelect == 1){
		/* Determine vibration threshold using nominal 
			values of motor based a ISO vibration table */

		float lim1ISO = 0;

		if(parameters->nomPower > 15 && parameters->nomPower < 300){
			if(parameters->Foundation == 0)
				lim1ISO = 4.5;
			else
				lim1ISO = 2.8;
		} else if (parameters->nomPower > 300 && parameters->nomPower < 50000){
			if(parameters->Foundation == 0)
				lim1ISO = 7.1;
			else
				lim1ISO = 4.5;
		}

		/* Detect Anomaly based on velocity vibration */
		if(vibration->vRMS > lim1ISO)
			vibration->Anomalie = (vibration->Anomalie | VIBRATION);
		else
		 	vibration->Anomalie = (vibration->Anomalie & VIBRATION);
	} else {

		/* Detect Anomaly based on velocity vibration */
		if (vibration->vRMS > ((float)parameters->maxVibration))
			vibration->Anomalie = (vibration->Anomalie | VIBRATION);
		else
		 	vibration->Anomalie = (vibration->Anomalie & ~VIBRATION);
	}

	/* Detect motor state */
	if(vibration->aRMS > MIN_POWER_ON_VIBRATION)
		vibration->MotorState = 1;
	else
		vibration->MotorState = 0;

	/* Identify main frequency */
	float mainFreq;
	float D1 = fabs(parameters->nomRPM / 60.0 - vibration->Peak1.freq); 
	float D2 = fabs(parameters->nomRPM / 60.0 - vibration->Peak2.freq); 
	float D3 = fabs(parameters->nomRPM / 60.0 - vibration->Peak3.freq);

	if ((D1 <= D2) && (D1 <= D3)) 
		mainFreq = vibration->Peak1.freq;
	else if ((D2 <= D1) && (D2 <= D3))
		mainFreq = vibration->Peak2.freq;
	else 
		mainFreq = vibration->Peak3.freq;

	/* Detect Anomaly Looseness */
	if(vibration->Peak1.freq > mainFreq)
		vibration->Anomalie = (vibration->Anomalie | LOOSENESS);
	else
		vibration->Anomalie = (vibration->Anomalie & ~LOOSENESS);

	/* Estimate speed using frequency near nominal RPM */
	if(vibration->MotorState == 1)
		vibration->Speed = (uint16_t)mainFreq*60;
	else
		vibration->Speed = 0;

	/* Sensor State */
	vibration->SensorState = getStatus();

	/* Detect Anomaly based on historical trend */
	if(getAnBehavior())
		vibration->Anomalie = (vibration->Anomalie | BEHAVIOR);
	else
		vibration->Anomalie = (vibration->Anomalie & ~BEHAVIOR);


}

uint8_t getStatus(void)
{
	return mb_status;
}

void setStatus(uint8_t sts)
{
	mb_status = sts;
}

// ============================================

uint8_t	getAnBehavior(void)
{
	return anBehavior;
}

void	setAnBehavior(uint8_t value)
{
	anBehavior = value;
}

uint8_t getEnLearning( void ){
	return enLearning;
}

void    setEnLearning( uint8_t value ){
	enLearning = 0;
}

// ============================================

uint8_t FFT_LibRevision(void)
{
	return(FFT_LIB_REV);
}

uint8_t Exponent(uint16_t value)
{
	uint8_t result = 0;
	while (((value >> result) & 1) != 1) result++;
	return(result);
}

void Swap(float *x, float *y)
{
	float temp = *x;
	*x = *y;
	*y = temp;
}

// FIR FILTER ==========================================
#define MAX_SAMPLES 1024
#define FILTER_ORDER 12

static float coefficients[FILTER_ORDER + 1] = { -7.86184503e-04, -1.39428843e-03, -3.05627359e-03, -5.32981930e-03,
											-7.60602981e-03, -9.27366244e-03, 9.90901435e-01, -9.27366244e-03,
											-7.60602981e-03, -5.32981930e-03, -3.05627359e-03, -1.39428843e-03,
											-7.86184503e-04 };

static float history[FILTER_ORDER + 1] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

void high_pass_filter_execute(float *input, int length)
{
    // Verificar que la longitud no supere el tamaño máximo permitido
    if (length > MAX_SAMPLES) {
        length = MAX_SAMPLES;
    }

    for (int i = 0; i < length; i++) {
        // Desplazar los valores del historial
        for (int j = FILTER_ORDER; j > 0; j--) {
            history[j] = history[j - 1];
        }

        // Agregar el nuevo valor al historial
        history[0] = input[i];

        // Realizar la convolución entre los coeficientes y el historial
        float filtered_value = 0.0;
        for (int j = 0; j <= FILTER_ORDER; j++) {
            filtered_value += coefficients[j] * history[j];
        }

        // Almacenar el valor filtrado directamente en el arreglo de salida (input = output)
        input[i] = filtered_value;
    }
}

void mean(float *signal, uint16_t len, float *mean) {
    float sum = 0.0;

    for (uint16_t i = 0; i < len; i++) {
        sum = sum + signal[i];
    }

    *mean = sum / len;
}

void    deviation(float *signal, uint16_t len, float mean, float *deviation){
	float sum = 0.0;

	for (uint16_t i = 0; i < len; i++) {
		sum = sum + (signal[i] - mean)*(signal[i] - mean);
	}

	sum = sqrtf(sum/len);
	*deviation = sum;

}