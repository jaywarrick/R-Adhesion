
// Number of milliseconds between samples
const int SAMPLE_INTERVAL = 10;
int sensorValue = 0;        // value read from the pot
int outputValue = 0;        // value output to the PWM (analog out)
int duration = 1000;
int start = 0;
double current = 0.0;
double voltage = 1;
double DACMAX = 3.3;
double t = 0;
double A = 0.8;            // Amplitude
double offset = 0; //double A plus a small amount
double tf = 300;         //# of sec
double pi = 3.14159;
double fi = .01;           //frequency initial
double ff = 1;         // frequency final
double F = (ff/fi);
double N = ((log(ff/fi))/(log(2)));
double R = (N/tf);
double P = 0.0;
double X = 0.0;


void setup() {
//  pinMode(DAC1, OUTPUT);
  double F = (ff/fi);
  double N = ((log(ff/fi))/(log(2)));
  double R = (N/tf);
  analogWriteResolution(12);
  Serial.begin(9600);
  start = millis();
  pinMode(3,INPUT);
  if (digitalRead(3)){
    fi = 1;
    ff = .01;
    F = (ff/fi);
    N = ((log(ff/fi))/(log(2)));
    R = (N/tf);
  }
}

void loop() {
if (t < 1000*tf)
{
  t = (millis() - start);
  P = double(pow(2, (R*((t)/1000))));
  X = (sin(2*pi*((fi*(-1+(P)))/(R*log(2)))));
  voltage= ((-A)*(asin(X)))+offset;
  outputValue = int((4095.0/DACMAX)*voltage+2048.0);
  Serial.println(outputValue);
  analogWrite(DAC1, outputValue);
  //delay(2);
}
else
{
 t = 1000*tf;
 P = double(pow(2, (R*((t)/1000))));
  X = (sin(2*pi*((fi*(-1+(P)))/(R*log(2)))));
  voltage= ((-A)*(asin(X)))+offset;
  outputValue = int((4095.0/DACMAX)*voltage+2048.0);
  Serial.println(outputValue);
  analogWrite(DAC1, outputValue);
}
}

