#define solenoide 13

void setup() {
  Serial.begin(9600);
  pinMode(solenoide, OUTPUT);
}

void loop() {
  if (Serial.available() > 0) {
    for (int i = 10; i > 0; i--) {
      Serial.println(i);
      delay(500);
    }
    digitalWrite(solenoide, HIGH);
    Serial.println("ON");
    delay(5000);
    digitalWrite(solenoide, LOW);
    Serial.println("OFF");
  }
}
