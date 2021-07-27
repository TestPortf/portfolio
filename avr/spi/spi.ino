#include <SPI.h>

int hexDigitValue[] = {0x03, /* 0 */
		       0x9F, /* 1 */
		       0x25, /* 2 */
		       0x0D, /* 3 */
		       0x99, /* 4 */
		       0x49, /* 5 */
		       0x41, /* 6 */
		       0x1F, /* 7 */
		       0x01, /* 8 */
		       0x09};/* 9 */

const int slaveSelectPin = 10;
int button_number;

void button_init(){
    pinMode(7, INPUT_PULLUP);
    pinMode(6, INPUT_PULLUP);
    pinMode(5, INPUT_PULLUP);
    pinMode(4, INPUT_PULLUP);
    pinMode(3, INPUT_PULLUP);
    pinMode(2, INPUT_PULLUP);
    pinMode(1, INPUT_PULLUP);
    pinMode(0, INPUT_PULLUP);
}

void setup() {
    // Инициализируем кнопки
    button_init();
    // Инициализируем SPI и ножку отвечающую за переход в режим хранения
    pinMode(slaveSelectPin, OUTPUT);
    SPI.begin();
}

void loop() {
    // Считываем значение, подаваемое на кнопки в виде 8-разрядного числа
    button_number = PIND;
    // Передаём данные
    digitalWrite(slaveSelectPin, HIGH);
    SPI.transfer(hexDigitValue[button_number % 10]);
    SPI.transfer(hexDigitValue[(button_number / 10) % 10]);
    SPI.transfer(hexDigitValue[button_number / 100]);
    digitalWrite(slaveSelectPin, LOW);
}
