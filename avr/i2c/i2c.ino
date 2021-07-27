#include <Wire.h>
#include <SevenSegPCF8574.h>
// Инициализация дисплеев
SevenSegPCF8574 high = SevenSegPCF8574(0x21, 7, 6, 5, 4, 3, 2, 1, 0);
SevenSegPCF8574 mid = SevenSegPCF8574(0x22, 7, 6, 5, 4, 3, 2, 1, 0);
SevenSegPCF8574 low = SevenSegPCF8574(0x24, 7, 6, 5, 4, 3, 2, 1, 0);

int panel;

void setup() {
    // Включает передачу по i2c
    Wire.begin();
    // Добавляет знак к каждому дисплею (по умолчанию - 0)
    high.addDigit(8);
    mid.addDigit(9);
    low.addDigit(10);
    // Инициализация панели с кнопками
    pinMode(7, INPUT_PULLUP);
    pinMode(6, INPUT_PULLUP);
    pinMode(5, INPUT_PULLUP);
    pinMode(4, INPUT_PULLUP);
    pinMode(3, INPUT_PULLUP);
    pinMode(2, INPUT_PULLUP);
    pinMode(1, INPUT_PULLUP);
    pinMode(0, INPUT_PULLUP);
}
void loop() {
    // Считывание значения с панели с кнопками
    panel = PIND;
    // Установка значений для каждого дислея (перевод в строку используется
    // чтобы выводить 0. Если выводить числами 0 почему-то не отображается)
    high.setString(String(panel/100));
    mid.setString(String((panel/10)%10));
    low.setString(String(panel%10));
    // Выводим значения дисплеев
    high.show();
    mid.show();
    low.show();
}
