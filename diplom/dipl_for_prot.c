#ifndef __AVR_ATmega328P__
#define __AVR_ATmega328P__
#endif

#define A_B  7
#define BUF  6
#define GA  5
#define SHDN  4

#define ACCUMULATOR_COMPARISON_SHIFT 24
#define COMPARISON_MASK (1 << ACCUMULATOR_COMPARISON_SHIFT) - 1
#define SIGNAL_DC 2048


#define SIN_LENGTH 256
#define PORTD3_REG_POSITION 8 
#define PORTD4_REG_POSITION 16
#define SIN_WAVE 2
#define TRIANGLE_WAVE 3
#define SAWTOOTH_WAVE 0
#define MEANDR_WAVE 1
#define KEY_PRESSED 1
#define NUMBER_OF_SAMPLES 15

#include <avr/io.h>
#include <util/delay.h>
#include <avr/interrupt.h>
int sin_wave[256] =  {0,    50,   100,   151,   201,   251,   301,   351,   400,   450,   499,   548,   596,   644,   692,   739, \
                        786,   832,   878,   923,   968,  1012,  1056,  1099,  1141,  1182,  1223,  1263,  1303,  1341,  1379,  1416, \
                       1452,  1487,  1521,  1554,  1587,  1618,  1649,  1678,  1706,  1734,  1760,  1785,  1809,  1832,  1854,  1875, \
                       1895,  1913,  1931,  1947,  1962,  1976,  1988,  2000,  2010,  2019,  2026,  2033,  2038,  2042,  2045,  2047, \
                       2047,  2046,  2044,  2040,  2036,  2030,  2023,  2014,  2005,  1994,  1982,  1969,  1954,  1939,  1922,  1904, \
                       1885,  1865,  1843,  1821,  1797,  1773,  1747,  1720,  1692,  1663,  1633,  1603,  1571,  1538,  1504,  1469, \
                       1434,  1397,  1360,  1322,  1283,  1243,  1203,  1162,  1120,  1077,  1034,   990,   946,   901,   855,   809, \
                        763,   716,   668,   620,   572,   523,   474,   425,   376,   326,   276,   226,   176,   126,    75,    25, \
                        -25,   -75,  -126,  -176,  -226,  -276,  -326,  -376,  -425,  -474,  -523,  -572,  -620,  -668,  -716,  -763, \
                       -809,  -855,  -901,  -946,  -990, -1034, -1077, -1120, -1162, -1203, -1243, -1283, -1322, -1360, -1397, -1434, \
                      -1469, -1504, -1538, -1571, -1603, -1633, -1663, -1692, -1720, -1747, -1773, -1797, -1821, -1843, -1865, -1885, \
                      -1904, -1922, -1939, -1954, -1969, -1982, -1994, -2005, -2014, -2023, -2030, -2036, -2040, -2044, -2046, -2047, \
                      -2047, -2045, -2042, -2038, -2033, -2026, -2019, -2010, -2000, -1988, -1976, -1962, -1947, -1931, -1913, -1895, \
                      -1875, -1854, -1832, -1809, -1785, -1760, -1734, -1706, -1678, -1649, -1618, -1587, -1554, -1521, -1487, -1452, \
                      -1416, -1379, -1341, -1303, -1263, -1223, -1182, -1141, -1099, -1056, -1012,  -968,  -923,  -878,  -832,  -786, \
                       -739,  -692,  -644,  -596,  -548,  -499,  -450,  -400,  -351,  -301,  -251,  -201,  -151,  -100,   -50, 0};

int meandr_wave[256] =  {0,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047, \
                          2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047, \
                          2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047, \
                          2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047, \
                          2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047, \
                          2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047, \
                          2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047, \
                          2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047,  2047, \
                         -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, \
                         -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, \
                         -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, \
                         -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, \
                         -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, \
                         -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, \
                         -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, \
                         -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047, -2047};

int sawtooth_wave[256] =  {0,    16,    32,    48,    64,    80,    96,   112,   128,   144,   160,   176,   192,   208,   224,   240, \
                             256,   273,   289,   305,   321,   337,   353,   369,   385,   401,   417,   433,   449,   465,   481,   497, \
                             513,   529,   546,   562,   578,   594,   610,   626,   642,   658,   674,   690,   706,   722,   738,   754, \
                             770,   786,   802,   818,   835,   851,   867,   883,   899,   915,   931,   947,   963,   979,   995,  1011, \
                            1027,  1043,  1059,  1075,  1092,  1108,  1124,  1140,  1156,  1172,  1188,  1204,  1220,  1236,  1252,  1268, \
                            1284,  1300,  1316,  1332,  1348,  1364,  1381,  1397,  1413,  1429,  1445,  1461,  1477,  1493,  1509,  1525, \
                            1541,  1557,  1573,  1589,  1605,  1621,  1637,  1654,  1670,  1686,  1702,  1718,  1734,  1750,  1766,  1782, \
                            1798,  1814,  1830,  1846,  1862,  1878,  1894,  1910,  1927,  1943,  1959,  1975,  1991,  2007,  2023,  2039, \
                           -2039, -2023, -2007, -1991, -1975, -1959, -1943, -1927, -1911, -1894, -1878, -1862, -1846, -1830, -1814, -1798, \
                           -1782, -1766, -1750, -1734, -1718, -1702, -1686, -1670, -1654, -1638, -1621, -1605, -1589, -1573, -1557, -1541, \
                           -1525, -1509, -1493, -1477, -1461, -1445, -1429, -1413, -1397, -1381, -1365, -1348, -1332, -1316, -1300, -1284, \
                           -1268, -1252, -1236, -1220, -1204, -1188, -1172, -1156, -1140, -1124, -1108, -1092, -1075, -1059, -1043, -1027, \
                           -1011,  -995,  -979,  -963,  -947,  -931,  -915,  -899,  -883,  -867,  -851,  -835,  -819,  -802,  -786,  -770, \
                            -754,  -738,  -722,  -706,  -690,  -674,  -658,  -642,  -626,  -610,  -594,  -578,  -562,  -545,  -529,  -513, \
                            -497,  -481,  -465,  -449,  -433,  -417,  -401,  -385,  -369,  -353,  -337,  -321,  -305,  -289,  -273,  -256, \
                            -240,  -224,  -208,  -192,  -176,  -160,  -144,  -128,  -112,   -96,   -80,   -64,   -48,   -32,   -16, 0};

int triangle_wave[256] =  {0,    32,    64,    96,   128,   160,   192,   224,   256,   289,   321,   353,   385,   417,   449,   481, \
                             513,   546,   578,   610,   642,   674,   706,   738,   770,   802,   835,   867,   899,   931,   963,   995, \
                            1027,  1059,  1091,  1124,  1156,  1188,  1220,  1252,  1284,  1316,  1348,  1381,  1413,  1445,  1477,  1509, \
                            1541,  1573,  1605,  1638,  1670,  1702,  1734,  1766,  1798,  1830,  1862,  1894,  1927,  1959,  1991,  2023, \
                            2039,  2007,  1975,  1943,  1911,  1878,  1846,  1814,  1782,  1750,  1718,  1686,  1654,  1621,  1589,  1557, \
                            1525,  1493,  1461,  1429,  1397,  1365,  1332,  1300,  1268,  1236,  1204,  1172,  1140,  1108,  1075,  1043, \
                            1011,   979,   947,   915,   883,   851,   819,   786,   754,   722,   690,   658,   626,   594,   562,   529, \
                             497,   465,   433,   401,   369,   337,   305,   273,   240,   208,   176,   144,   112,    80,    48,    16, \
                             -16,   -48,   -80,  -112,  -144,  -176,  -208,  -240,  -272,  -305,  -337,  -369,  -401,  -433,  -465,  -497, \
                            -529,  -562,  -594,  -626,  -658,  -690,  -722,  -754,  -786,  -818,  -851,  -883,  -915,  -947,  -979, -1011, \
                           -1043, -1075, -1108, -1140, -1172, -1204, -1236, -1268, -1300, -1332, -1364, -1397, -1429, -1461, -1493, -1525, \
                           -1557, -1589, -1621, -1654, -1686, -1718, -1750, -1782, -1814, -1846, -1878, -1910, -1943, -1975, -2007, -2039, \
                           -2023, -1991, -1959, -1927, -1894, -1862, -1830, -1798, -1766, -1734, -1702, -1670, -1638, -1605, -1573, -1541, \
                           -1509, -1477, -1445, -1413, -1381, -1348, -1316, -1284, -1252, -1220, -1188, -1156, -1124, -1091, -1059, -1027, \
                            -995,  -963,  -931,  -899,  -867,  -835,  -802,  -770,  -738,  -706,  -674,  -642,  -610,  -578,  -546,  -513, \
                            -481,  -449,  -417,  -385,  -353,  -321,  -289,  -256,  -224,  -192,  -160,  -128,   -96,   -64,   -32, 0};


unsigned long int phase_incriment_oct[4][13] = {{25538070, 27056536, 28665392, 30369908, 32175942, 34089155, 36116184, 38263863, 40539024, 42949672, 45503616, 48209446, 51076141},                                                    
                                                {51076141, 54113269, 57329028, 60740013, 64351885, 68178310, 72232369, 76527532, 81078245, 85899345, 91007233, 96418696, 102152087},                                                    
                                                {102152087, 108226342, 114661960, 121480026, 128703575, 136356621, 144464934, 153055259, 162156295, 171798691, 182014466, 192837589, 204303785},                                                    
                                                {204303785, 216452685, 229323921, 242960443, 257407151, 272712852, 288930259, 306110128, 324312980, 343597383, 364029714, 385674397, 408607570}};

uint16_t attack_rates[16] = {1, 2, 3, 4, 5, 8, 12, 20, 32, 37, 43, 51, 64, 85, 128, 255};
uint8_t decay_rates[16] = {2, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60};
uint8_t release_rates[16] = {1, 2, 3, 4, 5, 8, 12, 20, 32, 37, 43, 51, 64, 85, 128, 255};


uint8_t pushed_key;
volatile uint8_t current_note_number = 13;
volatile uint8_t last_note_number =0;
volatile uint8_t KEY_PRESSED_FLAG;
volatile int dac_sent_flag;
volatile uint16_t dac_data;

void spi_dac_sent(int value_to_sent){
    // Сообщаем о начале передачи по spi
    PORTB &= ~(1<<PORTB2);
    // Посылаем старшие 4 бита данных и 4 управляющих ЦАПом бита 
    SPDR = (value_to_sent >> 8) | (1 << GA) | (1 << SHDN);
    // Ждём окончания передачи
    while(!(SPSR & (1<<SPIF)));
    // Передаём младшие 8 бит данных
    SPDR = value_to_sent & 0x00FF;
    // Ждём окончания передачи
    while(!(SPSR & (1<<SPIF)));
    // Сообщаем об окончании передачи по spi
    PORTB |= (1<<PORTB2);
}

void spi_init(void){
    // Включаем ss, mosi и sck для spi
    DDRB |= ((1 << PORTB2) | (1 << PORTB3) | (1 << PORTB5)); 
    PORTB &= ~((1 << PORTB2) | (1 << PORTB3) | (1 << PORTB5));
    // Включим шину spi и обявим устройство ведущим
    SPCR = ((1 << SPE) | (1 << MSTR));
    SPSR |= SPI2X;
}

void adc_init(void){
    // Включение ацп и предделителя частоты его работы 128
    ADCSRA |= (1 << ADEN) | (1 << ADPS2) | (1 << ADPS1) | (1 << ADPS0);
    // Конфигурирование формата вывода значения ацп таким образом, чтобы
    // в регистре ADCH находилось 8 старших бит значений преобразования ацп
    ADMUX |= (1 << ADLAR);
}

void key_init(void){
    // Устанавливаем подтягивающие разисторы на выводы
    // подключенные к ёмкостным клавишам 
    PORTC |= (1 << PORT0) | (1 << PORT1);
    PORTB |= (1 << PORT0) | (1 << PORT1) | (1 << PORT4);
    PORTD = 0xFF;
}


uint8_t get_key_rise_time_PORTC(uint8_t key_to_check, uint8_t key_pressed_time){
    uint8_t time = 0;
    
    // Выводим логический ноль на выбранном порту ввода-вывода
    PORTC &= ~(key_to_check); 
    DDRC  |= key_to_check;
    asm("nop");
    asm("nop");
    // состоянием порта
    cli();
    // Конфигурируем порт для работы в режиме входа
    DDRC &= ~(key_to_check); 
    // Включаем подтяжку к питанию
    PORTC |= key_to_check;   
    // Проверяем, устанавливается ли единица за отведённое время
    for (time = 0; time < key_pressed_time; time++){
        if (PINC & key_to_check){
            // Включаем прерывания
            sei();
            // Если установилась, то кнопка не нажата
            return 0;
        }
    }
    // Если за отведённое время не установилась, то
    // кнопка нажата
    return 1;
}

uint8_t get_key_rise_time_PORTB(uint8_t key_to_check, uint8_t key_pressed_time){
    uint8_t time = 0;
    
    // Выводим логический ноль на выбранном порту ввода-вывода
    PORTB &= ~(key_to_check); 
    DDRB  |= key_to_check;
    asm("nop");
    asm("nop");
    // состоянием порта
    cli();
    // Конфигурируем порт для работы в режиме входа
    DDRB &= ~(key_to_check); 
    // Включаем подтяжку к питанию
    PORTB |= key_to_check;   
    // Проверяем, устанавливается ли единица за отведённое время
    for (time = 0; time < key_pressed_time; time++){
        if (PINB & key_to_check){
            // Включаем прерывания
            sei();
            // Если установилась, то кнопка не нажата
            return 0; 
        }
    }
    // Если за отведённое время не установилась, то
    // кнопка нажата
    return 1;
}

uint8_t get_key_rise_time_PORTD(uint8_t key_to_check, uint8_t key_pressed_time){
    uint8_t time = 0;
    
    // Выводим логический ноль на выбранном порту ввода-вывода
    PORTD &= ~(key_to_check); 
    DDRD  |= key_to_check;
    asm("nop");
    asm("nop");
    // состоянием порта
    cli();
    // Конфигурируем порт для работы в режиме входа
    DDRD &= ~(key_to_check); 
    // Включаем подтяжку к питанию
    PORTD |= key_to_check;   
    // Проверяем, устанавливается ли единица за отведённое время
    for (time = 0; time < key_pressed_time; time++){
        if (PIND & key_to_check){
            // Включаем прерывания
            sei();
            // Если установилась, то кнопка не нажата
            return 0;
        }
    }
    // Если за отведённое время не установилась, то
    // кнопка нажата
    return 1;
}

void timer_init(void){
    // Устанавливаем режим СТС (сброс по совпадению)
    TCCR2A |= (1<<WGM21); 
    // Устанавливаем бит разрешения прерывания 2ого счетчика по совпадению с OCR0A
    TIMSK2 |= (1<<OCIE2A); 
    // Записываем количество отсчётов таймера до генерации отсчёта сигнала
    // Оно выбрано так, чтобы имея ввиду делитель 8 частота дискретизации
    // была равна ~22000 Гц
    OCR2A = 0x72; 
    // Разрешаем прерывания
    sei();
    // Включаем таймер с предделителем тактового сигнала равного 8
    TCCR2B |= 1 << CS21;
}

ISR (TIMER2_COMPA_vect){
    // Принимается отсчёт сигнала
    uint16_t dac_data_isr = dac_data;
    // Посылается на ЦАП
    spi_dac_sent(dac_data_isr);
    // Устанавливается флаг синхронизации главной функции
    // с прерыванием
    dac_sent_flag = 1;
}

void key_pressed(int key_number){
    last_note_number = current_note_number;
    current_note_number = key_number;
    KEY_PRESSED_FLAG = KEY_PRESSED;
}

unsigned int interpolate(unsigned long int current_phase_index, int *table_pointer){
    uint8_t nearest_table_index;
    int y_0, y_1;
    unsigned long int phase_accum_and_table_dif;

    nearest_table_index = current_phase_index >> ACCUMULATOR_COMPARISON_SHIFT;
    y_0 = table_pointer[nearest_table_index];
    if (nearest_table_index < 255){
        y_1 = table_pointer[nearest_table_index + 1];
    }
    else y_1 = table_pointer[0];

    phase_accum_and_table_dif = (COMPARISON_MASK & current_phase_index); 
    return (y_0 * ((1 << ACCUMULATOR_COMPARISON_SHIFT) - phase_accum_and_table_dif) + y_1 * phase_accum_and_table_dif) >> ACCUMULATOR_COMPARISON_SHIFT;
}


int main(void){
    uint8_t adc_select = 2;

    uint8_t selected_octave;
    uint8_t signal_waveshape;
    uint8_t attack_duration;
    uint8_t decay_duration;
    uint8_t sustain_level;
    uint8_t release_duration;

    float volume;
    uint8_t ATTACK_END_FLAG;
    uint8_t DECAY_END_FLAG;

    uint8_t current_key;
    uint8_t touch_recognition_time = 4;
    KEY_PRESSED_FLAG = !KEY_PRESSED;

    unsigned long int phase_index = 0;
    int *wave_pointer;

    uint8_t dac_sent_flag_local;

    adc_init();
    spi_init();
    key_init();
    timer_init();

    while(1){
        KEY_PRESSED_FLAG = !KEY_PRESSED;

        // Опрос одной клавиши раз в цикл
        switch (current_key++){
            case 0:
                if (get_key_rise_time_PORTD((1 << PIN0), touch_recognition_time)) key_pressed(current_key - 1);
                break;

            case 1:
                if (get_key_rise_time_PORTD((1 << PIN1), touch_recognition_time)) key_pressed(current_key - 1);
                break;

            case 2:
                if (get_key_rise_time_PORTD((1 << PIN2), touch_recognition_time)) key_pressed(current_key - 1);
                break;

            case 3:
                if (get_key_rise_time_PORTD((1 << PIN3), touch_recognition_time)) key_pressed(current_key - 1);
                break;

            case 4:
                if (get_key_rise_time_PORTD((1 << PIN4), touch_recognition_time)) key_pressed(current_key - 1);
                break;

            case 5:
                if (get_key_rise_time_PORTD((1 << PIN5), touch_recognition_time)) key_pressed(current_key - 1);
                break;

            case 6:
                if (get_key_rise_time_PORTD((1 << PIN6), touch_recognition_time)) key_pressed(current_key - 1);
                break;

            case 7:
                if (get_key_rise_time_PORTD((1 << PIN7), touch_recognition_time)) key_pressed(current_key - 1);
                break;

            case 8:
                if (get_key_rise_time_PORTB((1 << PIN0), touch_recognition_time)) key_pressed(current_key - 1);
                break;

            case 9:
                if (get_key_rise_time_PORTB((1 << PIN1), touch_recognition_time)) key_pressed(current_key - 1);
                break;

            case 10:
                if (get_key_rise_time_PORTB((1 << PIN4), touch_recognition_time)) key_pressed(current_key - 1);
                break;

            case 11:
                if (get_key_rise_time_PORTC((1 << PIN0), touch_recognition_time)) key_pressed(current_key - 1);
                break;

            case 12:
                if (get_key_rise_time_PORTC((1 << PIN1), touch_recognition_time)) key_pressed(current_key - 1);
                current_key = 0;
                break;
        }

        // Опрос одного вывода ацп раз в цикл
        if (adc_select == 2) selected_octave = ADCH >> 6;
        if (adc_select == 3){
            signal_waveshape = ADCH >> 6;
            if (signal_waveshape == SIN_WAVE) wave_pointer = sin_wave;
            if (signal_waveshape == MEANDR_WAVE) wave_pointer = meandr_wave;
            if (signal_waveshape == SAWTOOTH_WAVE) wave_pointer = sawtooth_wave;
            if (signal_waveshape == TRIANGLE_WAVE) wave_pointer = triangle_wave;
        }
        if (adc_select == 4) attack_duration = attack_rates[ADCH >> 4];
        if (adc_select == 5) decay_duration = decay_rates[ADCH >> 4];
        if (adc_select == 6) sustain_level = ADCH;
        if (adc_select == 7) release_duration = release_rates[ADCH >> 4];
        adc_select++;
        if (adc_select > 7) adc_select = 2;
        ADMUX = (1 << ADLAR) | adc_select;
        ADCSRA |= ADLAR;

        // Если нажата другая клавиша, то все значения и флаги должны быть сброшены
        if (current_note_number != last_note_number){
            volume = 0;    
            ATTACK_END_FLAG = 0;
            DECAY_END_FLAG = 0;
        }

        // Если клавиша нажата, то происходит увелечение фазового аккумулятора на
        // фазовый инкремент соотвествующий выбранной октаве и ноте
        if (KEY_PRESSED_FLAG == KEY_PRESSED & current_note_number == last_note_number) {
            phase_index += phase_incriment_oct[selected_octave][current_note_number];

            // Если фазы атаки ещё не было, то сигнал постепенно нарастает
            // фаза атаки заканчивается при достижении максимальной громкости
            if (!ATTACK_END_FLAG){
                volume += 1/attack_duration;
                if (volume >= 1) ATTACK_END_FLAG = 1;
            }
            // Если фазы спада ещё не было, то убавляй сигнал
            // до уровня установленного значением задержки
            else if (!DECAY_END_FLAG){
                volume -= 1/decay_duration;
                if (volume < sustain_level/255) DECAY_END_FLAG = 1;
            }

        }

        // Если нажатая клавиша была отпущена, то постепенно убавляй громкость сигнала
        if (KEY_PRESSED_FLAG == !KEY_PRESSED & volume >= 0){
            phase_index += phase_incriment_oct[selected_octave][current_note_number];
            while(volume > 0){                    
                volume -= 1/release_duration;
            }
        }

        // интерполяция и амплитудная модуляция сигнала
        dac_data = volume * interpolate(phase_index, wave_pointer) + SIGNAL_DC;

        // Ожидание прерывания от таймера
        while (!dac_sent_flag_local){
            dac_sent_flag_local = dac_sent_flag;
        }

        // Обнуляем флаг, чтобы снова ждать прерывания
        // в следующем цикле
        dac_sent_flag = 0;
    }
}