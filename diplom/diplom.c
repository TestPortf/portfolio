#ifndef __AVR_ATmega328P__
#define __AVR_ATmega328P__
#endif

#define A_B  7
#define BUF  6
#define GA  5
#define SHDN  4

#define SIN_LENGTH 256
#define PORTD3_REG_POSITION 8 
#define PORTD4_REG_POSITION 16
#define SAWTOOTH_WAVE 0
#define SQUARE_WAVE 1
#define KEY_PRESSED 1
#define NUMBER_OF_SAMPLES 15

#include <avr/io.h>
#include <util/delay.h>
#include <avr/interrupt.h>

int sin_wave[SIN_LENGTH] = {2048, 2098, 2148, 2199, 2249, 2299, 2349, 2399, 2448, 2498, 2547, 2596, 2644, 2692, 2740, 2787,\
                            2834, 2880, 2926, 2971, 3016, 3060, 3104, 3147, 3189, 3230, 3271, 3311, 3351, 3389, 3427, 3464,\
                            3500, 3535, 3569, 3602, 3635, 3666, 3697, 3726, 3754, 3782, 3808, 3833, 3857, 3880, 3902, 3923,\
                            3943, 3961, 3979, 3995, 4010, 4024, 4036, 4048, 4058, 4067, 4074, 4081, 4086, 4090, 4093, 4095,\
                            4095, 4094, 4092, 4088, 4084, 4078, 4071, 4062, 4053, 4042, 4030, 4017, 4002, 3987, 3970, 3952,\
                            3933, 3913, 3891, 3869, 3845, 3821, 3795, 3768, 3740, 3711, 3681, 3651, 3619, 3586, 3552, 3517,\
                            3482, 3445, 3408, 3370, 3331, 3291, 3251, 3210, 3168, 3125, 3082, 3038, 2994, 2949, 2903, 2857,\
                            2811, 2764, 2716, 2668, 2620, 2571, 2522, 2473, 2424, 2374, 2324, 2274, 2224, 2174, 2123, 2073,\
                            2022, 1972, 1921, 1871, 1821, 1771, 1721, 1671, 1622, 1573, 1524, 1475, 1427, 1379, 1331, 1284,\
                            1238, 1192, 1146, 1101, 1057, 1013,  970,  927,  885,  844,  804,  764,  725,  687,  650,  613,\
                             578,  543,  509,  476,  444,  414,  384,  355,  327,  300,  274,  250,  226,  204,  182,  162,\
                             143,  125,  108,   93,   78,   65,   53,   42,   33,   24,   17,   11,    7,    3,    1,    0,\
                               0,    2,    5,    9,   14,   21,   28,   37,   47,   59,   71,   85,  100,  116,  134,  152,\
                             172,  193,  215,  238,  262,  287,  313,  341,  369,  398,  429,  460,  493,  526,  560,  595,\
                             631,  668,  706,  744,  784,  824,  865,  906,  948,  991, 1035, 1079, 1124, 1169, 1215, 1261,\
                            1308, 1355, 1403, 1451, 1499, 1548, 1597, 1647, 1696, 1746, 1796, 1846, 1896, 1947, 1997, 2047};

int square_wave[256] = {   0, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095,\
                        4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095,\
                        4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095,\
                        4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095,\
                        4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095,\
                        4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095,\
                        4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095,\
                        4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095,\
                           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,\
                           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,\
                           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,\
                           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,\
                           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,\
                           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,\
                           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,\
                           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0};

//int sawtooth_oct_length[12] = {470, 444, 419, 396, 374, 353, 333, 314, 296, 279, 264, 249};
int signal_oct_int_step[12] = {8, 9, 9, 10, 11, 11, 12, 13, 13, 14, 15, 16};
int signal_oct_float_step[12] = {7, 2, 8, 3, 0, 6, 3, 0, 8, 7, 5, 4};
int sawtooth_oct_length[12] = {235, 222, 209, 198, 187, 176, 166, 157, 148, 139, 132, 124};
//int signal_oct_int_step[12] = {4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8};
//int signal_oct_float_step[12] = {3, 1, 4, 1, 0, 3, 1, 0, 4, 3, 2, 2};

uint8_t pushed_key;
// Костыль
uint8_t current_note_number = 8;
int waveform_counter;
int waveform_step;
uint8_t selected_waveform = SQUARE_WAVE;
int waveform_period_counter;
int waveform_int_signal_level;
int waveform_float_signal_level;
int tmp_sig;
uint8_t KEY_PRESSED_FLAG;
volatile static uint8_t dac_sent = 0;
volatile static uint8_t attack_counter = 0;

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

void key_init(void){
    // Загрузчик мешает корректрой работе устройства
    // на arduino nano этот регистр нужно отчищать 
    UCSR0B = 0; 
    //// Конфигурируем 0, 1 и 2-й выводы мк как работающие на выход
    DDRD |= (1 << PORT0) | (1 << PORT1) | (1 << PORT2);
    //// Устанавливаем потягивающие резисторы на выводах 3 и 4
    PORTD |= (1 << PORT3) | (1 << PORT4);


    // Проверка
    DDRD |= (1 << PORT6);
}

void release_key(uint8_t key_to_check){
    
    // Выводим логический ноль на выбранном порту ввода-вывода
    PORTD &= ~(key_to_check); 
    DDRD  |= key_to_check;
    asm("nop"); 
    asm("nop");
    // Конфигурируем порт для работы в режиме входа
    DDRD &= ~(key_to_check); 
    // Включаем подтяжку к питанию
    PORTD |= key_to_check;   
    asm("nop"); 
    asm("nop");

    PORTD &= ~(key_to_check); 
    DDRD  |= key_to_check;
}

uint8_t get_key_rise_time(uint8_t key_to_check, uint8_t key_pressed_time){
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
            //PORTD &= ~(key_to_check);   
            return 0;
        }
    }
    // Если за отведённое время не установилась, то
    // кнопка нажата
    //PORTD &= ~(key_to_check);   
    return 1;
}

void timer_init(void){
    // Устанавливаем режим СТС (сброс по совпадению)
    TCCR2A |= (1<<WGM21); 
    // Устанавливаем бит разрешения прерывания 2ого счетчика по совпадению с OCR0A
    TIMSK2 |= (1<<OCIE2A); 
    // Записываем количество отсчётов таймера до генерации отсчёта сигнала
    // Оно выбрано так, чтобы имея ввиду делитель 8 частота дискретизации
    // была равна ~44100 Гц
    OCR2A = 0x4F; 
    // Разрешаем прерывания
    sei();
    // Включаем таймер с предделителем тактового сигнала равного 8
    TCCR2B |= 1 << CS21;
}

void pwm_init(void){
    PORTD &= ~(1 << PORT6);
    DDRD |= (1 << PORT6);
    OCR0A = 150;
    TCCR0A |= (1 << COM0A1) | (1 << COM0A0) | (1 << WGM00) | (1 << WGM01);
    TCCR0B |= (1 << CS01); 

}

ISR (TIMER2_COMPA_vect){
    if (KEY_PRESSED_FLAG){
            // Сообщаем о начале передачи по spi
            PORTB &= ~(1<<PORTB2);
            // Посылаем старшие 4 бита данных и 4 управляющих ЦАПом бита 
            SPDR = (waveform_int_signal_level >> 8) | (1 << GA) | (1 << SHDN);
            tmp_sig = waveform_int_signal_level & 0x00FF;
            // Ждём окончания передачи
        if (selected_waveform == SAWTOOTH_WAVE){
            waveform_int_signal_level += signal_oct_int_step[current_note_number];
            waveform_float_signal_level += signal_oct_float_step[current_note_number];
            if (waveform_float_signal_level > 10){
                waveform_int_signal_level++;
                waveform_float_signal_level = 0;
            }
        } 
        else if (selected_waveform == SQUARE_WAVE){
            waveform_int_signal_level = square_wave[waveform_period_counter];
        }
            while(!(SPSR & (1<<SPIF)));
            // Передаём младшие 8 бит данных
            SPDR = tmp_sig;
            // Ждём окончания передачи
        if (selected_waveform == SAWTOOTH_WAVE){
            waveform_period_counter++;
            if (waveform_period_counter > sawtooth_oct_length[current_note_number]){
                waveform_int_signal_level = 0;
                waveform_period_counter = 0;
            }
        }
        else if (selected_waveform == SQUARE_WAVE){
            waveform_period_counter++;
            if (waveform_period_counter > 255){
                waveform_period_counter = 0;
            } 
        }
            while(!(SPSR & (1<<SPIF)));
            // Сообщаем об окончании передачи по spi
            PORTB |= (1<<PORTB2);
    }
    else {
        waveform_int_signal_level = 0;
        waveform_period_counter = 0;
    }
}

int main(void){
    uint8_t current_key;
    uint8_t last_key;
    uint8_t prelast_key;
    uint8_t preprelast_key;
    uint8_t touch_recognition_time = 3;
    int pressed_array[NUMBER_OF_SAMPLES];
    int stat_cycles_counter;
    int not_pressed_cycles;
    int pressed_cycles;
    float pressed_mean = 0;
    uint8_t key_for_release;
    uint8_t wait_time = 4;
    spi_init();
    key_init();
    timer_init();

    while(1){
        for (current_key = 0; current_key < 13; current_key++) {
            if (current_key < 8){
                PORTD = (PORTD & ~7) | current_key;
                if ((get_key_rise_time(PORTD3_REG_POSITION, touch_recognition_time))){
                    selected_waveform = SAWTOOTH_WAVE; 
                    pressed_array[stat_cycles_counter] = current_key;
                    break;
                }
            }
            if (current_key >= 8){
                PORTD = (PORTD & ~7) | (current_key & 7);
                if ((get_key_rise_time(PORTD4_REG_POSITION, touch_recognition_time))){
                    pressed_array[stat_cycles_counter] = current_key;
                    selected_waveform = SAWTOOTH_WAVE; 
                    break;
                }
            }
            pressed_array[stat_cycles_counter] = 0;
        }
        stat_cycles_counter++;
        // Стабилизируем нажатие клавишы            
        if (stat_cycles_counter > NUMBER_OF_SAMPLES){
            for (stat_cycles_counter = 0; stat_cycles_counter < NUMBER_OF_SAMPLES; stat_cycles_counter++){                    
                if (pressed_array[stat_cycles_counter] == 0){
                    not_pressed_cycles++;
                }
                else {
                    pressed_cycles++;
                    pressed_mean += pressed_array[stat_cycles_counter];
                }
            }
            if (pressed_cycles > not_pressed_cycles){                    
                pressed_mean = pressed_mean / NUMBER_OF_SAMPLES;
                for (current_key = 0; current_key < 13; current_key++){                        
                    if (pressed_mean > current_key){
                        current_note_number = current_key;
                        KEY_PRESSED_FLAG = KEY_PRESSED;
                        PORTD &= ~(1 << PORTD6);
                    }
                }
            }
            else {
                KEY_PRESSED_FLAG = !KEY_PRESSED;
                PORTD |= (1 << PORTD6);
            }
            stat_cycles_counter = 0;
            not_pressed_cycles = 0;
            pressed_cycles = 0;
            pressed_mean = 0;
        } 

        //for (current_key = 0; current_key < 13; current_key++) {
        //    if (current_key < 8){
        //        PORTD = (PORTD & ~7) | current_key;
        //        if ((get_key_rise_time(PORTD3_REG_POSITION, touch_recognition_time))){
        //            preprelast_key = prelast_key;
        //            prelast_key = last_key;
        //            last_key = pushed_key;
        //            pushed_key = current_key;
        //            selected_waveform = SAWTOOTH_WAVE; 
        //            KEY_PRESSED_FLAG = KEY_PRESSED;      
        //            break;
        //        }
        //    }
        //    if (current_key >= 8){
        //        PORTD = (PORTD & ~7) | (current_key & 7);
        //        if ((get_key_rise_time(PORTD4_REG_POSITION, touch_recognition_time))){
        //            preprelast_key = prelast_key;
        //            prelast_key = last_key;
        //            last_key = pushed_key;
        //            pushed_key = current_key;
        //            selected_waveform = SAWTOOTH_WAVE; 
        //            KEY_PRESSED_FLAG = KEY_PRESSED;      
        //            break;
        //        }
        //    }
        //    //not_pressed_cycles++;
        //    KEY_PRESSED_FLAG = 0;
        //}
        //// Стабилизируем нажатие клавишы
        //if (pushed_key == last_key & last_key == prelast_key & preprelast_key == prelast_key & KEY_PRESSED_FLAG){
        //    current_note_number = pushed_key;
        //    not_pressed_cycles = 0;
        //}
        //else{
        //    not_pressed_cycles++;
        //}

        //// Прогоняем проверку другой клавиши, чтобы сбрасывать ёмкость
        //if (current_note_number == 0){
        //    get_key_rise_time(PORTD3_REG_POSITION, touch_recognition_time);
        //}

        //if (not_pressed_cycles > 1280){
        //    current_note_number = 13;
        //    selected_waveform = 3;
        //    waveform_int_signal_level = 0;
        //    not_pressed_cycles = 0;
        //}

        ////while (!dac_sent){
        ////    asm("nop");
        ////}

        //dac_sent = 0;
    }
    return 0;
}