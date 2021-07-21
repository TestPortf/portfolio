#include <stdio.h>
#include "C8051F120.h"

sbit LED = P1^6;
cnt = 0x10;

void Reset_Sources_Init()
{
    WDTCN     = 0xDE;
    WDTCN     = 0xAD;
}

void Timer_Init()
{
    SFRPAGE   = TIMER01_PAGE;
    TMOD      = 0x01;
}

void Port_IO_Init()
{
    SFRPAGE   = CONFIG_PAGE;
    P1MDOUT   = 0x40;
    XBR2      = 0x40;
}

void Interrupts_Init()
{
    IE        = 0x82;
}

// Функция инициализации для устройства,
// вызов Device () из основной программы
void Init_Device(void)
{
    Reset_Sources_Init();
    Port_IO_Init();
    Interrupts_Init();
		Timer_Init();
}

Timer_ISR(void) interrupt 1{
    TF0 = 0;    //чистка флага прерывания
}

void main(void){
		Init_Device();
		TR0 = 1; // Включаем таймер
		LED = 0; // Выключаем диод
		while (1){
				PCON |= 0x1; // Включаем режим ожидания
				cnt--;
				if (cnt == 0){
						LED=~LED;//Через несколько прерываний инвертируем состояние диода
						cnt = 0x10; // Заново устанавливаем задержку смены состояния диода
				}
		}
}

