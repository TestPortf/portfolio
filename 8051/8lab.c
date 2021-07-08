#include <stdio.h>
#include "C8051F120.h"

sbit LED=P1^6; //указание на светодиод, подключенный к выводу P1.6
sbit SW=P3^7; //указание на кнопку, подключенную к выводу P3.7
int i; 
i = 0;
int SW_SHD; 
int array [8] = {0xFF,0x80,0x40,0x20,0x10,0x8,0x4,0x2};

void Reset_Sources_Init() //отключение сторожевого таймера
{
    WDTCN     = 0xDE;
    WDTCN     = 0xAD;
}

void PCA_Init() работа с ПМС
{
    SFRPAGE   = PCA0_PAGE;
    PCA0CN    = 0x40; //включение ПМС
    PCA0CPH4  = array[i]; //устанавливаем скважность ШИМ
    PCA0CPM4  = 0x42; //установки ПМС в режим 8-разрядного ШИМ
		
}

void Port_IO_Init()
{
    SFRPAGE   = CONFIG_PAGE;
    P1MDOUT   = 0x40; //задание двухтактного выхода p1.6 (Push-pull output)
    XBR0      = 0x2F;
    XBR2      = 0x44;
}

// Функция инициализации для устройства,
// вызов Device () из основной программы
void Init_Device(void) //Набор вызовов функций для инициализации устройства
{
    Reset_Sources_Init();
    PCA_Init();
    Port_IO_Init();
}

void main(void) {
	Init_Device();
	LED = 0; //выключение светодиода
	SW_SHD = SW; //запоминаем состояние кнопки

	while(1){
	  if(SW != SW_SHD) { //проверяем изменилось ли состояние кнопки
			SW_SHD = SW; //если да, то запоминаем его
			if (SW_SHD==0){ //если кнопка нажата, то
			  i = (i+1)%8; //итерируемся по массиву индексов
			  SFRPAGE = PCA0_PAGE; //эта команда нужна, чтобы можно было записать что либо в  PCA0CPH4
			  PCA0CPH4 = array[i]; //устанавливаем скважность ШИМ
			}	
		}
	}	
}
