#include "stm32f4xx.h"

#define CTRL_REG1 														 ((uint8_t) 0x20)
#define OUTX														 			 ((uint8_t) 0x29)
#define OUTY																	 ((uint8_t) 0x2B)

#define CTRL_REG1_XEN													 ((uint8_t) 0x01)
#define CTRL_REG1_YEN													 ((uint8_t) 0x02)
#define CTRL_REG1_FS													 ((uint8_t) 0x20)
#define CTRL_REG1_PD													 ((uint8_t) 0x40)
#define CTRL_REG1_DR													 ((uint8_t) 0x80)

#define TRANSMIT_RW														 ((uint8_t) 0x80)

void SPI_Tx(uint8_t address, uint8_t data);
uint8_t SPI_Rx(uint8_t address);
