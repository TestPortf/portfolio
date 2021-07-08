#include "stm32f4xx.h"

#define PLL_M				8
#define PLL_N				336
#define PLL_P				2
#define PLL_Q				7
#define AHB_PRE			RCC_CFGR_HPRE_DIV1
#define APB1_PRE		RCC_CFGR_PPRE1_DIV4	
#define APB2_PRE		RCC_CFGR_PPRE2_DIV2		

void LEDs_ini(void);
void TIM4_ini(void);
void SPI1_ini(void);
