#include "main.h"
UART_HandleTypeDef huart2;
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_USART2_UART_Init(void);
__STATIC_INLINE void DelayMicro(__IO uint32_t micros)
{
	micros *=(SystemCoreClock / 1000000) / 5;
	while (micros--);
}
int main(void)
{
	uint8_t str[1];
  HAL_Init();
  SystemClock_Config();
  MX_USART2_UART_Init();
  while (1)
  {
		uint8_t state = HAL_UART_GetState(&huart2);
		// Если передача не идёт
		if( (state != HAL_UART_STATE_BUSY_RX) && (state != HAL_UART_STATE_BUSY_TX_RX) ) {

    while( HAL_UART_Transmit_IT(&huart2, str, 1) == HAL_BUSY ); //передаём данные
    HAL_UART_Receive_IT (&huart2, str, 1); //принимаем данные
		}
  }
}

