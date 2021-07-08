#include <at89c51xd2.h>

#define bit_12_8 P3_3
#define STS 		 P1_7
#define CS			 P3_4
#define A0			 P3_5
#define CE			 P3_6
#define R_C			 P3_7
#define sync 		 P3_2

#define extract_first_byte(value) (value & 0xF) << 4
#define extract_second_byte(value) value >> 4

unsigned int reverse(unsigned int value)
{
	unsigned int result = 0;
	if ( value & 0x200 ) result |= 0x001;
	if ( value & 0x100 ) result |= 0x002;
	if ( value & 0x080 ) result |= 0x004;
	if ( value & 0x040 ) result |= 0x008;
	if ( value & 0x020 ) result |= 0x010;
	if ( value & 0x010 ) result |= 0x020;
	if ( value & 0x008 ) result |= 0x040;
	if ( value & 0x004 ) result |= 0x080;
	if ( value & 0x002 ) result |= 0x100;
	if ( value & 0x001 ) result |= 0x200;
	return result;
}

void shift_transmit(unsigned int data_to_transmit){
	char first_byte;
	char second_byte;
	
	data_to_transmit = reverse(data_to_transmit);

	first_byte = extract_first_byte(data_to_transmit);
	second_byte = extract_second_byte(data_to_transmit);
	
	sync = 0;
	SBUF = first_byte;
	while(TI==0);
	TI = 0; TI = 0;
	
	SBUF = second_byte;
	while(TI==0);
	TI = 0; TI = 0;
	sync = 1;
}

void Delay(int nCount){
	while(nCount--);
}

unsigned int ad1674_read(void){
	R_C = 0;
	while(STS==1);
	R_C = 1; R_C = 1;
	return P0 | ((P2 & 0xF) << 8);
}

unsigned int filter_signal(double A[5], double B[5], unsigned int current_point, unsigned int past_point){
	return 
}

void main (void)
{
	unsigned int current_adc_data;
	unsigned int past_adc_data = 1;
	double A_coef[5] = {1, -3.8261594180878564230852134642191231250762939453125,
											5.64021677662502707306657612207345664501190185546875,
											-3.788568424693382841184075005003251135349273681640625,
											0.98044753188059352577710114928777329623699188232421875};
	double B_coef[5] = {0.0000482615212977725109556002835997645661336719058454036712646484375,
											0, -0.000096523042595545021911200567199529132267343811690807342529296875,
											0, 0.0000482615212977725109556002835997645661336719058454036712646484375};
	
	bit_12_8 = 1; CS = 0; A0 = 0; CE = 1;
	
	while(1){
		current_adc_data = ad1674_read() >> 2;
		shift_transmit(current_adc_data);
	}
}