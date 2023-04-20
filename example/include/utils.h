// #include <immintrin.h>

void get_4_digits(int input, int *output)
{
    if (input < 10)
    {
        output[0] = 0;
        output[1] = 0;
        output[2] = 0;
        output[3] = input;
    }
    else if (input < 100)
    {
        output[0] = 0;
        output[1] = 0;
        output[2] = input / 10;
        output[3] = input % 10;
    }
    else if (input < 1000)
    {
        output[0] = 0;
        output[1] = input / 100;
        output[2] = input / 10 - input / 100 * 10;
        output[3] = input % 10;
    }
    else if (input < 10000)
    {
        output[0] = input / 1000;
        output[1] = input / 100 - input / 1000 * 10;
        output[2] = input / 10 - input / 100 * 10;
        output[3] = input % 10;
    }
}