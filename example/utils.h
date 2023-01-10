// #include <mpi.h>
#include <immintrin.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

typedef float data_type;


bool verify_arrays(data_type *array1, data_type *array2, int n)
{
  data_type diff = 0.f;
  int i;

  for (i = 0; i < n; i++)
  {
    diff = fabs(array1[i] - array2[i]);
    if (diff > 1e-4)
    {
      printf("error. %f,%f,%d\n", array1[i], array2[i], i);
      return false;
    }
  }

  return true;
}

// Creates an array of random numbers. Each number has a value from 0 - 1
data_type *create_rand_nums(int num_elements)
{
  data_type *rand_nums = (data_type *)malloc(sizeof(data_type) * num_elements);
  assert(rand_nums != NULL);
  int i;
  for (i = 0; i < num_elements; i++)
  {
    rand_nums[i] = (rand() / (data_type)RAND_MAX);
  }
  return rand_nums;
}

data_type *create_fixed_nums(int num_elements, int world_rank)
{
  data_type *rand_nums = (data_type *)malloc(sizeof(data_type) * num_elements);
  assert(rand_nums != NULL);
  int i;
  for (i = 0; i < num_elements; i++)
  {
    rand_nums[i] = world_rank + i * 0.1;
  }
  return rand_nums;
}

data_type *inilize_arr(int num_elements)
{
  data_type *rand_nums = (data_type *)malloc(sizeof(data_type) * num_elements);
  assert(rand_nums != NULL);
  memset(rand_nums, 0, sizeof(data_type) * num_elements);
  return rand_nums;
}

void *inilize_arr_withoutset(int num_elements)
{
  void *rand_nums = (void *)malloc(sizeof(data_type) * num_elements);
  assert(rand_nums != NULL);
  return rand_nums;
}

int get_pof2(int num)
{   
    int pof2 = 1;
    while(num > 1)
    {
        num = num/2;
        pof2 *= 2;
    }
    return pof2;
}