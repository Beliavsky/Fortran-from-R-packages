

bool binvecGE(unsigned long long *x, unsigned long long *y, int size)
{
  int i = size - 1;
  for(; i >= 0 and x[i] == y[i]; --i);
  return i < 0 or x[i] > y[i];
}


// bool binvecGE(unsigned long long *x, unsigned long long *y, int size)
// {
//   int i = size - 1;
//   for(; i >= 0; --i)
//   {
//     if(x[i] > y[i]) return true;
//     if(x[i] != y[i]) return false;
//   }
//   return true;
// }




